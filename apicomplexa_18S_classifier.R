library(rentrez)
library(seqinr)
library(Biostrings)
library(tidyverse)
# library(ggbeeswarm)
library(patchwork)
library(ggthemes)
library(ggridges)
library(ggalluvial)
library(ggrepel)
library(lubridate)

library(randomForest)
#naivebayes

taxid <- entrez_search(db="taxonomy", term="apicomplexa")

# Downloads ----

# my own 'ez_rentrez' suite of functions that wrap rentrez:: functions
source('./ez_rentrez.R')

# ncbi api key
apikey <- 'b15a151c3a183199be1617b6de4f030ade08' 
# apicomplexa 18S, not longer than 10kbp, and not from a genome
searchexp <- '18S ribosomal rna[Title] AND "apicomplexa"[Organism] AND 300:2500[SLEN] AND biomol_genomic[PROP] NOT genome[TITL] AND (ddbj_embl_genbank[filter] OR refseq[filter])'

# # Get search metadata, summaries, and seqs in fasta format
# apicomplexa.search <- get_ncbi_ids(searchexp)
# apicomplexa.summary.df <- get_ESummary_df(searchexp, apikey)
# apicomplexa.fasta <- get_Efasta(searchexp, apikey)
#  
# # save downloaded data
# write_rds(apicomplexa.search, './Assignment2/data/apicomplexa_search_ncbi.rds')
# write_rds(apicomplexa.summary.df, './Assignment2/data/apicomplexa_summaries_ncbi.rds')
# write_rds(apicomplexa.fasta, './Assignment2/data/apicomplexa_fasta.rds')

rm(list = ls())

#-----

### DATASET ----

# read downloaded data
# search metadata
apicomplexa.search <- read_rds('./Assignment2/data/apicomplexa_search_ncbi.rds')
# summary data
apicomplexa.summary.df <- read_rds('./Assignment2/data/apicomplexa_summaries_ncbi.rds')
# sequence data
apicomplexa.fasta <- read_rds('./Assignment2/data/apicomplexa_fasta.rds')
# rows from summary + sequence are in same order match (accession.version + title) from .fasta to title in .summary.df

# make full sequences + summary dataframe
apicomplexa.df <- apicomplexa.summary.df %>%
  # add sequence data to summary.df
  bind_cols(apicomplexa.fasta) %>%
  # parse dates
  mutate(across(contains('date'), function(x) lubridate::ymd(x))) %>%
  # replace blanks with explicit NA
  mutate(across(where(is.character), ~na_if(.x, ''))) %>%
  # make strings into factors
  mutate(across(where(is.character), as.factor)) %>%
  mutate(taxid = factor(taxid)) %>%
  rename(title = title...4)  %>%
  # drop cols that are all NAs
  janitor::remove_empty('cols')

# just checking that dataframes were matched up in correct order
x <- apicomplexa.df %>%
  mutate(acc1 = str_remove(accessionversion, '\\.[0-9]+'),
         matched = ifelse(acc1 == caption, TRUE, FALSE)) %>%
  select(uid, caption, acc1, matched, accessionversion, everything()) %>%
  filter(matched==FALSE)
nrow(x)==0

glimpse(apicomplexa.df)
rm(apicomplexa.fasta, apicomplexa.search, apicomplexa.summary.df, x)

# What columns have useful data?
summary(apicomplexa.df %>%
          select(slen, biomol, geneticcode, topology, strand,
                 sourcedb, completeness, createdate, genome,
                 organism, taxid))

### FILTERING -----


apicomplexa.select.df <- apicomplexa.df %>%
  # not rna, not mitochondrial or ambiguous
  filter(!moltype == 'rna' & genome == 'genomic') %>%
  # not uncultured or ambiguous species name
  filter(!str_detect(organism, 'sp\\.$|uncultured')) %>%
  # keep only useful columns
  select(uid, caption, title, organism, taxid, slen, seq)

summary(apicomplexa.select.df %>% select(-seq, -title))  

# fix organism names & extract genera: 
apicomplexa.select.df <- apicomplexa.select.df %>%
  mutate(
    # 'malaria parasite P. <species>' -> Plasmodium species
    organism = str_replace(organism, 'malaria parasite P\\.', 'Plasmodium'),
    # extract genus names
    genus = str_extract(organism, '^[A-Za-z]+'))

summary(factor(apicomplexa.select.df$organism))
summary(factor(apicomplexa.select.df$genus))

# create genera dataset with >=20 representative sequences (65 gen.)
genera <- names(which(table(apicomplexa.select.df$genus)>20))
api.genera.df <- apicomplexa.select.df %>%
  filter(genus %in% genera) %>%
  filter(!genus =='Apicomplexa')

# create species dataset from taxids with >=20 rep seqs
species <- names(which(table(apicomplexa.select.df$organism)>20))
api.species.df <- apicomplexa.select.df %>%
  filter(organism %in% species)

rm(species, genera, apicomplexa.select.df, apicomplexa.df)

# EDA ----


# histogram showing distribution of sequences per taxon id
a <- api.species.df %>%
  count(organism) %>%
  ggplot(aes(x=n)) +
  geom_histogram(bins = 20) +
  scale_x_log10() +
  labs(subtitle = 'species dataset',
       x = 'representative sequences (n)', 
       y = 'taxa (by organism name)')  +
  theme_tufte(base_family = 'sans')
b <- api.genera.df %>%
  count(genus) %>%
  ggplot(aes(x=n)) +
  geom_histogram() +
  scale_x_log10() +
  labs(subtitle = 'genera dataset',
       x = 'representative sequences (n)', 
       y = 'taxa (by genus name)') +
  theme_tufte(base_family = 'sans')

a+b & plot_annotation(
  caption = 'Figure 1. Histograms showing the distribution of representative\n sequences per taxa at the species (left) and genus-rank (right)',
  theme = theme(plot.caption = element_text(hjust=0, size=rel(1))))

rm(a,b)

### FEATURE GENERATION ----

## Tidy up sequences data

# remove leading NNNs
api.genera.df <- api.genera.df %>%
  mutate(seq = str_remove(seq, "^[-N]+")) %>%
  mutate(seq = str_remove(seq, "[-N]+$")) %>%
  # mutate(seq = str_remove_all(seq, "-+")) %>%
  filter(str_count(seq, "N") <= (0.05 * str_count(seq)))

# plot sequence length density by genus
 api.genera.df %>% 
  ggplot(aes(x=slen, y = genus)) +
  geom_density_ridges(panel_scaling = TRUE, alpha = 0.7) +
  labs(x = '18S sequence length (bp)', y = 'Genus',
       caption = '\nFigure 2. 18S gene sequence length varies by genera (after trimming ambiguous bases') +
  scale_y_discrete(limits = rev(levels(api.genera.df$genus))) +
  theme_tufte(base_family = 'sans') +
  theme(plot.caption = element_text(hjust=0, size=rel(1)))

 
class(api.genera.df)

## Count kmers
api.genera.df <- as.data.frame(api.genera.df)
api.genera.df$seq <- DNAStringSet(api.genera.df$seq)

api.genera.features.df <- api.genera.df %>%
  # single nucleotide frequency
  cbind(
    letterFrequency(api.genera.df$seq, 
                        letters = c('A', 'C','G'),
                        as.prob = TRUE),
    letterFrequency(api.genera.df$seq, 
                    letters = c('CG'),
                    as.prob = TRUE),
    dinucleotideFrequency(api.genera.df$seq, 
                          as.prob = TRUE),
    trinucleotideFrequency(api.genera.df$seq,
                           as.prob = TRUE)
    # oligonucleotideFrequency(api.genera.df$seq, 
    #                          width = 4, 
    #                          as.prob = TRUE)
    )
  

# Splitting train/test data


api.genera.features.df <-  api.genera.features.df %>%
  mutate(seq = as.character(seq), 
         genus = as.factor(genus))
class(api.genera.df)
summary(factor(api.genera.df$genus))


train <- api.genera.features.df %>%
  select(-seq)
  # dplyr::filter(genus %in% c('Cyclospora', 'Goussia', 'Toxoplasma'))

summary(factor(train$genus))
glimpse(train)

test <- api.genera.features.df %>%
  group_by(genus) %>%
  mutate(n=n()) %>%
  sample_n(floor(n/4)) %>%
  select(uid, genus, 9:92)
glimpse(test)
summary(factor(test$genus))

train <- api.genera.features.df %>%
  filter(!uid %in% test$uid) %>%
  select(genus, 9:92)
summary(factor(train$genus))

test$uid <- NULL
genus.table <- train %>%
  count(genus) %>%
  select(genus, n) %>% rename(count=n)

## RANDOM FOREST -----
set.seed(1)
randomforest.genus <- 
  randomForest(x = train[, 2:85],
               y = train$genus,
               ntree = 500, importance = TRUE)

randomforest.genus

rf.confusion <- as.data.frame(randomforest.genus$confusion)
rf.confusion$actual <- rownames(randomforest.genus$confusion)

ggplot(rf.confusion) +
  geom_point(aes(y=actual, x = class.error))
glimpse(rf.confusion)

genus.table <- genus.table %>%
  bind_cols(error = rf.confusion$class.error)
genus.table



rf.confusion %>%
  select(-class.error) %>%
  select(actual, everything())
#We can also dive into the specifics, such as by looking at the following. See the documentation for the randomForest() function for more information.
rf.importance <- data.frame(gene_classifier$importance)
gene_classifier$oob.times

summary(gene_classifier)
gene_classifier$classes
## Naive bayes classifier




# create classifier
m <- naiveBayes(genus ~ ., data = iris)
## alternatively:
m <- naiveBayes(iris[,-5], iris[,5])
m
table(predict(m, iris), iris[,5])
# }



# create classifier
library(e1071)
m <- naiveBayes(train[,-1], train[,1])
m
predicted <- predict(m, test)
# }
table(predicted, train[,1])
glimpse(predicted)






### plots ----

# plot log.error ~ log.count
ggplot(genus.table, aes(y=log(error), x = log(count))) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_text_repel(aes(label = genus)) +
  # ylim(0, .7) +
  geom_smooth(method = lm, formula = y~x+I(x^2) + I(x^3), se=FALSE) + 
  labs(caption = 'The relationship between classification error rate and training sequence count') +
  theme_tufte(base_family = 'sans')