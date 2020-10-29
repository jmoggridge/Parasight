
# Kmer-based classification of Apicomplexa 18S rRNA gene #

#libraries
library(tidyverse)
library(rentrez)
library(janitor)
library(seqinr)
library(Biostrings)
library(patchwork)
library(ggthemes)
library(ggridges)
library(caret)
library(klaR) # naive bayes classifier
library(randomForest)
library(doParallel)
library(beepr)
rm(list=ls())



## Downloads ----
# 
# # suite of functions that wrap rentrez::entrez_* functions
# source('./ez_rentrez.R')
# 
# # ncbi api key
# apikey <- '889fdb9786a14019a0a1257196a09ba4ba08' 
# # apicomplexa 18S, not longer than 10kbp, and not from a genome
# searchexp <- '18S ribosomal rna[Title] AND "apicomplexa"[Organism] AND 300:2500[SLEN] AND biomol_genomic[PROP] NOT genome[TITL] AND (ddbj_embl_genbank[filter] OR refseq[filter])'
# 
# # Get search metadata, summaries, sequences, and taxonomy info
# apicomplexa.search <- get_ncbi_ids(searchexp, db = 'nuccore')
# apicomplexa.summary.df <- get_ESummary_df(searchexp, 'nuccore', apikey)
# apicomplexa.fasta <- get_Efasta(searchexp, apikey)
# apicomplexa.taxa <- get_nuccore_taxonomy(apicomplexa.search$web_history)
# 
# # save downloaded data
# write_rds(apicomplexa.search, './data/apicomplexa_search_ncbi.rds')
# write_rds(apicomplexa.summary.df, './data/apicomplexa_summaries_ncbi.rds')
# write_rds(apicomplexa.fasta, './data/apicomplexa_fasta.rds')
# write_rds(apicomplexa.taxa, './data/apicomplexa_taxonomy.rds')
# # read data
# apicomplexa.search <- read_rds('./data/apicomplexa_search_ncbi.rds')
# apicomplexa.summary.df <- read_rds('./data/apicomplexa_summaries_ncbi.rds')
# apicomplexa.fasta <- read_rds('./data/apicomplexa_fasta.rds')
# apicomplexa.taxa <- read_rds('./data/apicomplexa_taxonomy.rds')
# 
# # link all dataframes
# apicomplexa.df <- apicomplexa.summary.df %>%
#   # add sequence data to summary.df
#   bind_cols(apicomplexa.fasta) %>%
#   # link taxonomy information using 'taxid' as key
#   mutate(taxid = as.character(taxid)) %>%
#   left_join(., apicomplexa.taxa, by = 'taxid')
# 
# write_rds(apicomplexa.df, './data/apicomplexa.rds')
# rm(list=ls())

#-----

## Tidying data ----

# read downloaded data
apicomplexa.df <- read_rds('./data/apicomplexa.rds')


# make full sequences + summary dataframe
apicomplexa.df <- apicomplexa.df %>%
  # make strings into factors; rename title
  mutate(across(where(is.character), as.factor)) %>%
  mutate(taxid = factor(taxid)) %>%
  dplyr::rename(title = title...4)  %>%
  # replace blanks with explicit NA
  mutate(across(where(is.character), ~na_if(.x, ''))) %>%
  # drop any cols that are all NAs; drop unwanted cols
  janitor::remove_empty('cols') %>%
  dplyr::select(
    -c(biosample, strain, statistics, geneticcode, tech,
       contains('sub'), projectid, properties, oslt, title...34,
       extra, contains('assembly')))

# just checking that dataframes were matched up in correct order
x <- apicomplexa.df %>%
  mutate(acc1 = str_remove(accessionversion, '\\.[0-9]+'),
         matched = ifelse(acc1 == caption, TRUE, FALSE)) %>%
  dplyr::select(uid, caption, acc1, matched, accessionversion,
                everything()) %>%
  filter(matched==FALSE)
nrow(x)==0

glimpse(apicomplexa.df)
rm(x)



# Tidying sequences, taxonomy info
apicomplexa.df <- apicomplexa.df %>%
  # trim lineages, fix organism names
  mutate(
    Lineage = str_remove(
      Lineage,'cellular organisms; Eukaryota; Sar; Alveolata; Apicomplexa; '),
    organism = str_replace(organism, 'malaria parasite P\\.', 'Plasmodium'),
    organism = ifelse(organism =='Babesia canis canis', 'Babesia canis', organism),
    # extract genus names from organism name
    genus = str_extract(organism, '^[A-Za-z]+'),
  # Tidy sequences by removing head/tail Ns and any gaps
    seq = str_remove(seq, regex("^[-N]+|[-N]+$")),
    seq = str_remove_all(seq, "-+")) 


## Filtering -----

apicomplexa.filter.df <- apicomplexa.df %>%
  # not rna, not mitochondrial or ambiguous
  filter(!moltype == 'rna' & genome == 'genomic') %>%
  # not uncultured prganisms or ambiguous species names
  filter(!str_detect(organism, 'sp\\.$|uncultured'))  %>%
  # remove sequences with more than 0.5% internal N
  filter(str_count(seq, "N") <= (0.005 * str_count(seq))) %>%
  # keep only useful columns
  dplyr::select(id, organism, genus, seq, title, createdate,
                taxid, slen, ScientificName, ParentTaxId, 
                Rank, Lineage)

summary(apicomplexa.filter.df %>% dplyr::select(-seq, -title))  

## remove any duplicated sequences, keep first instance
apicomplexa.filter.df <- apicomplexa.filter.df %>%
  group_by(seq) %>%
  sample_n(1) %>%
  ungroup()


# create genera dataset with >=300 representative sequences (65 genera)
table(apicomplexa.filter.df$genus) %>% sort(decreasing = TRUE)
genera <- names(which(table(apicomplexa.filter.df$genus)>300))
api.genera.df <- apicomplexa.filter.df %>%
  filter(genus %in% genera) %>%
  filter(!genus =='Apicomplexa') %>%
  mutate(genus = as.factor(genus))

table(api.genera.df$genus)

## down-sampling all genera to match genus with smallest number of sequences
api.genera.df <- downSample(api.genera.df, factor(api.genera.df$genus))
table(api.genera.df$genus)

rm(genera, apicomplexa.filter.df, apicomplexa.df)


# EDA ----

# histogram showing distribution of sequences per taxon id
# a <- api.species.df %>%
#   count(organism) %>%
#   ggplot(aes(x=n)) +
#   geom_histogram(bins = 20) +
#   scale_x_log10() +
#   labs(subtitle = 'species dataset',
#        x = 'representative sequences (n)', 
#        y = 'taxa (by organism name)')  +
#   theme_tufte(base_family = 'sans')

# b <- api.genera.df %>%
#   dplyr::count(genus) %>%
#   ggplot(aes(x=n)) +
#   geom_histogram() +
#   scale_x_log10() +
#   labs(subtitle = '',
#        caption = 'Figure 1. Histogram showing the distribution of genbank sequences\n in apicomplexan genera',
#        x = 'Representative sequences (n)', 
#        y = 'Genera') +
#   theme_tufte(base_family = 'sans') +
#   geom_rangeframe() +
#   theme(plot.caption = element_text(hjust=0, size=rel(0.9)))

# plot sequence length density by genus
api.genera.df %>% 
  ggplot(aes(x=slen, y = genus)) +
  geom_density_ridges(panel_scaling = TRUE, alpha = 0.7) +
  labs(x = '18S sequence length (bp)', y = 'Genus',
       caption = 'Figure 2. 18S gene sequence length varies by genera (after trimming ambiguous bases') +
  scale_y_discrete(limits = rev(levels(api.genera.df$genus))) +
  theme_tufte(base_family = 'sans') +
  theme(plot.caption = element_text(hjust=1.4, size=rel(0.8)))


### FEATURE GENERATION ----

# generates sequence feature columns for classification
generate_kmer_features <- function(df){
  df <- as.data.frame(df)
  df$seq <- DNAStringSet(df$seq)
  features.df <- df %>%
    #  (1-4)-mer frequency
    cbind(
      letterFrequency(df$seq, letters = c('A', 'C','G'),
                      as.prob = TRUE),
      dinucleotideFrequency(df$seq, as.prob = TRUE),
      trinucleotideFrequency(df$seq, as.prob = TRUE),
      oligonucleotideFrequency(df$seq, 4, as.prob = TRUE)
      # ,
      # oligonucleotideFrequency(df$seq, 5, as.prob = TRUE)
    )
  return(features.df)
}

api.genera.df <-  api.genera.df %>%
  generate_kmer_features(.) %>%
  dplyr::select(-c(seq, Class))

glimpse(api.genera.df)
write_rds(api.genera.df, './data/genus_rank data for classification.rds')


## Splice data into train + test ----

set.seed(1)
genus.data <- api.genera.df
# create a stratified partition of data: 50/50 split train:test
inTrainingSet <- createDataPartition(y = genus.data$genus, p=.5, list=FALSE)

# split data rows into train, test
genusTrain <- genus.data[inTrainingSet, ] %>%
  dplyr::select(genus, 13:length(genus.data))
genusTest <- genus.data[-inTrainingSet, ] %>%
  dplyr::select(genus, 13:length(genus.data))

# we can see that sampling prob is even across genera
table(genusTrain$genus)
table(genusTest$genus)
rm(inTrainingSet, genus.data)

### CLASSIFIERS ----

# column index for all k-mer proportions
features <- setdiff(names(genusTrain), 'genus')


## Naive Bayes 1----

# specify bootstrap with reps for better estimation of model accuracy
ctrl <- trainControl(
  method = "boot", # boot, cv, repeatedcv,...
  number = 5,
  classProbs = TRUE,
  allowParallel = TRUE
)

# set tuning parameters search grid:
# parameter   class                label
# 1        fL numeric   Laplace Correction
# 2 usekernel logical    Distribution Type
# 3    adjust numeric Bandwidth Adjustment
search_grid <- expand.grid(
  fL = 0:2,
  usekernel = c(TRUE, FALSE), # false doesn't work
  adjust = seq(0,4,2) # 0 doesnt work
)

# setup cluster of processors (4 on my macbook)
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)

# keep track of runtime, for reference
timed <- Sys.time()
#
nb.model <- train(x = genusTrain[, features],
                  y = genusTrain$genus,
                  method ='nb',
                  trControl = ctrl,
                  tuneGrid = search_grid,
                  verbose = FALSE
                  )
# tell me when it finishes, if ever
print('done')
beep(8)
timed <- Sys.time() - timed
stopCluster(cl)
write_rds(nb.model, './models/naive_bayes.model.rds')

nb.model <- read_rds('./models/naive_bayes.model.rds')
nb.model
nb.model$results
nb.model$bestTune


ggplot(nb.model$results,
       aes(x = adjust)) +
  geom_linerange(aes(ymax=Accuracy + AccuracySD,
                    ymin = Accuracy - AccuracySD),
                stat='identity', colour = 'red4', alpha = 0.6) +
  geom_line(aes(y = Accuracy), colour = 'red4', alpha = 0.6) +
  geom_point(aes(y = Accuracy), colour = 'red4', alpha = 0.6) +
  geom_linerange(aes(ymax = Kappa + KappaSD,
                    ymin = Kappa - KappaSD), colour = 'blue3', 
                 alpha = 0.6) +
  geom_line(aes(y = Kappa), color = 'blue3', alpha = 0.6) +
  geom_point(aes(y = Kappa), color = 'blue3', alpha = 0.6) +
  ylim(0.6,0.9)

caret::plot.train(nb.model)
plot(nb.model, metric = 'Kappa')

test.nb.model <- predict(nb.model$finalModel, genusTest[, features])
nb1.confusion <- confusionMatrix(table(test.nb.model$class, genusTest$genus))
beep(7)
rm(ctrl, search_grid)


## Naive Bayes 2: try to optimize parameters more (by doing less initially) ----

# using default search from caret::train

# setup cluster of processors (4 on my macbook)
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)

# keep track of runtime, for reference
timed <- Sys.time()
features <- setdiff(names(genusTrain), 'genus')
#
nb.model2 <- train(x = genusTrain[, features],
                   y = genusTrain$genus,
                   method ='nb')
# tell me when it finishes, if ever
print('done')
beep(8)
timed <- Sys.time() - timed
stopCluster(cl)

write_rds(nb.model2, './models/naive_bayes.model2.rds')
nb.model2 <- read_rds('./models/naive_bayes.model2.rds')

nb.model2
nb.model2$results
nb.model2$finalModel$tuneValue

# Test classifier on test data
test.nb.model2 <- predict(nb.model2$finalModel, genusTest[, features])
nb2.confusion <- confusionMatrix(table(test.nb.model2$class, genusTest$genus))
beep(1)


# Naives Bayes 3: with preprocessing  (center, scale) ----

# setup cluster of processors (4 on my macbook)
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)

# keep track of runtime, for reference
timed <- Sys.time()

# train the naive bayes model with preprocessing of kmer proportions
nb.model3 <- caret::train(
  x = genusTrain[, features],
  y = genusTrain$genus,
  method ='nb',
  preProc= c("center", "scale")
)

# tell me when it finishes, if ever
print('done')
beep(7)
timed <- Sys.time() - timed
stopCluster(cl)
write_rds(nb.model3, './models/naive_bayes.model3.rds')

nb.model3 <- read_rds('./models/naive_bayes.model3.rds')
nb.model3
nb.model3$results


# Test classifier on test data: 
# don't specify final model, s.t. test data are preprocessed in the same way as the training data were in the nb.model3
test.nb.model3 <- predict(nb.model3, 
                          newdata = genusTest[, features])
summary(test.nb.model3)
nb3.confusion <- confusionMatrix(table(test.nb.model3, genusTest$genus))
nb3.confusion
beep(1)

rm(ctrl, timed, grid, cl)

## Random Forest 1-------

# specify bootstrap with reps for better estimation of model accuracy
ctrl <- caret::trainControl(
  method = "boot", # boot, cv, repeatedcv,...
  number = 5,
  classProbs = TRUE,
  allowParallel = TRUE
)

# set tuning parameters search grid 
# mtry -> number of vars sampled at each split
search_grid <- expand.grid(.mtry = 8:13)

# setup cluster of processors (4 on my macbook)
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)

# keep track of runtime, for reference
timed <- Sys.time()

# train model with search grid and resample procedure
rf.model <- caret::train(x = genusTrain[, features],
                  y = genusTrain$genus,
                  method ='rf', 
                  ntree = 1000,
                  trControl = ctrl,
                  tuneGrid = search_grid,
                  verbose = FALSE
)
# tell me when it finishes, if ever
print('done')
beepr::beep(8)
timed <- Sys.time() - timed
timed
stopCluster(cl)
write_rds(rf.model, './models/random_forest.model.rds')
rf.model <- read_rds('./models/random_forest.model.rds')

# check out the training results
rf.model
rf.model$results
plot(rf.model)
rf.model$finalModel # OOB estimate of  error rate: 1.9%
plot(rf.model$finalModel)
rf.model$finalModel$confusion

# Classify test dataset
test.rf.model <- predict(rf.model$finalModel, 
                         genusTest[, features])
rf.confusion <- confusionMatrix(table(test.rf.model, genusTest$genus))
rm(ctrl, search_grid)



# 
# ## Neural net classifier ----
# 
# # setup cluster of processors (4 on my macbook)
# cl <- parallel::makeCluster(4, setup_strategy = "sequential")
# registerDoParallel(cl)
# # keep track of runtime, for reference
# timed <- Sys.time()
# # train the naive bayes model with preprocessing of kmer proportions
# nnet.model <- caret::train(
#   x = genusTrain[, features],
#   y = genusTrain$genus,
#   method ='nnet',
#   preProc= c("center", "scale")
# )
# # tell me when it finishes, if ever
# print('done')
# beep(7)
# timed <- Sys.time() - timed
# stopCluster(cl)
# write_rds(nnet.model, './models/nnet.model.rds')
# nnet.model <- read_rds('./models/nnet.model.rds')
# 
# # wow, this model is terrible hahaha!
# nnet.model
# nnet.model$results
# ggplot(nnet.model)



## gbm - Stochastic Gradient Boosting ----

# setup cluster of processors (4 on my macbook)
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)
# keep track of runtime, for reference
timed <- Sys.time()
# train the naive bayes model with preprocessing of kmer proportions
gbm.model <- caret::train(
  x = genusTrain[, features],
  y = genusTrain$genus,
  method ='gbm',
  preProc= c("center", "scale")
)
# tell me when it finishes, if ever
print('done')
beep(7)
timed <- Sys.time() - timed
stopCluster(cl)
write_rds(gbm.model, './models/gbm.model1.rds')


gbm.model <- read_rds('./models/gbm.model1.rds')
ggplot(gbm.model)
gbm.model
gbm.model$results

# Classify test dataset
test.gbm.model <- predict(gbm.model, genusTest[, features])
gbm1.confusion <- confusionMatrix(table(test.gbm.model, genusTest$genus))
beep(1)
gbm1.confusion



## gbm2 ----
ctrl <- trainControl(
  method = "boot", # boot, cv, repeatedcv,...
  number = 10,
  classProbs = TRUE,
  allowParallel = TRUE
)
grid <- expand.grid(
  n.trees = seq(100, 1000, by = 150),
  interaction.depth = seq(1,7, by=2),
  shrinkage = c(.01, .1),
  n.minobsinnode = 10
)

# setup cluster of processors (4 on my macbook)
cl <- parallel::makeCluster(4, setup_strategy = "sequential")
registerDoParallel(cl)
# keep track of runtime, for reference
timed <- Sys.time()

# train the gbm model without preprocessing of kmer proportions
gbm.model2 <- caret::train(
  x = genusTrain[, features],
  y = genusTrain$genus,
  method ='gbm', 
  trControl = ctrl,
  tuneGrid = grid
)
d# tell me when it finishes, if ever
print('done')
beep(7)
Sys.time() - timed
stopCluster(cl)
write_rds(gbm.model2, './models/gbm.model2.rds')

gbm.model2 <- read_rds('./models/gbm.model2.rds')
ggplot(gbm.model2)
gbm.model2
gbm.model2$results
gbm.model2$finalModel


# Classify test dataset
test.gbm.model2 <- predict(gbm.model2, genusTest[, features])
gbm2.confusion <- confusionMatrix(table(test.gbm.model2,
                                        genusTest$genus))
beep(1)
gbm2.confusion



## Plots -----
ggplot(rf.model)
ggplot(nb.model)
ggplot(nb.model2)
ggplot(nb.model3)


## Compare Models ----

# 
# nb.confusion <- caret::confusionMatrix(table(test.nb.model$class, genusTest$genus))
# rf.confusion <- caret::confusionMatrix(table(test.rf.model, genusTest$genus))

confusion_overall <- tibble(
  name = names(nb.confusion$overall), 
  `GBM2` = gbm2.confusion$overall,
  `GBM1` = gbm1.confusion$overall,
  `Naive Bayes3` = nb3.confusion$overall,
  `Naive Bayes2` = nb2.confusion$overall,
  `Naive Bayes1` = nb1.confusion$overall,
  `Random Forest` = rf.confusion$overall,
  ) %>%
  # transpose table
  pivot_longer(GBM2:`Random Forest`, names_to = 'Classifier', values_to = 'value') %>%
  pivot_wider()
confusion_overall

ggplot(confusion_overall,
       aes(y = Classifier, x = Accuracy, xmin = AccuracyLower, xmax = AccuracyUpper)) +
  geom_pointrange() +
  scale_y_discrete(limits = rev(levels(factor(confusion_overall$Classifier)))) +
  labs()



# classConfusion_to_DF


nb <- as.data.frame(nb.confusion$byClass) %>%
  mutate(Genus = rownames(nb.confusion$byClass),
         Model='Naive Bayes',
         Genus = str_remove(Genus, 'Class: ')) %>%
  group_by(Genus, Model) %>%
  pivot_longer(where(is.numeric))

rf <- as.data.frame(rf.confusion$byClass) %>%
  mutate(Genus = rownames(rf.confusion$byClass),
         Model='Random Forest',
         Genus = str_remove(Genus, 'Class: ')) %>%
  group_by(Genus, Model) %>%
  pivot_longer(where(is.numeric))

confusion_class <- bind_rows(nb, rf) 
confusion_class <- confusion_class %>%
  filter(!str_detect(name, 'Prevalence|Rate|Pred|F1')) 
# filter(!name %in% c('Prevalence', 'F1', 'Neg Pred Value', 'Pos Pred Value'))



ggplot(confusion_class, aes(y=Genus, x = value, colour = Model)) +
  geom_point(alpha = 0.7) +
  geom_rangeframe(color = 'darkgrey') +
  scale_x_continuous(n.breaks = 3, ) +
  facet_wrap(name~., scales = 'free_x',ncol = 3) +
  labs(x = '') +
  theme_tufte(base_family = 'sans') +
  theme(legend.position = 'top', 
        legend.title = element_blank(),
        panel.grid.major.y = element_line(size = 0.2,
                                          colour = 'lightgray'),
        axis.text.x = element_text(size = 8))


rm(nb, rf, nb.confusion, rf.confusion)





