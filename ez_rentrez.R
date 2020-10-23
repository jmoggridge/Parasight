## EZrentrez: functions to handle large queries to NCBI nuccore

require(rentrez)
require(tidyverse)
require(Biostrings)

# link <- rentrez::entrez_link(
#   dbfrom = 'nuccore', 
#   db = 'taxonomy', 
#   web_history = basic_search$web_history,
#   )
# link <- tibble(enframe(link))



#' @title get_ncbi_ids
#'
#' @param searchexp the search expression string with terms, booleans, etc.
#' @param db the ncbi database to search
#' @return the search results for the given expression with all available ids, and a webhistory token for entrez_summary or entrez_fetch
#' 
get_ncbi_ids <- function(searchexp, db){
  
  # find out how many ids available from 1st search
  Esearch <- rentrez::entrez_search(
    db = db, 
    term = searchexp, 
    retmax = 100)
  # use 'count' from first search to get all ids; get web history with 'use_history = TRUE'
  Esearch2 <- rentrez::entrez_search(
    db = db, 
    term = searchexp, 
    retmax = Esearch$count,
    use_history = TRUE)
  message('Returning metadata from basic search...')
  print(Esearch2)
  return(Esearch2)
}
# ## vignette:
# apikey <- 'b1a183199be1617b6de4f030ade08'
# searchexp <- '18S ribosomal rna[Title] AND "apicomplexa"[Organism] AND 300:2500[SLEN] AND biomol_genomic[PROP] NOT genome[TITL] AND (ddbj_embl_genbank[filter] OR refseq[filter])'
# basic_search <- get_ncbi_ids(searchexp, db = 'nuccore')
# basic_search$web_history
# basic_search$count
# 
# 
# searchexp <- 'apicomplexa'
# basic_search <- get_ncbi_ids(searchexp, db = 'taxonomy')



#' @title get_ESummary_df
#'
#' @details A wrapper for rentrez::entrez_search and rentrez::entrez_summary that takes your search expression, gets all the summaries, then takes the list of summaries and returns it as a tidy tibble. You'll need to provide your own api key from NCBI.
#' @param searchexp the search expression for ncbi
#' @param apikey your apikey from ncbi
#'
#' @return a tibble of all summary records available for given search expression
#'
#' @examples
#' apikey <- 'b1a183199be1617b6de4f030ade08' # use your own apikey
#' searchexp <- '18S AND apicomplexa[ORGN] AND 0:10000[SLEN] NOT (genome[TITL])'
#' summary.df <- get_ESummary_df(searchexp, apikey)

get_ESummary_df <- function(searchexp, db, apikey){
  
  # perform basic search to get webhistory and records count
  basic_search <- get_ncbi_ids(searchexp, db)
  web_history <- basic_search$web_history
  id.count <- basic_search$count
  
  # init list to gather downloads 
  Esummary_list <- list()
  # init df to compile records from all downloads
  df <- tibble(id = basic_search$ids)
  
  # display downloads progress
  message('Getting summaries....')
  progress_bar = txtProgressBar(min=0, max=id.count, style = 1, char="=")
  
  # iterate until all summary records obtained
  for (i in seq(1, id.count, 500)) {
    # print(paste0(round(i/id.count*100), '%'))
    # add each query to growing Esummary list
    Esummary_list <-  c(
      Esummary_list,
      entrez_summary(db = db, web_history = web_history,
                     retstart = i-1, retmax=500, api_key = apikey,
                     always_return_list = TRUE, retmode = 'json')
    )
    setTxtProgressBar(progress_bar, value = i)
    Sys.sleep(0.11)
  }
  
  message('\ndone downloads....\nflattening lists....')
  df <- df %>%
    mutate(listcol = Esummary_list) %>%
    unnest_wider(listcol)
  return(df)
}

# 
# # ncbi api key, obviously you'll need to put your own (don't copy this fake one)
# 
# apikey <- 'b15a151c3a183199be1617b6de4f030ade08' 
# searchexp <- '18S ribosomal rna[Title] AND "apicomplexa"[Organism] AND 300:2500[SLEN] AND biomol_genomic[PROP] NOT genome[TITL] AND (ddbj_embl_genbank[filter] OR refseq[filter])'
# 
# searchexp <- 'apicomplexa'
# summary.df <- get_ESummary_df(searchexp, apikey, db = 'taxonomy')
# # 


#' get_Efasta
#' 
#' gets all sequences from all search results in fasta format, using webhistory from get_ncbi_ids(searchexp)
#' @param searchexp the search expression that you want all fasta records for
#' @param apikey your ncbi apikey for rapid downloads
#'
#' @return returns a dataframe with columns for 'accession + title' and sequence
#' 
#'
#' @examples
#' # do initial search
#' apicomplexa.search <- get_ncbi_ids(searchexp)
#' # keep the web history link and records count
#' webhist <- apicomplexa.search$web_history
#' id.count <- apicomplexa.search$count
#' apicomplexa.fasta <- get_Efasta(webhist, id.count, apikey)
#' 
get_Efasta <- function(searchexp, apikey) {
  basic_search <- get_ncbi_ids(searchexp)
  webhist <- basic_search$web_history
  id.count <- basic_search$count
  
  message(paste('Fetching', id.count, 'fasta records'))
  progress_bar = txtProgressBar(min=0, max=id.count, style = 1, char="=")
  
  fasta <- ''
  for (i in seq(1, id.count, 500)){
    print(paste0(round(i/id.count*100), '%'))
    temp <- entrez_fetch(db = "nuccore", web_history = webhist,
                         retstart = i-1, retmax=500, api_key = apikey,
                         rettype = 'fasta')
    fasta <- paste0(fasta, temp)
    setTxtProgressBar(progress_bar, value = i)
  }
  print('100%')
  message('Victory!\nTidying data....')
  # split fasta into lines and write a temp file
  write(fasta, "temp.fasta", sep = "\n")
  # read lines as a DNAStringSet object using Biostrings package
  fasta <- Biostrings::readDNAStringSet("temp.fasta")
  # delete temp file
  file.remove('temp.fasta')
  # arrange fasta into df
  fasta.df <- data.frame(title = names(fasta), seq = paste(fasta))
  message('Done.')
  return(fasta.df)
}

