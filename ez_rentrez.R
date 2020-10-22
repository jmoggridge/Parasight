## EZrentrez: functions to handle large queries to NCBI nuccore

require(rentrez)
require(tidyverse)
require(Biostrings)

#' @title get_ncbi_ids
#'
#' @param searchexp the search expression string with terms, booleans, etc.
#'
#' @return the search results for the given expression with all available ids, and a webhistory token for entrez_summary or entrez_fetch
get_ncbi_ids <- function(searchexp){
  
  # find out how many ids available from 1st search
  Esearch <- entrez_search(db = "nuccore", term = searchexp, retmax = 100)
  # use 'count' from first search to get all ids; get webenv with 'use_history = TRUE'
  Esearch2 <- entrez_search(db = "nuccore", 
                            term = searchexp, retmax = Esearch$count,
                            use_history = TRUE)
  message('Returning ids from Entrez nuccore search:')
  print(Esearch2)
  return(Esearch2)
}


#' @title get_Esummaries
#'
#' @details subroutine run by get_Esummary_df(); wrapper for entrez_summary() to get large volumes of summary records for a given search expression in 500 records/query steps.
#' 
#' @param web_history a web history token returned by entrez_search
#' @param id.count the total number of available records from entrez_search
#' @param apikey an apikey is necessary to do multiple downloads/second
#'
#' @return a list of summary records. Same as output from entrez_summary but with multiple queries' results concatenated into a larger list
#' 
#' @examples Esummary_list <- get_Esummaries(webhist, id.count, apikey)
#' 
get_Esummaries <- function(web_history, id.count, apikey) {
  # init list to gather downloads
  Esummary_list <- list()
  # iterate until all summary records obtained
  for (i in seq(1, id.count, 500)) {
    # display downloads progress
    print(paste0(round(i/id.count*100), '%'))
    # add each query to growing Esummary list
    Esummary_list <-  c(
      Esummary_list,
      entrez_summary(db = "nuccore", web_history = web_history,
                     retstart = i-1, retmax=500, api_key = apikey,
                     always_return_list = TRUE, retmode = 'json')
    )
    Sys.sleep(0.11)
  }
  return(Esummary_list)
}

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
#' searchexp <- '18S AND apicomplexa[ORGN] AND 0:10000[SLEN]) NOT (genome[TITL])'
#' summary.df <- get_ESummary_df(searchexp, apikey)

get_ESummary_df <- function(searchexp, apikey){
  basic_search <- get_ncbi_ids(searchexp)
  webhist <- basic_search$web_history
  id.count <- basic_search$count
  df <- tibble(id = basic_search$ids)
  print(glimpse(df))
  message('Getting summaries....')
  Esummary_list <- get_Esummaries(webhist, id.count, apikey)
  message('done downloads....\nflattening lists....')
  df <- df %>%
    mutate(listcol = Esummary_list) %>%
    unnest_wider(listcol)
  return(df)
  }

# ncbi api key, obviously you'll need to put your own (don't copy this fake one)
# apikey <- 'b1a183199be1617b6de4f030ade08' 
# searchexp <- '18S AND apicomplexa[ORGN] AND 0:10000[SLEN]) NOT (genome[TITL])'
# summary.df <- get_ESummary_df(searchexp, apikey)
# 



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
  fasta <- ''
  for (i in seq(1, id.count, 500)){
    print(paste0(round(i/id.count*100), '%'))
    temp <- entrez_fetch(db = "nuccore", web_history = webhist,
                         retstart = i-1, retmax=500, api_key = apikey,
                         rettype = 'fasta')
    fasta <- paste0(fasta, temp)
  }
  print('100%')
  message('Victory! You downloaded those seqs like a champ :)\n\nTidying data....')
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

