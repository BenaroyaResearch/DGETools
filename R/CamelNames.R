#
# camel Case operations
# Scott Presnell, SPresnell@benaroyaresearch.org
# January 28th, 2019
#
if(getRversion() >= "3.1.0") utils::globalVariables(c("."))
#' Convert column names to camelCase
#'
#' Convert column names to camelCase by taking them all to lower case, and looking for word boundaries and makeing those upper case.
#' Then trim off the white space.
#'
#' @param names list of names as from dataframe colnames().
#'
#' @return A list of names converted to camel case
#'
#' @author Scott R Presnell, \email{SPresnell@@benaroyaresearch.org}
#'
#' @export
#' @import stringr magrittr
#'
CamelNames <- function(names){
  newNames <- names %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all(pattern="-", replacement="_") %>%
    # the dot operator in a closure after %>% means use the piped material here.
    { gsub("^(.)", "\\L\\1", ., perl=T) } %>%
    { gsub(" (.)", "\\U\\1", ., perl=T) } %>%
    { gsub("_(.)", "\\U\\1", ., perl=T) } %>%
    # str_replace_all doesn't support PCRE escapes
    #    str_to_lower() %>% # change variable names to lower case
    #    str_replace_all(pattern = "[:punct:]", replacement="_") %>% #Remove punctuation marks
    #    str_replace_all(pattern = " ", replacement="_") %>% # remove whitespace
    #    str_replace_all(pattern = "__", replacement = "_") %>% #Remove double underscores
    #    str_replace_all(pattern = "_$", replacement = "") %>% #Remove underscores at the end of variables
    #    make.unique(sep = "_") %>% # de-dup variable names
    stringr::str_trim() #trim whitespace

  #  dframe <- as.data.frame(dframe) #Make a data frame as tibbles sometimes have unexpected downstream behavior

  return(newNames)
}

