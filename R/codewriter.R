################################################################################
## Write out options to template Rmd file
################################################################################

analysisAPTS <- function(datafile, exposure = NULL, response = NULL,
                         timevar = NULL) {
        con <- file("input_parameters.Rmd", "w")
        on.exit(close(con))
        op <- options()
        options(useFancyQuotes = FALSE)
        on.exit(options(op), add = TRUE)
        writeLines("```{r input parameters}", con)
        writeLines(sprintf("datafile <- %s", dQuote(datafile)), con)
        writeLines(sprintf("exposure <- %s", ifelse(is.null(exposure), "NULL", dQuote(exposure))), con)
        writeLines(sprintf("response <- %s", ifelse(is.null(response), "NULL", dQuote(response))), con)
        writeLines(sprintf("timevar <- %s", ifelse(is.null(timevar), "NULL", dQuote(timevar))), con)
        writeLines("```", con)
}

writeInputParameters <- function(paramlist) {

}
