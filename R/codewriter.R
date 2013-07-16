################################################################################
## Write out options to template Rmd file
################################################################################

makeList <- function(...) {
        args <- substitute(list(...))
        nms <- sapply(args[-1], deparse)
        vals <- list(...)
        names(vals) <- nms
        vals
}

analysisAPTS <- function(datafile, exposure = NULL, response = NULL,
                         timevar = NULL) {
        writeInputParameters(datafile, exposure, response, timevar)
}

writeInputParameters <- function(...) {
        vals <- makeList(...)
        parnames <- names(vals)

        con <- file("input_parameters.Rmd", "w")
        on.exit(close(con))
        op <- options()
        options(useFancyQuotes = FALSE)
        on.exit(options(op), add = TRUE)
        writeLines("```{r input parameters}", con)
        for(i in seq_along(vals)) {
                if(is.null(vals[[i]]))
                        msg <- sprintf("%s <- NULL", parnames[i])
                else if(is.character(vals[[i]]))
                        msg <- sprintf("%s <- \"%s\"", parnames[i], vals[[i]])
                else if(is.numeric(vals[[i]]))
                        msg <- sprintf("%s <- %f", parnames[i], vals[[i]])
                else 
                        msg <- sprintf("%s <- %s", parnames[i],
                                       as.character(vals[[i]]))
                writeLines(msg, con)
        }
        writeLines("```", con)
}
