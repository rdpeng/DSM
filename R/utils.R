## Read R Markdown file and split out sections into modules


is.even <- function(x) {
        !(x %% 2)
}

splitModules <- function(file) {
        if(!require(knitr))
                stop("need 'knitr' package")
        txt <- readLines(file)
        isheader <- which(grepl("^##", txt, perl = TRUE))
        ischunk <- which(grepl("^```", txt, perl = TRUE))
        int <- findInterval(isheader, ischunk)
        headers <- isheader[is.even(int)]
        outdir <- "modules"
        outpath <- sub("^##\\s+", "", txt[headers], perl = TRUE)
        outpath <- gsub("\\s+", "_", outpath, perl = TRUE)
        ord <- seq_along(outpath)
        prefix <- formatC(ord, flag = "0", width = max(nchar(ord)))
        outpath <- paste(prefix, outpath, sep = "_")
        folders <- file.path(outdir, outpath)
        r <- sapply(folders, dir.create, recursive = TRUE, showWarnings = FALSE)
        if(!all(r))
                warning("not all module directories could be created")
        preamble <- txt[seq_len(min(headers) - 1)]
        preambledir <- file.path(outdir, "00_Preamble")
        dir.create(preambledir, recursive = TRUE)
        writeLines(preamble, file.path(preambledir, "code.Rmd"))
        
        ## Split file into separate directories
        len <- c(diff(headers), length(txt) - headers[length(headers)] + 1)
        for(i in seq_along(headers)) {
                output <- txt[seq(headers[i], headers[i] + len[i] - 1)]
                writeLines(output, file.path(folders[i], "code.Rmd"))
        }

        ## Run modules and gather input/output data
        env <- new.env()
        input <- ls.str(env)
        knit(text = preamble, envir = env, quiet = TRUE)
        output <- ls.str(env)
        saveRDS(input, file.path(preambledir, "input.rds"))
        saveRDS(output, file.path(preambledir, "output.rds"))
        input <- output
        for(i in seq_along(folders)) {
                message(folders[i])
                saveRDS(input, file.path(folders[i], "input.rds"))
                inputtxt <- txt[seq(headers[i], headers[i] + len[i] - 1)]
                knit(text = inputtxt, envir = env, quiet = TRUE)
                output <- ls.str(env)
                saveRDS(output, file.path(folders[i], "output.rds"))
                input <- output
        }
}

assembleModules <- function(moduledir, file = NULL) {
        folders <- dir(moduledir, recursive = FALSE, full.names = TRUE)
        txtlist <- lapply(folders, function(f) {
                readLines(file.path(f, "code.Rmd"))
        })
        txt <- unlist(txtlist)
        if(!is.null(file))
                writeLines(txt, file)
        invisible(txt)
}
