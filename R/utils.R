## Read R Markdown file and split out sections into modules


is.even <- function(x) {
        !(x %% 2)
}

splitModules <- function(file) {
        txt <- readLines(file)
        isheader <- which(grepl("^##", txt, perl = TRUE))
        ischunk <- which(grepl("^```", txt, perl = TRUE))
        int <- findInterval(isheader, ischunk)
        headers <- isheader[is.even(int)]
        outdir <- "modules"
        outpath <- sub("^##\\s+", "", txt[headers], perl = TRUE)
        outpath <- gsub("\\s+", "_", outpath, perl = TRUE)
        folders <- file.path(outdir, outpath)
        dir.create(folders, recursive = TRUE)
}
