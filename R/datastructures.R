################################################################################
## Data structures for environmental epi example (air pollution time series)
################################################################################

################################################################################
## Random utilities

library(xlsx)

## Read an Excel file with some help

read.excel <- function(file, ...) {
        d0 <- read.xlsx(file, ...)
        for(i in seq_len(ncol(d0))) {
                cl <- class(d0[[i]])
                if(cl == "factor" || cl == "character") {
                        v <- as.character(d0[[i]])
                        d0[[i]] <- type.convert(v)
                }
        }
        d0
}

################################################################################
## Data Classes

setOldClass("lm")
setOldClass("glm")
setOldClass("data.frame")
setClass("exposureModel",
         representation(model = "lm",
                        exposure = "character"))
setClass("sensitivityResults",
         representation(models = "list",
                        orig = "exposureModel"))
setClass("adjustTimeDF",
         representation(dfseq = "integer",
                        timeVec = "character",
                        "sensitivityResults"))
setClass("APTSData",
         representation(response = "character",
                        exposure = "character",
                        timevar = "character",
                        "data.frame"))
setClass("APTSSummary",
         representation(miss = "numeric", skew = "list", disp = "numeric",
                        highlev = "character", outlier = "character",
                        data = "APTSData"))
setClass("APTSModel",
         representation(timevar = "character",
                        summary = "APTSSummary",
                        "exposureModel"))

setClass("optimalTimeDF",
         representation(dfmin = "numeric", aic = "numeric",
                        dfseq = "numeric", method = "character"))

################################################################################
## Generic Functions

setGeneric("modelEffect", function(object, ...) standardGeneric("modelEffect"))
setGeneric("modelCoef", function(object, ...) standardGeneric("modelCoef"))
setGeneric("getNames", function(x, ...) standardGeneric("getNames"))
setGeneric("modelType", function(object, ...) standardGeneric("modelType"))
setGeneric("describe", function(object, ...) standardGeneric("describe"))
setGeneric("sensitivityAnalysis", function(object, ...)
           standardGeneric("sensitivityAnalysis"))
setGeneric("plot")
setGeneric("summary")
setGeneric("checkData", function(object, ...) standardGeneric("checkData"))
setGeneric("influential", function(object, ...) standardGeneric("influential"))
setGeneric("selectTimeDF", function(object, ...) standardGeneric("selectTimeDF"))
setGeneric("influenceAnalysis", function(object, ...) standardGeneric("influenceAnalysis"))

setGeneric("nyears", function(object, ...) standardGeneric("nyears"))

################################################################################
## Methods
################################################################################

################################################################################
## Utilities

setAs("optimalTimeDF", "numeric", function(from) {
        from@dfmin
})

setMethod("plot", "optimalTimeDF", function(x, y, pch = 19, ...) {
        plot(x@dfseq, x@aic, xlab = "Degrees of freedom per year",
             ylab = "AIC", pch = pch, ...)
})

setMethod("show", "optimalTimeDF", function(object) {
        show(object@dfmin)
})

setMethod("getNames", "exposureModel",
          function(x, ...) {
                  frm <- formula(x@model)
                  others <- attr(terms(x@model), "term.labels")
                  i <- grep(x@exposure, others)
                  list(response = as.character(frm)[2],
                       exposure = others[i],
                       other = others[-i])
          })

char2comma <- function(x) {
        n <- length(x)
        if(n <= 1L)
                x
        else if(n == 2L)
                paste(sQuote(x), collapse = " and ")
        else
                paste0(paste(sQuote(x[-n]), collapse = ", "), ", and ",
                       sQuote(x[n]))
}

## For now read an RDS file
readAPTSData <- function(file, ..., response = NULL, exposure = NULL,
                         timevar = NULL) {
        ## d0 <- read.csv(file, ...)
        d0 <- readRDS(file, ...)
        nms <- names(d0)
        if(is.null(response))
                response <- nms[1]
        if(is.null(exposure))
                exposure <- nms[2]
        if(is.null(timevar)) {
                i <- grep(c("date|time"), tolower(nms), perl = TRUE)
                if(length(i) > 0) {
                        timevar <- nms[i]
                        message(sprintf("using '%s' for 'timevar'", timevar))
                }
                else
                        stop("need to specify 'timevar'")
        }
        stopifnot(response %in% nms)
        stopifnot(exposure %in% nms)
        stopifnot(timevar %in% nms)
        new("APTSData", d0, response = response, exposure = exposure,
            timevar = timevar, d0)
}

setMethod("nyears", "Date", function(object, ...) {
        diff(range(as.POSIXlt(object)$year)) + 1
})

setMethod("nyears", "APTSModel", function(object, ...) {
        timevar <- object@timevar
        dat <- as(object@summary@data, "data.frame")
        nyears(dat[, timevar])
})

################################################################################
## EDA

library(moments)

## Check for missing data in exposure, missing overall, skewness in
## predictors, overdispersion in outcome, outliers; Generate a
## report/summary
setMethod("checkData", "APTSData", function(object, ...) {
        nms <- names(object)
        ## Missing data
        miss.exp <- mean(is.na(object[, object@exposure]))
        miss.response <- mean(is.na(object[, object@response]))
        comp <- mean(complete.cases(object))
        miss <- c(response = miss.response, exposure = miss.exp,
                  complete = comp)
        ## Right skewness in predictors
        use <- setdiff(nms, object@response)
        skewvar <- lapply(use, function(vname) {
                try(suppressWarnings(agostino.test(object[, vname],
                                                   alt = "less")),
                    silent = TRUE)
        })
        u <- !sapply(skewvar, inherits, what = "try-error")
        skewvar <- skewvar[u]
        names(skewvar) <- use[u]
        ## Overdispersion
        disp <- var(object[, object@response]) / mean(object[, object@response])
        ## High leverage predictors (continuous)
        use <- setdiff(nms, object@response)
        tempform <- reformulate(use)
        mm <- model.matrix(tempform, as(object, "data.frame"))
        lev <- diag(mm %*% tcrossprod(solve(crossprod(mm)), mm))
        highlev <- names(which(lev >= (mean(lev) * 4)))
        ## Response outliers
        dummyform <- reformulate(sprintf("ns(%s, %d)", object@timevar,
                                         nyears(object[, object@timevar]) * 4L),
                                 object@response)
        dummy <- glm(dummyform, data = as(object, "data.frame"),
                     family = poisson, na.action = na.exclude)
        res <- rstandard(dummy)
        out <- names(which(res > 6))

        new("APTSSummary", miss = miss, skew = skewvar, disp = disp,
            highlev = highlev, outlier = out, data = object)
})

## Return obs. IDs (row #s) of influential observations
setMethod("influential", "APTSSummary", function(object, ...) {
        unique(c(object@highlev, object@outlier))
})


setMethod("describe", "APTSSummary", function(object, ...) {
        miss.msg <- sprintf("missing: %.2f%% of response, %.2f%% of exposure", 100*object@miss["response"], 100*object@miss["exposure"])
        cc.msg <- sprintf("complete: %.2f%% of observations were complete cases", 100*object@miss["complete"])
        p <- p.adjust(sapply(object@skew, "[[", "p.value"), "BH")
        p <- p[!is.na(p)]
        sk <- p < 0.05
        skew.msg <- sprintf("right-skew: %s", ifelse(all(!sk), "NA", paste(names(which(sk)), collapse = ", ")))
        lev.msg <- sprintf("high leverage (row #): %s", paste(object@highlev, collapse = ", "))
        out.msg <- sprintf("outlier (row #): %s", paste(object@outlier, collapse = ", "))
        c(miss = miss.msg, complete = cc.msg, skew = skew.msg,
          highlev = lev.msg, outlier = out.msg)
})

setMethod("show", "APTSSummary", function(object) {
        d <- describe(object)
        x <- strwrap(paste0("* ", d), exdent = 2)
        writeLines(x)
})

################################################################################
## Model Building

setMethod("selectTimeDF", "APTSModel", selectTimeDF)

selectTimeDF <- function(object, dfseq = 2:16, ...) {
        nms <- getNames(object)
        mod <- object@model
        pos <- grep(object@timevar, nms$other)
        dat <- as(object@summary@data, "data.frame")
        nyr <- nyears(dat[, object@timevar])
        dfseq.total <- dfseq * nyr
        timeVec <- paste0("ns(", object@timevar, ", ", dfseq.total, ")")
        results <- lapply(timeVec, function(tv) {
                newFormula <- reformulate(c(tv, nms$other[-pos]),
                                          response = object@exposure)
                try(update(mod, formula = newFormula, family = gaussian),
                    silent = TRUE)
        })
        aic <- sapply(results, function(x) try(AIC(x), silent = TRUE))
        names(aic) <- dfseq
        dfmin <- dfseq[which.min(aic)]
        new("optimalTimeDF", dfmin = dfmin, aic = aic, dfseq = dfseq)
}

################################################################################
## Model Reporting

setMethod("modelType", "lm",
          function(object, ...) {
                  "linear regression"
          })

setMethod("modelType", "glm",
          function(object, ...) {
                  fam <- family(object)
                  if(fam$family == "gaussian")
                          "linear regression"
                  else if(fam$family == "binomial") 
                          paste(fam$link, "regression")
                  else if(fam$family == "poisson")
                          "Poisson regression"
                  else 
                          paste(fam$family, "regression")                  
          })


setMethod("describe", "exposureModel", function(object, ...) {
        nms <- getNames(object)
        sprintf("a %s model where %s was regressed on exposure %s along with %s",
                modelType(object@model), sQuote(nms$response),
                sQuote(nms$exposure), char2comma(nms$other))
})

setMethod("modelCoef", "exposureModel", function(object, incr = 1, ...) {
        fam <- family(object@model)$family
        nms <- getNames(object)
        ci <- suppressMessages(confint.default(object@model, nms$exposure))
        ci <- ci * incr
        cf <- coef(object@model)
        i <- grep(object@exposure, names(cf))
        b <- cf[i]
        bb <- b * incr

        if(fam == "binomial" || fam == "poisson") {
                bb <- exp(bb)
                ci <- exp(ci)
        }
        ci.txt <- prettyNum(ci)
        sprintf("%s (95%% CI: %s, %s)", prettyNum(abs(bb)),
                ci.txt[1], ci.txt[2])
})

modelEffect.method <- function(object, incr = 1, pvalue = FALSE, ...) {
        fam <- family(object@model)$family
        nms <- getNames(object)
        cf <- coef(object@model)
        i <- grep(object@exposure, names(cf))
        si <- sign(cf[i])

        if(fam == "binomial") {
                msg <- sprintf("a %s unit %s in %s was associated with a %s %s odds of %s",
                               abs(incr), ifelse(incr > 0, "increase", "decrease"),
                               sQuote(nms$exposure), modelCoef(object, incr),
                               ifelse(si > 0, "increased", "decreased"),
                               sQuote(nms$response))
        }
        else if(fam == "poisson") {
                msg <- sprintf("a %s unit %s in %s was associated with a %s %s rate of %s",
                               abs(incr), ifelse(incr > 0, "increase", "decrease"),
                               sQuote(nms$exposure), modelCoef(object, incr),
                               ifelse(si > 0, "increased", "decreased"),
                               sQuote(nms$response))
        }
        else {
                msg <- sprintf("a %s unit %s in %s was associated with a %s unit %s in %s",
                               abs(incr), ifelse(incr > 0, "increase", "decrease"),
                               sQuote(nms$exposure), modelCoef(object, incr),
                               ifelse(si > 0, "increase", "decrease"),
                               sQuote(nms$response))
        }
        if(pvalue) {
                p <- summary(object@model)$coefficients[nms$exposure, 4]
                msg <- sprintf("%s (p = %s)", msg, prettyNum(p))
        }
        msg
}


setMethod("modelEffect", "exposureModel", modelEffect.method)

setMethod("show", "exposureModel",
          function(object) {
                  out <- paste0("* ", describe(object))
                  writeLines(strwrap(out, exdent = 2))
                  invisible(object)
          })


################################################################################
## Sensitivity Analysis

setMethod("influenceAnalysis", "APTSModel", function(object, ...) {
        ## Influential observations
        fit <- object@model
        idx <- influential(object@summary)
        dat <- as(object@summary@data, "data.frame")
        expo <- getNames(object)$exposure
        modellist <- lapply(idx, function(id) {
                i <- match(id, row.names(dat))
                newdat <- dat[-i, ]
                m <- try(update(fit, data = newdat), silent = TRUE)
                new("exposureModel", model = m, exposure = object@exposure)
        })
        rm.all <- try({
                i <- match(idx, row.names(dat))
                new("exposureModel", model = update(fit, data = dat[-i, ]),
                    exposure = object@exposure)
        }, silent = TRUE)
        results <- append(modellist, list(rm.all))
        names(results) <- c(idx, paste(idx, collapse = "."))
        new("sensitivityResults", models = results, orig = object)
})

setMethod("sensitivityAnalysis", "APTSModel", function(object, dfseq, ...) {
        dfseq <- as.integer(dfseq)
        dfseq.total <- dfseq * nyears(object)
        nms <- getNames(object)
        timevar <- object@timevar
        pos <- grep(timevar, nms$other)
        origTimeVar <- nms$other[pos]
        timeVec <- paste("ns(", timevar, ", ", dfseq.total, ")", sep = "")
        results <- vector("list", length = length(dfseq.total))

        for(i in seq(along = dfseq.total)) {
                newFormula <- reformulate(c(nms$exposure, nms$other[-pos],
                                            timeVec[i]),
                                          response = nms$response)
                out <- update(object@model, formula = newFormula)
                expm <- new("exposureModel", model = out,
                            exposure = object@exposure)
                results[[i]] <- expm
        }
        names(results) <- as.character(dfseq)
        new("adjustTimeDF", models = results, dfseq = dfseq,
            timeVec = timeVec, orig = object)
})

setMethod("describe", "sensitivityResults", function(object, ...) {
        nms <- getNames(object@orig)
        b <- sapply(object@models, function(x) {
                coef(x@model)[nms$exposure]
        })
        imin <- which.min(b)
        imax <- which.max(b)
        msg1 <- sprintf("a sensitivity analysis for %s",
                        describe(object@orig))
        msg2 <- modelEffect(object@orig)
        msg3 <- sprintf("the coefficient for %s varied from %s to %s",
                        sQuote(nms$exposure),
                        modelCoef(object@models[[imin]]),
                        modelCoef(object@models[[imax]]))
        list(model.orig = msg1, coef.orig = msg2, sens.range = msg3)
})

setMethod("show", "sensitivityResults",
          function(object) {
                  out <- lapply(describe(object), function(x) {
                          strwrap(paste0("* ", x), exdent = 2)
                  })
                  with(out, cat(model.orig, coef.orig, sens.range, sep = "\n"))
                  invisible(object)
          })

setMethod("plot", "sensitivityResults", function(x, y, xlab = "Models", ylab = "Coefficient", pch = 19, legend = TRUE, ...) {
        nms <- getNames(x@orig)
        b <- sapply(x@models, function(x) coef(x@model)[nms$exposure])
        ci <- sapply(x@models, function(x) {
                confint.default(x@model, nms$exposure)
        })
        xpts <- seq(along = x@models)
        ypts <- b
        plot(xpts, ypts, ylim = range(ci), xlab = xlab, ylab = ylab,
             pch = pch, xaxt = "n", ...)
        segments(xpts, ci[1, ], xpts, ci[2, ])
        abline(h = coef(x@orig@model)[nms$exposure], lty = 2)
        if(legend)
                legend("topright", legend = "Original Estimate", lty = 2,
                       bty = "n")
        axis(1, xpts, names(x@models))
        invisible()
})

setMethod("plot", "adjustTimeDF", function(x, y, xlab = "Degrees of freedom for time", ylab = "Coefficient", pch = 19, legend = TRUE, ...) {
        nms <- getNames(x@orig)
        b <- sapply(x@models, function(x) coef(x@model)[nms$exposure])
        ci <- sapply(x@models, function(x) {
                confint.default(x@model, nms$exposure)
        })
        xpts <- x@dfseq
        ypts <- b
        plot(xpts, ypts, ylim = range(ci), xlab = xlab, ylab = ylab,
             pch = pch, ...)
        segments(xpts, ci[1, ], xpts, ci[2, ])
        abline(h = coef(x@orig@model)[nms$exposure], lty = 2)
        if(legend)
                legend("topright", legend = "Original Estimate", lty = 2,
                       bty = "n")
        invisible()
})

