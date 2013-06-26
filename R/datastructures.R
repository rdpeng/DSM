################################################################################
## Data structures for environmental epi example (air pollution time series)
################################################################################

setOldClass("lm")
setOldClass("glm")

setClass("exposureModel",
         representation(model = "lm",
                        exposure = "character"))

setGeneric("getNames", function(x, ...) standardGeneric("getNames"))
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

setGeneric("modelType", function(object, ...) standardGeneric("modelType"))

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

setGeneric("describe", function(object, ...)
           standardGeneric("describe"))

setMethod("describe", "exposureModel", function(object, ...) {
        nms <- getNames(object)
        sprintf("a %s model where %s was regressed on exposure %s along with %s",
                modelType(object@model), sQuote(nms$response),
                sQuote(nms$exposure), char2comma(nms$other))
})

setGeneric("modelEffect", function(object, ...) standardGeneric("modelEffect"))

setGeneric("modelCoef", function(object, ...) standardGeneric("modelCoef"))

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


setClass("airpollutionTS",
         representation(timevar = "character",
                        "exposureModel"))
newAPTS <- function(model, exposure, timevar = NULL) {
        expm <- new("exposureModel", model = model, exposure = exposure)
        if(is.null(timevar)) {
                nms <- getNames(expm)
                i <- grep(c("date|time"), tolower(nms$other), perl = TRUE)
                if(length(i) > 0) {
                        timevar <- nms$other[i]
                        message(sprintf("using '%s' for 'timevar'", timevar))
                }
                else
                        stop("need to specify 'timevar'")
        }
        new("airpollutionTS", timevar = timevar, expm)
}

setGeneric("sensitivityAnalysis", function(object, ...)
           standardGeneric("sensitivityAnalysis"))

setClass("adjustTimeDF",
         representation(models = "list",
                        dfseq = "integer",
                        timeVec = "character",
                        orig = "airpollutionTS"))

adjustTimeDF <- function(object, dfseq, timevar = "date", smoothType = "ns") {
        dfseq <- as.integer(dfseq)
        nms <- getNames(object)
        pos <- grep(timevar, nms$other)
        origTimeVar <- nms$other[pos]
        timeVec <- paste("ns(", timevar, ", ", dfseq, ")", sep = "")
        results <- vector("list", length = length(dfseq))

        for(i in seq(along = dfseq)) {
                newFormula <- reformulate(c(nms$exposure, nms$other[-pos],
                                            timeVec[i]),
                                          response = nms$response)
                out <- update(object@model, formula = newFormula)
                expm <- new("exposureModel", model = out,
                            exposure = object@exposure)
                results[[i]] <- new("airpollutionTS", timevar = timevar, expm)
        }
        new("adjustTimeDF", models = results, dfseq = dfseq,
            timeVec = timeVec, orig = object)
}

setMethod("sensitivityAnalysis", "airpollutionTS", adjustTimeDF)


setMethod("describe", "adjustTimeDF", function(object, ...) {
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

setMethod("show", "adjustTimeDF",
          function(object) {
                  out <- lapply(describe(object), function(x) {
                          strwrap(paste0("* ", x), exdent = 2)
                  })
                  with(out, cat(model.orig, coef.orig, sens.range, sep = "\n"))
                  invisible(object)
          })

setGeneric("plot")
setGeneric("summary")

setMethod("plot", "adjustTimeDF", function(x, y, xlab = "Degrees of freedom for time", ylab = "Coefficient", ...) {
        nms <- getNames(x@orig)
        b <- sapply(x@models, function(x) coef(x@model)[nms$exposure])
        ci <- sapply(x@models, function(x) {
                confint.default(x@model, nms$exposure)
        })
        xpts <- x@dfseq
        ypts <- b
        plot(xpts, ypts, ylim = range(ci), xlab = xlab, ylab = ylab, ...)
        segments(xpts, ci[1, ], xpts, ci[2, ])
        abline(h = coef(x@orig@model)[nms$exposure], lty = 2)
})

################################################################################
## EDA

