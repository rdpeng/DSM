################################################################################
## Data structures for environmental epi example
################################################################################

setGeneric("units")
setGeneric("units<-")

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
                  else if(fam$family == "binomial") {
                          paste(fam$link, "regression")
                  }
                  else 
                          paste(fam$family, "regression")                  
          })

setGeneric("describeModel", function(object, ...)
           standardGeneric("describeModel"))

setMethod("describeModel", "modelOutput",
          function(object, ...) {
                  nms <- getNames(object)
                  sprintf("in this %s model, %s was regressed on %s",
                          modelType(object@model), sQuote(nms$response),
                          char2comma(nms$other))
          })

setGeneric("modelEffect", function(object, ...) standardGeneric("modelEffect"))

setMethod("modelEffect", "modelOutput",
          function(object, incr = 1, ...) {
                  b <- coef(object@model)[object@exposure]
                  bb <- b * incr
                  sprintf("a %s unit %s in %s was associated with a %s unit %s in %s",
                          abs(incr), ifelse(incr > 0, "increase", "decrease"),
                          sQuote(object@exposure), prettyNum(abs(bb)),
                          ifelse(bb > 0, "increase", "decrease"),
                          sQuote(getNames(object)$response))
          })

startcaps <- function(x) {

}

setMethod("show", "modelOutput",
          function(object) {
                  out <- describeModel(object)
                  cat(out, ".\n", sep = "")
                  invisible(object)
          })
