##**********************************************************
## control parameters

setClass("OPTcontrol1",
         representation(iter.max = "integer",
                        tol      = "numeric"),
         prototype(iter.max = -1L,
                   tol      = sqrt(.Machine$double.eps)))

setClass("TopicModelcontrol1",
         representation(seed          = "integer",
                        verbose       = "integer",
                        prefix        = "character",
                        save          = "integer",
                        nstart        = "integer",
                        best          = "logical",
                        keep          = "integer",
                        estimate.beta = "logical",
                        "VIRTUAL"),
         prototype(verbose = 0L,
                   save = 0L,
                   best = TRUE,
                   keep = 0L,
                   estimate.beta = TRUE))

init_TopicModelcontrol <- function(.Object, prefix, seed, nstart, ...) {
  if (missing(prefix)) prefix <- tempfile()
  if (nchar(prefix) > 200) {
      stop("prefix is too long, please use a shorter prefix")
  }
  if (missing(seed)) {
    if (missing(nstart)) {
      nstart <- 1L
      seed <- as.integer(Sys.time())
    } else {
      seed <- sample(seq_len(10^6), nstart)
    }
  } else if (missing(nstart)) nstart <- length(seed)

  list(.Object = .Object, prefix = prefix, seed = seed, nstart = nstart, ... = ...)
}

setMethod("initialize", "TopicModelcontrol1", function(.Object, prefix, seed, nstart, ...) {
  args <- init_TopicModelcontrol(.Object, prefix, seed, nstart, ...)
  .Object <- do.call("callNextMethod", args)
  invisible(.Object)
})

setClass("VEMcontrol1",
         representation(var    = "OPTcontrol1",
                        em     = "OPTcontrol1",
                        initialize = "character",
                        "VIRTUAL"),
         prototype(var        = new("OPTcontrol1", iter.max = 500L, tol = 10^-6),
                   em         = new("OPTcontrol1", iter.max = 1000L, tol = 10^-4),
                   initialize = "random"))

setMethod("initialize", "VEMcontrol1", function(.Object, initialize = "random", ...) {
  initialize <- match.arg(initialize, c("random", "seeded", "model"))
  args <- init_TopicModelcontrol(.Object, ...)
  .Object <- do.call("callNextMethod",
                     c(args, initialize = initialize))
  invisible(.Object)
})

setClass("LDAcontrol1",
         representation(alpha = "numeric",
                        "VIRTUAL"),
         contains = "TopicModelcontrol1")

setClass("LDA_VEMcontrol1",
         representation(estimate.alpha   = "logical"),
         contains = c("LDAcontrol1", "VEMcontrol1"),
         prototype(estimate.alpha = TRUE))

setMethod("initialize", "LDA_VEMcontrol1", function(.Object, prefix, initialize = "random", ...) {
  if (missing(prefix)) prefix <- tempfile()
  .Object <- callNextMethod(.Object = .Object, prefix = prefix, initialize = initialize, ...)
  invisible(.Object)
})

setClass("LDA_Gibbscontrol1",
    representation(delta  = "numeric",
                   iter   = "integer",
                   thin   = "integer",
                   burnin = "integer",
                   initialize = "character"),
         contains = "LDAcontrol1",
         prototype(delta   = 0.1,
                   verbose = 0L,
                   iter    = 2000L,
                   burnin  = 0L,
                   nstart  = 1L,
                   best    = TRUE,
                   initialize = "random"))

setMethod("initialize", "LDA_Gibbscontrol1", function(.Object, initialize = "random", seed = as.integer(NA), ...) {
  initialize <- match.arg(initialize, c("random", "beta", "z"))
  .Object <- callNextMethod(.Object = .Object, initialize = initialize, seed = seed, ...)
  if (length(.Object@thin) == 0) .Object@thin <- .Object@iter
  invisible(.Object)
})

##**********************************************************
## Topic Models Objects

setClass("TopicModel1",
   representation(call            = "call",
                  Dim             = "integer",
                  control         = "TopicModelcontrol1",
                  k               = "integer",
                  terms           = "ANY",
                  documents       = "ANY",
                  beta            = "matrix",
                  gamma           = "matrix",
                  wordassignments = "ANY",
                  loglikelihood   = "numeric",
                  iter            = "integer",
                  logLiks         = "numeric",
                  n               = "integer",
                  "VIRTUAL"))

setClass("VEM1",
         contains = "TopicModel1",
         representation("VIRTUAL"))

setClass("RefTM.LDA",
         representation(alpha = "numeric",
                        "VIRTUAL"),
         contains = "TopicModel1")

setClass("LDA_VEM1",
         representation(),
         contains = c("RefTM.LDA", "VEM1"),
         prototype(control = new("LDA_VEMcontrol1")))

setClass("Gibbs1",
         contains = "TopicModel1",
         representation("VIRTUAL"))

setClass("LDA_Gibbs1",
         representation(seedwords = "ANY",
                        z = "integer"),
         contains = c("RefTM.LDA", "Gibbs1"),
         prototype(control = new("LDA_Gibbscontrol1")))

setClass("Gibbs_list1",
         representation(fitted = "list"))


setMethod("show", signature(object = "TopicModel1"), function(object) {
  cat("A", class(object), "topic model with", object@k, "topics.\n")
})

