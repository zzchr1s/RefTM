##*******************************************************
## Classes TMcontrol1, LDAcontrol1
##
## control parameters for the lda functions
## + superclass (TMcontol1)


##**********************************************************
## coercion
setAs("NULL", "LDA_VEMcontrol1", function(from, to) new(to))
setAs("NULL", "LDA_Gibbscontrol1", function(from, to) new(to))
setAs("NULL", "OPTcontrol1", function(from, to) new(to))

setAs("list", "LDA_VEMcontrol1", function(from, to) .list2control(from, to))
setAs("list", "LDA_Gibbscontrol1", function(from, to) .list2control(from, to))

setAs("list", "OPTcontrol1", function(from, to) .list2control(from, to))

.list2control <- function(from, to) {
  n = names(from)
  s = slotNames(to)
  p = pmatch(n, s)
  if(any(is.na(p)))
    stop(paste("Invalid slot name(s) for class",
               to, ":", paste(n[is.na(p)], collapse=" ")))

  prototype <- new(to)
  slotTypes <- getClass(to)@slots[p]
  from <- lapply(seq_along(from), function(i) {
    if (slotTypes[[i]] == "OPTcontrol1") {
      defaults <- slot(prototype, s[p][i])
      miss <- which(!slotNames(defaults) %in% names(from[[i]]))
      if (length(miss) > 0) {
        from[[i]] <- c(from[[i]],
                       sapply(slotNames(defaults)[miss],
                              function(s) slot(defaults, s), simplify = FALSE))
      }
    }
    as(from[[i]], slotTypes[[i]])
  })

  names(from) = s[p]
  do.call("new", c(from, Class=to))
}


