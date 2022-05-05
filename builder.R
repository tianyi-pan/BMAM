
.namesr <- function(id, usecoef, usegroup) {
  n <- expand.grid(ID = id,
                   Coef = usecoef,
                   Group = usegroup)
  n <- as.data.table(n)
  n[, sprintf("r_%s[%d,%s]",
              Group, ID, Coef)]
}



.namesz <- function(id, usecoef, usegroup) {
  n <- expand.grid(ID = id,
                   Coef = seq_along(usecoef),
                   Group = seq_along(usegroup))
  n <- as.data.table(n)
  n[, sprintf("z_%d[%d,%d]",
              Group, Coef, ID)]
}


.buildr <- function(data, id, usecoef,usegroup) {
  stopifnot(is.data.table(data))
  id <- unique(id)
  n <- .namesr(id, usecoef, usegroup)
  out <- as.data.frame(data[, ..n])
  as.list(out)
}

.buildz <- function(data, id, usecoef,usegroup) {
  stopifnot(is.data.table(data))
  id <- unique(id)
  n <- .namesz(id, usecoef, usegroup)
  out <- as.data.frame(data[, ..n])
  as.list(out)
}


.buildJ <- function(data, block) {
  n <- sprintf("J_%d", block)
  data[n]
}
