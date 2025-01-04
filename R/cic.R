
#' Pipe operate (credit magrittr package)
#' @importFrom magrittr `%>%`
`%>%` = magrittr::`%>%`

# to compute the counterfactual CDFs we need
# 1. the empirical CDF / inverse CDF for every group and year and error type
# 2. identifiers associated with each empirical CDF for
#   a. the outcome (treatment)
#   b. the year
#   c. the person
# 3. then we need the inversion logic for constructing the cf
# 4. finally we need a classifier that can map the cdfs to outcomes

#' Returns smallest x
min_or_zero = function(x) {
  if(is.null(x) | length(x) == 0) {
    0
  } else {
    min(x)
  }
}

max_or_zero = function(x) {
  if(is.null(x) | length(x) == 0) {
    0
  } else {
    max(x)
  }
}

#' Takes in three vectors and computes the CIC model for the
#' initial one
#' @param y_01 vector of outcomes for group 1 at time 0
#' @param y_00 vector of outcomes for group 0 at time 0
#' @param y_01 vector of outcomes for group 0 at time 1
#' @details assume that group 1 is the "treated" group.
#' We are interested in the distribution of outcomes
#' (F^N_11) that would have occurred at time 1 for the
#' treated group had the intervention "treatment" not occurred
cic = function(F_01, F_00, F_10) {
  cdf_11 = F_10(quantile(F_00, F_01(y_10)))
  k_cic = quantile(F_01(F_00(y_10)))
  list(
    ecdfs = list(
      F_01, F_00, F_10
    ),
    y_10 = y_10,
    cf_cdf = cdf_11,
    k_cic = k_cic
  )
}

#' takes in a list of data and computes empirical CDFs
#' @param l list of numeric vectors to pass to ecdf
#' @param ... other arguments to pass to ecdf
compute_cdfs = function(l, ...) {
  lapply(l, \(x) ecdf(x, ...))
}

#' Assert that N vectors are the same length
#' @param ... an arbitrary number of vectors
assert_same_length = function(...) {
  l = sapply(list(...), length)
  len = l[1]
  same_length = Reduce(`&`, l == len)
  stopifnot(same_length)
}

#' Identify counterfactuals from data
#' @param i a scalar / single identifier for which group to consider
#' @param t_i time for outcome y = 1
#' @param j index of identifiers for groups
#' @param t index of years
#' @param y outcomes for each group x year
#' @rdname id_cfs
#' @importFrom data.table
#' @export
id_cfs = function(i, t_i, j, t, y,
                  t_window = NA,
                  group = NA) {
  if(!is.na(t_window)) {
    window = which(t >= t_i & t >= t_i - t_window)
    j = j[window]
    t = t[window]
    y = y[window]
  }

  # for a group that is matched on observables
  if(!is.na(group)) {
    g_i = group[j == i & t == t_i]
    g_idx = which(group == g_i)
    j = j[g_idx]
    t = t[g_idx]
    y = y[g_idx]
  }

  pre_period_own = f_10(i, t_i, j, t, y)
  pre_period_other = f_00(i, t_i, j, t, y)
  same_period_other = f_01(i, t_i, j, t, y)
  same_period_own = f_11(i, t_i, j, t, y)

  all_cfs = rbind(pre_period_own,
                  pre_period_other,
                  same_period_other)
  all_cfs[, available_cf := .N,
          by = c("i", "t_i", "s", "j")]
  all_cfs = all_cfs[available_cf == 3]
  all_cfs[, available_cf := NULL]
  all_cfs[]
}

#' @rdname id_cfs
f_00 = function(i, t_i, j, t, y) {
  assert_same_length(j, t, y)
  idx = which(j != i & t < t_i & y == 0)
  data.table::data.table(i = i, t_i = t_i, j = j[idx],
                         s = t[idx], f = "00")
}

#' @rdname id_cfs
f_01 = function(i, t_i, j, t, y) {
  assert_same_length(j, t, y)
  idx = which(j != i & t == t_i & y == 1)
  data.table::data.table(i = i, t_i = t_i, j = j[idx],
                         s = t[idx], f = "01")
}

#' @rdname id_cfs
f_10 = function(i, t_i, j, t, y) {
  assert_same_length(j, t, y)
  idx = which(j == i & t < t_i & y == 0)
  data.table::data.table(i, t_i, j = j[idx],
                         s = t[idx], f = "10")
}

#' @rdname id_cfs
f_11 = function(i, t_i, j, t, y) {
  assert_same_length(j, t, y)
  idx = which(j == i & t == t_i & y == 1)
  data.table::data.table(i, t_i, j = j[idx], s = t[idx], f = "11")
}

#' Get all the counterfactuals for every person / time
#' @param i index of individuals
#' @param t index of time
#' @param y index of outcome
#' @param t_window time window to filter on
#' @param group group to match on (for matching on observables)
get_all_cfs = function(i, t, y, t_window = NA, group = NA) {
  ui = unique(i)
  cf_list = list()
  for(ii in 1:length(ui)) {
    t_ii_y1 = unique(t[i == ui[ii] & y == 1])
    if(length(t_ii_y1) == 0) {
      l = data.table::data.table()
    } else {
      l = list()
      for(tt in 1:length(t_ii_y1)) {
        l[[tt]] = id_cfs(ui[ii], t_ii_y1[tt],
                         i, t, y, t_window, group)
      }
    }

    names(l) = as.character(t_ii_y1)
    cf_list[[ii]] = l
  }
  names(cf_list) = ui
  cf_list
}
