
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

#' Takes in three vectors and computes the CIC model for the
#' initial one
#' @param y_01 vector of outcomes for group 1 at time 0
#' @param y_00 vector of outcomes for group 0 at time 0
#' @param y_01 vector of outcomes for group 0 at time 1
#' @details assume that group 1 is the "treated" group.
#' We are interested in the distribution of outcomes
#' (F^N_11) that would have occurred at time 1 for the
#' treated group had the intervention "treatment" not occurred
cic_from_wecdfs = function(F_01, F_00, F_10) {
  # restrict support to y_01 per corr 3.1
  y_01 = get_x(F_01)
  y_10 = get_x(F_10)
  y_10 = y_10[y_10 < max(y_01) & y_10 > min(y_01)]
  subset_ids = which(y_10 < max(y_01) & y_10 > min(y_01))
  cdf_cf_11 = F_10(quantile(F_00, F_01(y_10)))
  q = seq(0, 1, by = 0.01)
  if(length(cdf_cf_11) > 1) {
    y_11_N = stepfun(c(cdf_cf_11), c(y_10, max(y_10)), right = TRUE)(q)
  } else {
    y_11_N = NULL
  }

  k_cic = quantile(F_01, F_00(y_10))
  structure(
    list(
      wecdfs = list(
        F_01 = F_01, F_00 = F_00, F_10 = F_10
      ),
      y_10 = y_10,
      y_11_N = y_11_N,
      cf_cdf = cdf_cf_11,
      k_cic = k_cic,
      subset_ids = subset_ids
    ),
    class = "cic"
  )
}

#' Set parallel settings
#' @param n_workers number of workers
#' @param ... other arguments passed to future::plan
#' @importFrom future plan
#' @importFrom future sequential
#' @importFrom future multisession
set_parallel = function(n_workers, ...) {
  if(n_workers <= 1) {
    future::plan(future::sequential(...))
  } else {
    future::plan(future::multisession(workers = n_workers, ...))
  }
}

#' takes in a list of data and computes empirical CDFs
#' @param l list of numeric vectors to pass to wecdf
#' @param ... other arguments to pass to wecdf
#' @importFrom future.apply future_lapply
compute_cdfs = function(l, ...) {
  future.apply::future_lapply(l, \(x) wecdf(x, ...))
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
#' @param g_i group of the individual i for matching on observables
#' @param j index of identifiers for groups
#' @param t index of years
#' @param y outcomes for each group x year
#' @rdname id_cfs
#' @importFrom data.table `:=`
#' @details Note that y_i = 1 by convention, this function assumes
#' that indexing.
id_cfs = function(i, t_i, g_i,
                  j, t, g, y,
                  t_window = NA) {
  if(!is.na(t_window)) {
    window = which(t <= t_i & t >= (t_i - t_window))
    j = j[window]
    t = t[window]
    y = y[window]
    g = g[window]
  }

  pre_period_own = js_10(i, t_i, g_i, j, t, y, g)
  pre_period_other = js_00(i, t_i, g_i, j, t, y, g)
  same_period_other = js_01(i, t_i, g_i, j, t, y, g)
  same_period_own = js_11(i, t_i, g_i, j, t, y, g)

  # only consider comparisons where s is available
  # for both the donor person and the target person
  own_s = pre_period_own$s
  if(nrow(pre_period_other) > 0) {
    pre_period_other = pre_period_other[s %in% own_s]
  }

  # only consider the contemp period for individuals
  # who are also available in the pre-period, s
  pre_period_j = pre_period_other$j

  if(nrow(same_period_other) > 0) {
    same_period_other = same_period_other[j %in% pre_period_j]
  }

  all_comparison_ids = rbind(pre_period_own,
                             pre_period_other,
                             same_period_other,
                             same_period_own)
  all_comparison_ids[]
}

#' @rdname id_cfs
js_00 = function(i, t_i, g_i, j, t, y, g) {
  assert_same_length(j, t, y)
  idx = which(j != i & t < t_i & y == 0 & g == g_i)
  if(length(idx) > 0) {
    data.table::data.table(i = i, t_i = t_i, j = j[idx],
                           s = t[idx], f = "00")
  } else {
    data.table::data.table()
  }

}

#' @rdname id_cfs
js_01 = function(i, t_i, g_i, j, t, y, g) {
  assert_same_length(j, t, y)
  idx = which(j != i & t == t_i & y == 0 & g == g_i)
  if(length(idx) > 0) {
    data.table::data.table(i = i, t_i = t_i, j = j[idx],
                           s = t[idx], f = "01")
  } else {
    data.table::data.table()
  }
}

#' @rdname id_cfs
js_10 = function(i, t_i, g_i, j, t, y, g) {
  assert_same_length(j, t, y)
  idx = which(j == i & t < t_i & y == 0 & g == g_i)
  if(length(idx) > 0) {
    data.table::data.table(i, t_i, j = j[idx],
                           s = t[idx], f = "10")
  } else {
    data.table::data.table()
  }
}

#' @rdname id_cfs
js_11 = function(i, t_i, g_i, j, t, y, g) {
  assert_same_length(j, t, y)
  idx = which(j == i & t == t_i & y == 1 & g == g_i)
  if(length(idx) > 0) {
    data.table::data.table(i, t_i,
                           j = j[idx], s = t[idx],
                           f = "11")
  } else {
    data.table::data.table()
  }
}

#' Get the observations for i and t such that y = 1
#' @param i person index
#' @param t time index
get_itg_y1 = function(i, t, g, y) {
  data.table::data.table(
    i = i[y == 1],
    t = t[y == 1],
    g = g[y == 1]
  )
}

#' Get all the counterfactuals for every person / time
#' @param i index of individuals
#' @param t index of time
#' @param y index of outcome
#' @param t_window time window to filter on
#' @param group group to match on (for matching on observables)
#' @importFrom data.table rbindlist
get_all_cfs = function(i, t, g, y, t_window = NA) {

  # data.table of all person / time / group observations
  # where y = 1
  itg_y1 = get_itg_y1(i, t, g, y)
  # for each of these combinations we get counterfactuals
  cf_list = list()
  for(n in 1:nrow(itg_y1)) {
    i_n = itg_y1[n, i]
    t_n = itg_y1[n, t]
    g_n = itg_y1[n, g]
    cf_list[[n]] =
      id_cfs(i_n, t_n, g_n,
             i, t, g, y,
             t_window)
  }
  cf_list = data.table::rbindlist(cf_list)

  cf_list[]
}

#' Gets a list of CDFs for every group
#' @param i index of individual
#' @param t index for time
#' @param g index for group
#' @param x variable to construct the CDF of
#' @importFrom data.table data.table
build_cdfs = function(i, t, x, ...) {
  data.table::data.table(i, t, x) %>%
    split(by = c("i", "t")) %>%
    lapply(\(l) l$x) %>%
    compute_cdfs(...)
}

#' Build CIC Data
#' @param i index of individual
#' @param t index for time
#' @param g index for group
#' @param y index for binary outcome
#' @param x variable to construct the CDF of
#' @importFrom data.table merge.data.table
#' @importFrom data.table data.table
#' @importFrom data.table `:=`
build_cic_data = function(i, t, g, y, x, t_window = NA) {
  cdf_data = data.table::data.table(
    i = i, t = t, x = x
  )

  training_data = cdf_data
  training_data[, y := y]
  training_data = training_data[y == 1]

  cf_data = unique(
    data.table::data.table(
      i, t, g, y
    )
  )

  cfs = get_all_cfs(cf_data$i, cf_data$t, cf_data$g, cf_data$y, t_window)

  data_ids = unique(cfs[,.(i = j, t = s)])

  relevant_cdfs = merge(
    cdf_data,
    data_ids,
    by = c("i","t")
  )

  cdfs = build_cdfs(
    relevant_cdfs$i,
    relevant_cdfs$t,
    relevant_cdfs$x
  )
  list(
    counterfactuals = cfs,
    cdfs = cdfs,
    training_data = training_data
  )
}

#' Run CIC given the build_cic_data output
#' @param cic_data list output from cic_data
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor
#' @importFrom progressr with_progress
cic_from_cicdata = function(cic_data, n_boots) {
  cf_data = cic_data$counterfactuals

  # only consider cases where we have all four necessary groups
  # must be a faster way to do this lol
  cf_data[, cf_available := paste0(sort(unique(f)), collapse = ", "),
          by = c("i", "t_i", "j")]
  cf_data = cf_data[cf_available == "00, 01" | cf_available == "10, 11"]
  cf_data[, all_avail := paste0(sort(unique(f)), collapse = ", "),
          by = c("i", "t_i")]
  cf_data = cf_data[all_avail == "00, 01, 10, 11"]
  unique_i = unique(cf_data$i)
  # loop over people
  p = progressr::progressor(along = unique_i)

  l = future.apply::future_lapply(1:length(unique_i), \(ui) {
    ii = unique_i[ui]
    list_index = 1
    p()
    cf_data_i = cf_data[i == ii]

    # loop over their promotion dates
    unique_t_i = unique(cf_data_i$t_i)

    il = list()
    for(uti in 1:length(unique_t_i)) {
      tt = unique_t_i[uti]
      cf_data_it = cf_data[i == ii & t_i == tt]
      # there's only one of these for each i
      f11 = get_cdf(cic_data$cdfs,
                    cf_data_it[f=="11"]$j,
                    cf_data_it[f=="11"]$s)

      unique_s_10 = unique(cf_data_it[f == "10"]$s)
      unique_s_00 = unique(cf_data_it[f == "00"]$s)
      unique_s = intersect(unique_s_10, unique_s_00)
      # loop over counterfactual periods
      for(us in 1:length(unique_s)) {
        ss = unique_s[us]

        # only one of these for each i / s
        f10 = get_cdf(cic_data$cdfs,
                      cf_data_it[f=="10" & s == ss]$j,
                      cf_data_it[f=="10" & s == ss]$s)
        unique_j_01 = unique(cf_data_it[f == "01" & s == tt]$j)
        unique_j_00 = unique(cf_data_it[f == "00" & s == ss]$j)
        unique_j = intersect(unique_j_01, unique_j_00)

        for (uj in 1:length(unique_j)) {
          jj = unique_j[uj]
          # untreated group in the promotion year
          f01_idx = cf_data_it[s == t_i & j == jj & f == "01"]

          # untreated group before the promotion year
          f00_idx = cf_data_it[s == ss & j == jj & f == "00"]

          if(nrow(f00_idx) == 1 & nrow(f01_idx) == 1) {

            f01 = get_cdf(cic_data$cdfs,
                          f01_idx$j,
                          f01_idx$s)

            f00 = get_cdf(cic_data$cdfs,
                          f00_idx$j,
                          f00_idx$s)

            cl = cic_from_wecdfs(f01, f00, f10)

            y_11 = get_x(f11)

            q = seq(0, 1, by = 0.01)
            y_11_I = quantile(y_11, q)
            cl_dt = data.table(q = q,
                               x_11_I = y_11_I,
                               x_11_N = cl$y_11_N,
                               i = ii,
                               t = tt,
                               s = ss,
                               j = jj)
            il[[list_index]] = list(
              source_data = cl,
              table_summary = cl_dt
            )

            boot_dt = list()
            if(n_boots > 0) {
              for(b in 1:n_boots) {
                f00_b = randweight_cdf(f00)
                f01_b = randweight_cdf(f01)
                f10_b = randweight_cdf(f10)
                cl_b = cic_from_wecdfs(f01, f00, f10)
                boot_dt[[b]] =
                  data.table(q = q,
                             x_11_I = y_11_I,
                             x_11_N = cl$y_11_N,
                             i = ii,
                             t = tt,
                             s = ss,
                             j = jj,
                             b = b)
              }
              boot_dt = rbindlist(boot_dt)
              il[[list_index]]$bootstraps = boot_dt
            }
            list_index = list_index + 1

          }
        }
      }
    }
    il
  }, future.seed = T)

  l = unlist(l, recursive = F)

  tbls = lapply(l, \(x) x$table_summary) %>%
    rbindlist()


  source_list = lapply(l, \(x) x$source_data)


  if(n_boots > 0) {

    bootstraps = lapply(l, \(x) x$bootstraps) %>%
      rbindlist()

    structure(
      list(
        cic_data = tbls,
        source_data = source_list,
        training_data = cic_data$training_data,
        bootstraps = bootstraps
      ),
      class = "cic"
    )
  } else {

    structure(
      list(
        cic_data = tbls,
        source_data = source_list,
        training_data = cic_data$training_data
      ),
      class = "cic"
    )
  }

}

#' Simulate from dirichlet distribution
#' @param n number of simulations
#' @param alpha alpha to be passed (assumed constant)
rdirichlet = function(n, alpha = 1) {
  x = rgamma(n, alpha)
  x / sum(x)
}

#' Reweight CDFs with random exponentials
#' @param cdf cumulative distribution to randomly weight
#' @details Used for bootstrapping but nothing else
randweight_cdf = function(cdf) {
  x = get_x(cdf)
  w = rdirichlet(length(x))
  wecdf(x, w, normalize_weights = T)
}

#' Get the training data for an empirical CDF
#' @param e empirical distribution function to get
#' stuff from
get_x = function(e) {
  get("x", envir = environment(e))
}

#' Helper function, gets a cdf from a list given an index
#' @param cdflist list of cdfs from construct_cdfs
#' @param i person index
#' @param t time index
#' @param g group index
get_cdf = function(cdflist, i, t) {
  id = paste0(i, ".", t)
  cdflist[[id]]
}

#' Compute Changes in Changes model
#' @param i person index
#' @param t time index
#' @param y outcome index, y = 1 assumed to be what we proxy for
#' @param g optional group index
#' @param t_window time-window to restrict potential counterfactuals
#' @export
cic = function(i, t, y, x, g = NULL, t_window = NULL,
               n_boots = 50) {
  if(is.null(g)) {
    g = rep(1, length(i))
  }
  if(is.null(t_window)) {
    t_window = NA
  }

  data = build_cic_data(i, t, g, y, x, t_window)
  l = cic_from_cicdata(
    data, n_boots
  )

  structure(l, class = "cic")
}


#' Print function for cic
#' @param x fitted cic model
#' @param ... other arguments to print
#' @export
print.cic = function(x, ...) {
  print(x$cic_data, ...)
}






