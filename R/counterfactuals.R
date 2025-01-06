#' Generate combined counterfactual
#' @param cic cic object retrned from function cic
#' @export
#' @details If we think of the set of counterfactuals, {F_itsj}
#' where a cdf is indexed by a comparison group j and time j,
#' as being samples from the set of CDFs that would lead to
#' each class being obtained under the same distribution of
#' unobservables, then E(F_y1(x)) is the quantile of the expected
#' CDF that leads to outcome y1 and we can approximate it
#' with the empirical analogue.
#' @importFrom data.table setnames
#' @importFrom data.table `:=`
#' @importFrom data.table setorder
#' @importFrom data.table shift
combined_cf = function(cic) {

  out_data = cic$cic_data
  training_data = cic$training_data
  bs_data = cic$bootstraps
  training_data = training_data[y == 1]

  cf_cdf_data = out_data[,.(CDF = wecdf(k_cic)(k_cic),
                            x = k_cic),
                         by = c("i","t","s","j")]

  data.table::setorder(cf_cdf_data, "i","t","s", "j", "x")

  # aggregate first by it pair before integrating over it
  group = c("i","t","s", "j")
  cf_cdf_data[, PDF := CDF - data.table::shift(CDF), by = group]
  cf_cdf_data[is.na(PDF), PDF := CDF]
  data.table::setorder(cf_cdf_data, x)

  # for mass points which come from winsorizing, for example
  cf_cdf_data = cf_cdf_data[,.(PDF = sum(PDF)),
                            by = c("x", "i", "t")]

  cf_cdf_data[,CDF_it := cumsum(PDF / sum(PDF)), by = c("i", "t")]
  cf_cdf_data[,CDF_it := cumsum(PDF / sum(PDF)), by = c("i", "t")]
  cf_cdf_data =
    cf_cdf_data[,.(CDF_it, x, type = "counterfactual"),
                by = c("i", "t")]

  group = c("i", "t")
  cf_cdf_data[, PDF := CDF_it - data.table::shift(CDF_it), by = group]
  cf_cdf_data[is.na(PDF), PDF := CDF_it]
  data.table::setorder(cf_cdf_data, x)
  cf_cdf_data = cf_cdf_data[,.(PDF = sum(PDF)),
                            by = c("x")]
  cf_cdf_data[,CDF_total := cumsum(PDF / sum(PDF))]
  cf_cdf_data[,CDF_total := cumsum(PDF / sum(PDF))]
  cf_cdf_data =
    cf_cdf_data[,.(CDF = CDF_total, x, type = "counterfactual")]


  actual_cdf_data = training_data[,.(CDF_it = wecdf(x)(x),
                                     x), by = c("i", "t")]
  data.table::setorder(actual_cdf_data, CDF_it, i, t, x)
  actual_cdf_data[, PDF := CDF_it - data.table::shift(CDF_it),
                  by = c("i", "t")]
  actual_cdf_data[is.na(PDF), PDF := CDF_it]
  data.table::setorder(actual_cdf_data, x)
  actual_cdf_data = actual_cdf_data[,.(PDF = sum(PDF)),
                                    by = c("x")]
  actual_cdf_data[,CDF_total := cumsum(PDF / sum(PDF))]
  actual_cdf_data[,CDF_total := cumsum(PDF / sum(PDF))]
  actual_cdf_data =
    actual_cdf_data[,.(CDF = CDF_total, x, type = "actual")]


  # for the bootstrap replicates now
  # need to refactor this function but it's fine for now i think
  if(!is.null(bs_data)) {

    bs_cf_cdf_data = bs_data[,.(CDF = wecdf(k_cic)(k_cic),
                              x = k_cic),
                           by = c("i","t","s","j","b")]

    data.table::setorder(bs_cf_cdf_data, "i","t","s","j","b","x")

    # aggregate first by it pair before integrating over it
    group = c("i","t","s","b","j")
    bs_cf_cdf_data[, PDF := CDF - data.table::shift(CDF), by = group]
    bs_cf_cdf_data[is.na(PDF), PDF := CDF]
    data.table::setorder(bs_cf_cdf_data, x)

    # for mass points which come from winsorizing, for example
    bs_cf_cdf_data = bs_cf_cdf_data[,.(PDF = sum(PDF)),
                              by = c("x", "i","b","t")]

    bs_cf_cdf_data[,CDF_it := cumsum(PDF / sum(PDF)), by = c("i", "t")]
    bs_cf_cdf_data[,CDF_it := cumsum(PDF / sum(PDF)), by = c("i", "t")]
    bs_cf_cdf_data =
      bs_cf_cdf_data[,.(CDF_it, x, type = "counterfactual"),
                  by = c("i", "b", "t")]

    group = c("i", "t", "b")
    bs_cf_cdf_data[, PDF := CDF_it - data.table::shift(CDF_it), by = group]
    bs_cf_cdf_data[is.na(PDF), PDF := CDF_it]
    data.table::setorder(bs_cf_cdf_data, x)
    bs_cf_cdf_data = bs_cf_cdf_data[,.(PDF = sum(PDF)),
                              by = c("x","b")]
    bs_cf_cdf_data[,CDF_total := cumsum(PDF / sum(PDF)), by = b]
    bs_cf_cdf_data[,CDF_total := cumsum(PDF / sum(PDF)), by = b]
    bs_cf_cdf_data =
      bs_cf_cdf_data[,.(CDF = CDF_total, x,
                        type = paste0("counterfactual_",b))]

    plot_data = rbind(actual_cdf_data, cf_cdf_data, bs_cf_cdf_data)
  } else {
    plot_data = rbind(actual_cdf_data, cf_cdf_data)
  }

  plot_data[]
}
