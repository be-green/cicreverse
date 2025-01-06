
#' Check for NA values
#' @param x vector to check
#' @param na.rm whether to remove NAs
check_na = function(x, na.rm) {
  if(na.rm) {
    x[!is.na(x)]
  } else {
    if(anyNA(x)) {
      stop("x must not contain NA values")
    }
    x
  }
}

#' Method for the quantile function for weighted empriical CDFs
#' @param x empirical cdf function of class weighted_cdf
#' @param probs vector of probabilities to evaluate quantile function for
#' @param na.rm whether to remove NAs
#' @param ... other things to pass to weighted quantile
#' @export
#' @details This just returns the inverse cdf instead of all of those
#' fancy quantile versions.
quantile.wecdf = function (x, probs, ...) {
  data = environment(x)$x
  weights = environment(x)$weights

  # already sorted by wecdf
  cum_weights = cumsum(weights)
  stepfun(cum_weights, c(data, data[length(data)]), right = TRUE)(probs)
}


#' Weighted empirical cdf function
#' @param x vector to calculate the weighted empirical CDF for
#' @param weights vector of weights, must not contain NAs
#' @param na.rm whether to remove NAs from x and associated weight vector
#' @param normalize_weights whether to normalize weights to sum to 1
#' @details This is _heavily_ based on the code from the
#' ggdist package. Thanks to Matthew Kay and Brenton Wiernik
#' for writing a lovely package. I would have imported it
#' but I didn't want a whole dependency and wanted some additional
#' checks.
wecdf = function(x, weights = NULL, na.rm = FALSE,
                 normalize_weights = TRUE) {
  if (is.null(weights)) {
    weights = rep(1, length(x))
  }
  if(abs(sum(weights) - 1) > 1e-10) {
    if(normalize_weights) {
      weights = weights / sum(weights)
    } else {
      stop("Weights must sum to one if normalize_weights is FALSE")
    }
  }
  sort_order = order(x)
  x = x[sort_order]
  weights = weights[sort_order]
  p = cumsum(weights)/sum(weights)
  cdf = approxfun(x, p, yleft = 0, yright = 1,
                  ties = "ordered",
                  method = "constant")
  class(cdf) = c("wecdf", "stepfun", class(cdf))
  assign("weights", weights, envir = environment(cdf))
  attr(cdf, "call") = sys.call()
  cdf
}
