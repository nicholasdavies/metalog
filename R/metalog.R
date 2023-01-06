#' Build a metalog distribution
#'
#' Fits a metalog distribution to specified quantiles or data.
#'
#' Not all sets of coefficients \ifelse{html}{\out{<i>a<sub>i</sub></i>}}{\eqn{a_i}}
#' result in a feasible metalog distribution. Metalog distributions are designed
#' to be very flexible so that they can approximate a wide variety of different
#' probability distributions. But the cost of this flexibility is that it is
#' possible to specify a set of coefficients \ifelse{html}{\out{<i>a<sub>i</sub></i>}}{\eqn{a_i}}
#' resulting in a PDF that goes negative in some areas, or equivalently, with a
#' CDF that is decreasing in some areas. This makes for an invalid probability
#' distribution, which means that the coefficients are not "feasible" in the
#' terminology of Keelin (2016).
#'
#' In particular, the OLS (ordinary least squares) method of fitting metalog
#' distributions, which attempts to fit the specified quantiles as closely as
#' possible (or exactly in the case that the number of terms \code{n} is
#' equal to the number of specified quantiles), is not guaranteed to result in a
#' feasible distribution, though it often does, especially when the number of
#' terms \code{n} is lower. The LP (linear programming) method tries to ensure
#' that the resulting distribution is feasible, but it does this by compromising
#' on the exactness of the distribution's fit to the specified data, and the
#' success of this method cannot be guaranteed in general.
#'
#' See Keelin (2016) for more details on feasibility.
#'
#' @param x A numeric vector of length 2 or greater. Interpretation depends on
#'  the value of \code{y}.
#' @param y One of:
#' \itemize{
#'  \item{A numeric vector the same length as \code{x}: treated as probabilities
#'  corresponding to each value in \code{x}, defining the cumulative
#'  distribution function to be fitted. Both \code{x} and \code{y} must be
#'  specified in increasing order.}
#'  \item{A single number: treated as a confidence level, e.g. with 0.95 meaning
#'  a 95\% confidence interval. In this case, \code{x} should be either 2 or 3
#'  elements long, specifying either \code{c(ci_lower, ci_upper)} or
#'  \code{c(ci_lower, median, ci_upper)}.}
#'  \item{\code{"data"}: values in \code{x} are treated as empirical data
#'  (i.e. draws from an underlying distribution) to which the metalog
#'  distribution is fitted.}
#'  \item{\code{"bin"}: like \code{"data"}, but the data are binned first using
#'  base R's \code{hist} function. This is faster for very large data sets, but
#'  using \code{"data"} seems to generally give better results.}
#'  \item{\code{"coeff"}: the \code{x}s directly specify the coefficients
#'  \ifelse{html}{\out{<i>a<sub>i</sub></i>}}{\eqn{a_i}} of the metalog
#'  distribution.}
#' }
#' Character values for \code{y} can be abbreviated.
#' @param bounds Lower and upper bounds for the distribution. Can be one of:
#' \itemize{
#'  \item{\code{c(-Inf, Inf)}: Unbounded distribution, with support over the
#'  entire real line. The default.}
#'  \item{\code{c(bl, Inf)}: Semi-bounded distribution, with support over
#'  (\code{bl}, \eqn{\infty}).}
#'  \item{\code{c(-Inf, bu)}: Semi-bounded distribution, with support over
#'  (\eqn{-\infty}, \code{bu}).}
#'  \item{\code{c(bl, bu)}: Bounded distribution, with support over
#'  (\code{bl}, \code{bu}).}
#' }
#' @param n Number of terms in the metalog distribution. Must be 2 or greater. A
#' metalog distribution with \eqn{n} terms can, in principle, fit \eqn{n}
#' quantiles exactly.
#' @param method Which method to use in fitting the metalog distribution (when
#' \code{y} is something other than \code{"coeff"}). One of:
#' \itemize{
#'  \item{\code{"any"}: Try OLS first; if this fails to identify a feasible
#'  distribution, then try LP.}
#'  \item{\code{"ols"}: Ordinary least squares. This is the faster of the two
#'  methods, but will fail if the resulting metalog distribution is not
#'  feasible (see details).}
#'  \item{\code{"lp"}: Linear programming, using the R package \code{lpSolve}.
#'  This attempts to find a feasible metalog distribution that approximates the
#'  provided data as closely as possible.}
#'  \item{\code{"nl"}: Nonlinear optimization, using the R package \code{nloptr}.
#'  This attempts to find a feasible metalog distribution that approximates the
#'  provided data as closely as possible.}
#'  \item{\code{"cvx"}: Convex optimization, using the R package \code{CVXR}.
#'  This attempts to find a feasible metalog distribution that approximates the
#'  provided data as closely as possible.}
#'  \item{\code{"rf"}: Root finding, using the R package \code{rootSolve}.
#'  This attempts to find a feasible metalog distribution that approximates the
#'  provided data as closely as possible.}
#'  \item{\code{"exact"}: Exact search. First, this will try to use the
#'  symmetric-percentile triplet method if that is appropriate. Next, this will
#'  try to find a metalog distribution fitting the supplied quantiles and
#'  moments exactly, starting with an \code{n}-term metalog to an
#'  \code{n + 5}-term metalog, using ordinary least squares. If that fails, it
#'  will try again using root finding.}
#' }
#' Character values for \code{method} can be abbreviated.
#' @param mu If non-NULL, a vector specifying the moments of the distribution to
#' be fitted. The first entry is the mean, the second is the standard deviation,
#' and the \eqn{k}th entry (for \eqn{k > 2}) is the \eqn{k}th standardized moment.
#' If a moment is \code{NA}, then it is 'skipped'; e.g. mu = c(NA, 1) specifies
#' a standard deviation of 1 but leaves the mean unspecified.
#' @param refine Whether to refine the cached version of the metalog quantile
#' function and PDF. This is enabled in order to slightly speed up calls to
#' \link{pmlog} and \link{dmlog}, although testing has revealed that this
#' sometimes has the opposite effect, and this parameter may be removed in the
#' future.
#'
#' @return An \code{mlog} object with elements:
#' \item{a}{The coefficients \ifelse{html}{\out{<i>a<sub>i</sub></i>}}{\eqn{a_i}}
#' of the fitted metalog distribution.}
#' \item{b}{The bounds of the fitted metalog distribution; same as the parameter
#' \code{bounds} supplied to the function.}
#' \item{cache}{A \code{data.frame} containing a cached version of the metalog
#' quantile function \eqn{M} and PDF \eqn{m}, as specified in eqs. 6 and 9 in
#' Keelin (2016).}
#' \item{method}{A character string describing the method that was actually
#' used to determine the coefficients \ifelse{html}{\out{<i>a<sub>i</sub></i>}}{\eqn{a_i}},
#' either \code{"ordinary least squares"}, \code{"linear programming"}, or
#' \code{"specified coefficients"}.}
#'
#' @references T.W. Keelin (2016) The metalog distributions. \emph{Decision
#' Analysis} 13: 243-277. \url{https://doi.org/10.1287/deca.2016.0338}.
#'
#' @examples
#' # UNBOUNDED DISTRIBUTIONS
#'
#' # These are all equivalent:
#' ml1 = mlog(c(20, 30), c(0.05, 0.95))
#' ml2 = mlog(c(20, 25, 30), c(0.05, 0.5, 0.95))
#' ml3 = mlog(c(20, 30), 0.9)
#' ml4 = mlog(c(20, 25, 30), 0.9)
#'
#' # A skewed distribution:
#' ml5 = mlog(c(20, 28, 30), 0.9)
#'
#' # An empirically fitted distribution:
#' ml6 = mlog(quakes$lat)
#'
#' @export
mlog = function(x, y = "data", bounds = c(-Inf, Inf), n = NULL,
    mu = NULL, method = "auto", refine = TRUE, check_validity = TRUE)
{
    # Requirements: x and mu must be numeric
    if ((!is.null(x) && !is.numeric(x)) || (!is.null(mu) && !is.numeric(mu))) {
        stop("x and mu must be numeric.")
    }

    # Requirements on bounds
    if (!is.numeric(bounds) || length(bounds) != 2 || !all(is.infinite(bounds) | is.finite(bounds)) || bounds[1] >= bounds[2]) {
        stop("bounds must be a vector with 2 ordered numeric entries.")
    }

    # Extract some arguments
    if (is.character(y)) {
        y = match.arg(tolower(y), c("bin", "data", "coeff", "moments"));
    }
    bl = bounds[1];
    bu = bounds[2];
    nmu = sum(!is.na(mu));
    npts = length(x) + nmu;

    # Requirements: x and mu together must have at least 2 points
    if (npts < 2) {
        stop("Must supply at least 2 points to fit to.")
    }

    # Transform quantiles for fitting
    x0 = x;
    if (!identical(y, "coeff")) {
        x = metalog:::mlog_trans(x, bl, bu);
    }

    # Processing of x and y
    if (is.character(y)) {
        if (identical(y, "bin")) {
            # TODO change when bounds implemented to ensure bins inside bounds.
            breaks = hist(x, right = FALSE, plot = FALSE)$breaks;
            binwidth = breaks[2] - breaks[1];
            if (any(x == breaks[length(breaks)])) {
                breaks = c(breaks, breaks[length(breaks)] + binwidth);
            }
            histogram = hist(x, breaks = breaks, include.lowest = FALSE, right = FALSE, plot = FALSE);

            x = as.numeric(histogram$breaks);
            y = c(0, cumsum(histogram$density)) / sum(histogram$density);

            x = x[y > 0 & y < 1];
            y = y[y > 0 & y < 1];
            npts = length(x) + sum(!is.na(mu));

            # Sensible default for n
            if (is.null(n)) {
                n = min(npts, 6);
            }

            # Sensible default for method
            if (method == "auto") {
                method = "any";
            }
        } else if (identical(y, "data")) {
            tab = metalog:::xtabulate(x);
            x = tab$x;
            y = cumsum(tab$n) / (sum(tab$n) + 1);

            # Sensible default for n
            if (is.null(n)) {
                n = min(npts, 6);
            }

            # Sensible default for method
            if (method == "auto") {
                method = "any";
            }
        } else if (identical(y, "coeff")) {
            if (sum(!is.na(mu)) > 0) {
                stop("cannot specify moments when y = \"coeff\"")
            }
            a = x;
            if (is.null(n)) {
                n = length(x);
            } else if (n != length(x)) {
                stop("supplied ", length(x), " coefficients, but specified n = ", n);
            }
        } else if (identical(y, "moments")) {
            # TODO this option should be auto-detected
            y = NULL;

            # Sensible default for n
            if (is.null(n)) {
                n = npts;
            }

            # Sensible default for method
            if (method == "auto") {
                method = "exact";
            }
        }
    } else if (is.numeric(y)) {
        # Requirement: y must be between 0 and 1
        if (any(y <= 0 | y >= 1)) {
            stop("y must be strictly between 0 and 1.")
        }

        # One y value: represents confidence level
        if (length(y) == 1) {
            alpha = 1 - y;
            if (length(x) == 2) {
                y = c((1 - y) / 2, 1 - (1 - y) / 2);
            } else if (length(x) == 3) {
                y = c((1 - y) / 2, 0.5, 1 - (1 - y) / 2);
            } else {
                stop("cannot interpret y as confidence interval unless x has 2 or 3 elements")
            }
        } else if (length(x) != length(y)) {
            # Requirements: x and y must be the same length (if y not length 1)
            stop("x and y must be the same length.")
        }

        # Requirement: x and y must both be increasing
        if (is.unsorted(x) || is.unsorted(y)) {
            stop("x and y must be sorted.")
        }

        # Sensible default for n
        if (is.null(n)) {
            n = npts;
        }

        # Sensible default for method
        if (method == "auto") {
            method = "exact";
        }
    } else {
        stop("cannot interpret y")
    }

    # Requirements on n
    n = as.integer(n);
    if (n < 2) {
        stop("n must be 2 or greater.")
    }

    # m: length of CDF vector
    m = length(x);

    if (!identical(y, "coeff")) {
        a = NULL;
        method = match.arg(tolower(method), c("any", "ols", "lp", "nl", "cvx", "rf", "exact"));

        if (identical(method, "any")) {
            a = OLS_fit_metalog(x, y, mu, n);
            method_str = "ordinary least squares";
            # TODO this needs to check validity.
            if (is.null(a)) {
                a = LP_fit_metalog(x, y, mu, n);
                method_str = "linear programming";
            }
        } else if (identical(method, "ols")) {
            a = OLS_fit_metalog(x, y, mu, n);
            method_str = "ordinary least squares";
        } else if (identical(method, "lp")) {
            a = LP_fit_metalog(x, y, mu, n);
            method_str = "linear programming";
        } else if (identical(method, "nl")) {
            a = NL_fit_metalog(x, y, mu, n, bounds);
            method_str = "nonlinear optimization";
        } else if (identical(method, "cvx")) {
            a = CVX_fit_metalog(x, y, mu, n, bounds);
            method_str = "convex optimization";
        } else if (identical(method, "rf")) {
            a = RF_fit_metalog(x, y, mu, n, bounds);
            method_str = "root finding";
        } else if (identical(method, "exact")) {
            if (nmu == 0 && length(y) == 3 && y[2] == 0.5 && 1 - y[1] == y[3]) {
                a = metalog:::SPT_a(x0[1], x0[2], x0[3], y[1], bounds);
                method_str = "symmetric-percentile triplet";
            } else if (all(is.na(mu[-1])) && identical(bounds, c(-Inf, Inf))) {
                # TODO n + 5 customizable
                # TODO this has already done feasibility constraints, and ignores check_validity
                a = search_exact(x, y, bounds, mu = mu, n_highest = n + 5, method = "ols")$a;
                method_str = "exact search by ordinary least squares";
            } else {
                # TODO n + 5 customizable
                # TODO this has already done feasibility constraints, and ignores check_validity
                a = search_exact(x, y, bounds, mu = mu, n_highest = n + 5, method = "rf")$a;
                method_str = "exact search by root finding";
            }
        }

        n = length(a); # Fitting methods may adjust the value of n
    } else {
        method_str = "specified coefficients";
        y = 0.5; # for minp and maxp below in case of y = "coeff"
    }

    # Probability values for cached quantile function and PDF.
    # TODO make 1e-4 a parameter of the function
    minp = min(1e-4, y[1]);
    maxp = max(1 - 1e-4, y[length(y)]);
    yq = seq(minp, maxp, length.out = 100);

    # Create cached quantile function (M)
    while (TRUE) {
        M = metalog:::metalog_M(yq, a, n, bl, bu);

        # The way yq has been initialised, the cached points are spaced evenly
        # in probability space. We can use the cached quantile function to
        # interpolate new values of yq such that the cached points will be
        # spaced more evenly in quantile space, which then requires the cached
        # quantile function to be recalculated.
        if (refine) {
            yq = plogis(approx(M, qlogis(yq), n = 100)$y);
            refine = FALSE;
            next;
        } else {
            break;
        }
    }

    # Create cached PDF (m)
    m = metalog:::metalog_m(yq, a, n, bl, bu);

    # Check validity
    if (any(m < 0) && check_validity) {
        stop("Invalid metalog distribution.")
    }

    # Create metalog object
    ml = list(
        a = a,
        b = bounds,
        cache = data.frame(
            y = yq,
            M = M,
            m = m
        ),
        method = method_str
    );
    class(ml) = c("mlog", class(ml));

    return (ml)
}

mlog_trans = function(x, bl, bu)
{
    if (is.finite(bl)) {
        if (is.finite(bu)) {
            # c(bl, bu)
            return (log((x - bl) / (bu - x)));
        } else {
            # c(bl, Inf)
            return (log(x - bl));
        }
    } else if (is.finite(bu)) {
        # c(-Inf, bu)
        return (-log(bu - x));
    }
    return (x);
}
