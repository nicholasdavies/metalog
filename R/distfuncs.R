# Metalog quantile function
metalog_M = function(y, a, n, bl, bu)
{
    # Eq. 6, Keelin 2016:
    log_odds_y = log(y / (1 - y));
    ym0.5 = y - 0.5;
    M = a[1] + a[2] * log_odds_y;
    if (n >= 3) {
        M = M + a[3] * ym0.5 * log_odds_y;
    }
    if (n >= 4) {
        M = M + a[4] * ym0.5;
    }
    if (n >= 5) {
        k = 2;
        for (j in seq(5, n, by = 2)) {
            factor = ym0.5 ^ k;
            M = M + a[j] * factor;
            if (j + 1 <= n) {
                M = M + a[j + 1] * factor * log_odds_y;
            }
            k = k + 1;
        }
    }

    # Adjust for bounds
    # Keelin 2016, eqs. 11 and 14
    if (is.finite(bl)) {
        if (is.finite(bu)) {
            eM = exp(M);
            M = (bl + bu*eM) / (1 + eM);
        } else {
            M = bl + exp(M);
        }
    } else if (is.finite(bu)) {
        M = bu - exp(-M);
    }

    M[y == 0] = bl;
    M[y == 1] = bu;

    return (M)
}

# Metalog quantile function on complement of y
metalog_Mc = function(yc, a, n, bl, bu)
{
    # Eq. 6, Keelin 2016:
    # TODO should I use log1p forms?
    # log_odds_y = log1p(-yc) - log(yc);
    log_odds_y = log((1 - yc) / yc);
    ym0.5 = 0.5 - yc;
    M = a[1] + a[2] * log_odds_y;
    if (n >= 3) {
        M = M + a[3] * ym0.5 * log_odds_y;
    }
    if (n >= 4) {
        M = M + a[4] * ym0.5;
    }
    if (n >= 5) {
        k = 2;
        for (j in seq(5, n, by = 2)) {
            factor = ym0.5 ^ k;
            M = M + a[j] * factor;
            if (j + 1 <= n) {
                M = M + a[j + 1] * factor * log_odds_y;
            }
            k = k + 1;
        }
    }

    # Adjust for bounds
    # Keelin 2016, eqs. 11 and 14
    if (is.finite(bl)) {
        if (is.finite(bu)) {
            eM = exp(M);
            M = (bl + bu*eM) / (1 + eM);
        } else {
            M = bl + exp(M);
        }
    } else if (is.finite(bu)) {
        M = bu - exp(-M);
    }

    M[yc == 1] = bl;
    M[yc == 0] = bu;

    return (M)
}

# Metalog PDF
metalog_m = function(y, a, n, bl, bu)
{
    # Eq. 9, Keelin 2016:
    log_odds_y = log(y / (1 - y));
    ym0.5 = y - 0.5;
    yc = y * (1 - y);
    if (n == 2) {
        m = a[2] / yc;
    } else if (n == 3) {
        m = a[2] / yc + a[3] * (ym0.5/yc + log_odds_y);
    } else if (n >= 4) {
        m = a[2] / yc + a[3] * (ym0.5/yc + log_odds_y) + a[4];
    }
    if (n >= 5) {
        k = 1;
        for (j in seq(5, n, by = 2)) {
            factor = ym0.5 ^ k;
            m = m + a[j] * (j-1)/2 * factor;
            if (j + 1 <= n) {
                m = m + a[j+1] * (factor * ym0.5 / yc + (k+1) * factor * log_odds_y);
            }
            k = k + 1;
        }
    }
    m = 1 / m;

    # Adjust for bounds: Keelin 2016, eqs. 13 and 15
    if (is.finite(bl)) {
        if (is.finite(bu)) {
            eM = exp(metalog_M(y, a, n, -Inf, Inf));
            m = m * (1 + eM)^2 / ((bu - bl) * eM);
        } else {
            m = m * exp(-metalog_M(y, a, n, -Inf, Inf));
        }
    } else if (is.finite(bu)) {
        m = m * exp(metalog_M(y, a, n, -Inf, Inf));
    }

    return (m)
}

#' The metalog distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the metalog distribution specified by \code{ml}.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param ml Specification of the metalog distribution. Normally this is
#' obtained via a call to \code{\link{mlog}}.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]},
#' otherwise \eqn{P[X > x]}.
#'
#' @return \code{dmlog} gives the density, \code{pmlog} gives the distribution
#' function, \code{qmlog} gives the quantile function, and \code{rmlog}
#' generates random deviates.
#'
#' The length of the result is determined by \code{n} for \code{rmlog}, and by
#' the length of the first argument (\code{x}, \code{q}, or \code{p}) for
#' the other functions.
#'
#' @references T.W. Keelin (2016) The metalog distributions. \emph{Decision
#' Analysis} 13: 243-277. \url{https://doi.org/10.1287/deca.2016.0338}.
#'
#' @examples
#' ml = mlog(quakes$lat);
#'
#' x = seq(min(quakes$lat), max(quakes$lat), length.out = 100);
#' d = dmlog(x, ml);
#' p = pmlog(x, ml);
#'
#' plot(x, d);
#' hist(rmlog(1000, ml), breaks = 100);
#'
#' plot(x, p);
#' plot(x, qmlog(p, ml));
#'
#' @export
dmlog = function(x, ml, log = FALSE)
{
    d = bmlog(pmlog(x, ml), ml);
    if (log) {
        return (log(d))
    } else {
        return (d)
    }
}

#' @rdname dmlog
#' @export
pmlog = function(q, ml, lower.tail = TRUE, log.p = FALSE)
{
    # There is no closed-form CDF for metalog distributions, so we need to find
    # p such that qmlog(p) == q using numerical methods; here using Newton-Raphson.

    # Begin with a guess for p by interpolating the cached quantile function.
    p = approx(ml$cache$M, ml$cache$y, q, rule = 2)$y;

    # Iterate with Newton-Raphson until guesses stop improving, guess is within
    # tolerance or maximum iterations are reached.
    max_iter = 20;
    iter = 0;
    prev_error = Inf;
    error = max(abs(q - qmlog(p, ml)));

    while (error > .Machine$double.eps & iter < max_iter & error < prev_error) {
        b = bmlog(p, ml);
        p = pmax(1e-16, pmin(p + (q - qmlog(p, ml)) * b, 1 - 1e-16));
        prev_error = error;
        error = max(abs(q - qmlog(p, ml)));
        iter = iter + 1;
    }

    # Respect bounds
    p[q <= ml$b[1]] = 0;
    p[q >= ml$b[2]] = 1;

    # Complement probabilities if requested
    if (!lower.tail) {
        p = 1 - p;
    }
    # Take logarithm of probabilities if requested
    if (log.p) {
        p = log(p);
    }

    return (p)
}

#' @rdname dmlog
#' @export
qmlog = function(p, ml, lower.tail = TRUE, log.p = FALSE)
{
    # Transform inputs if needed
    if (log.p) {
        p = exp(p);
    }
    if (!lower.tail) {
        p = 1 - p;
    }

    # Check input is valid
    NaNs = p < 0 | p > 1;
    p = pmax(0, pmin(p, 1));

    # Calculate quantile function
    M = metalog_M(p, ml$a, length(ml$a), ml$b[1], ml$b[2]);

    # Handle edge cases (p <= 0 or p >= 1)
    M[p == 0] = ml$b[1];
    M[p == 1] = ml$b[2];
    M[NaNs] = NaN;

    # Emit warning if any NaNs produced
    if (any(NaNs)) {
        warning("NaNs produced")
    }

    return (M)
}

# inverse derivative of qmlog = PDF of metalog
bmlog = function(p, ml)
{
    # Calculate PDF
    metalog_m(p, ml$a, length(ml$a), ml$b[1], ml$b[2]);
}

#' @rdname dmlog
#' @export
rmlog = function(n, ml)
{
    U = runif(n);
    return (qmlog(U, ml))
}
