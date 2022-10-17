# Give all combinations of m whole numbers such that 1, 2, and degree are always
# included. It is an error to have degree < 2, m < 2, or m < degree; or, when
# m == 2, any degree other than 2.
combx = function(degree, m)
{
    stopifnot("m must be less than or equal to degree" = degree >= m,
        "degree must be 2 or greater" = degree >= 2,
        "m must be 2 or greater" = m >= 2,
        "degree cannot be greater than 2 when m is 2" = !(degree > 2 && m == 2));

    if (degree == 2) {
        return (matrix(1:2, nrow = 2))
    } else {
        # Double negation below is to short-circuit the behaviour of combn
        # in which the first argument being a positive integer n means 1:n
        # instead of just n.
        return (unname(rbind(1, 2, -combn(-(3:(degree-1)), m - 3), degree)));
    }
}

# Give the number of columns that would be returned by a call to combx
combxn = function(degree, m)
{
    stopifnot("m must be less than or equal to degree" = degree >= m,
        "degree must be 2 or greater" = degree >= 2,
        "m must be 2 or greater" = m >= 2,
        "degree cannot be greater than 2 when m is 2" = !(degree > 2 && m == 2));

    if (degree == 2) {
        return (1)
    } else {
        return (choose(degree - 3, m - 3));
    }
}

# coefficients of symmetrical distributions
symmc = function(k) {
    ifelse(k <= 3, c(1,2,4)[k], (4 * ((k - 4) %/% 2) + ((k - 4) %% 2) + 6))
}

# Find exactly matching metalog distribution, if possible
search_exact = function(x, y, bounds, mu, n_highest, method)
{
    npts = length(x) + sum(!is.na(mu)); # number of CDF points / moments to fit
    mu0 = mu;

    if (npts < 2) {
        stop("Cannot fit fewer than 2 points.")
    }

    # TODO constraints on n_highest

    # Create 'master' matrix Z in order to solve for metalog coefficients.
    # If fitting to the mean has been requested, put into Z if possible.
    if (length(mu) > 0 && !is.na(mu[1]) && identical(bounds, c(-Inf, Inf))) {
        Z = metalog:::Y_matrix(length(x) + 1, y, mu[1], n_highest);
        x = c(x, mu[1]);
        mu[1] = NA;
    } else {
        Z = metalog:::Y_matrix(length(x), y, NULL, n_highest);
    }

    # Check for invalid requests relating to method
    method = match.arg(method, c("ols", "rf"));
    if (method == "ols" && any(!is.na(mu))) {
        stop("Cannot fit moments using OLS, except mean when distribution is unbounded.");
    }

    # Check all possible combinations of coefficients.
    # Each must contain the first two terms, for mean and variance.

    success = FALSE;
    # The outer loop is for checking whether any lower degree mlogs would suit -
    # this can happen if e.g. symmetric metalog distribution is requested with 3
    # points, because then a 2-degree mlog would be suitable.
    for (m in 2:npts) {
        # Quantile and moment indices
        qi = seq_len(min(m, length(x)));
        mi = which(cumsum(!is.na(mu)) <= m - length(qi));

        for (degree in m:ifelse(m == 2, 2, n_highest)) {
            coeffs = combx(degree, m);
            # coeffs = symmc(coeffs);
            for (p in 1:ncol(coeffs)) {
                result = try({
                    if (method == "ols") {
                        a = metalog:::OLS_fit_metalog(x[qi], NULL, NULL, m, Z[qi, coeffs[, p]]);
                    } else {
                        a = metalog:::RF_fit_metalog(x[qi], NULL, mu[mi], m, bounds, Z[qi, coeffs[, p]]);
                    }

                    # Assemble distribution and check whether distribution is valid...
                    ml = mlog(a, "coeff", bounds = bounds);
                    # TODO doing this seems weird - I think it is confusing that
                    # the x values get transformed in mlog.
                    mlu = mlog(a, "coeff", bounds = c(-Inf, Inf));

                    # And whether fit is exact
                    target_xmu = c(x, mu);
                    target_xmu = target_xmu[!is.na(target_xmu)];
                    realized_xmu = c(qmlog(y, mlu),
                        metalog:::standard_moments(ml, which(!is.na(mu0))));
                    if (mean(abs(target_xmu - realized_xmu)) > 1e-8) {
                        stop("Inexact fit")
                    }
                }, silent = TRUE);

                if (!is(result, "try-error")) {
                    success = TRUE;
                    break;
                }
            }
            if (success) {
                break;
            }
        }
        if (success) {
            break;
        }
    }

    if (!success) {
        stop("Could not find appropriate coefficients.")
    }

    ml
}
