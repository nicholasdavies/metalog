
combx(3, 2)
combx(4, 2)
combx(5, 2)
combx(6, 2)

for (m in 3:12)
{
    for (degree in m:12)
    {
        ncol1 = combxn(degree, m)
        ncol2 = ncol(combx(degree, m))
        print(paste(ncol1, ncol2))

        if (ncol1 != ncol2) {
            stop(degree, m)
        }
    }
}

combx(2, 2)

x = c(0.05, 0.80, 0.95)
y = c(0.05, 0.50, 0.95)
bounds = c(-Inf, Inf)
mu = NULL
n_highest = 12


x = c(0.05, 0.50, 0.95)
y = c(0.05, 0.50, 0.95)
bounds = c(-Inf, Inf)
mu = 0.50
n_highest = 12
method = "rf"

x = c(-1.1, 1.1)
y = c(0.05, 0.95)
bounds = c(-Inf, Inf)
mu = c(NA, NA)
n_highest = 12
method = "rf"
ml = something(x, y, bounds, mu, n_highest, method)
plot(ml)

something = function(x, y, bounds, mu, n_highest, method)
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
            for (p in 1:ncol(coeffs)) {
                result = try({
                    if (method == "ols") {
                        a = metalog:::OLS_fit_metalog(x[qi], NULL, NULL, m, Z[qi, coeffs[, p]]);
                    } else {
                        a = metalog:::RF_fit_metalog(x[qi], NULL, mu[mi], m, bounds, Z[qi, coeffs[, p]]);
                    }

                    # Assemble distribution and check whether distribution is valid...
                    ml = mlog(a, "coeff");

                    # And whether fit is exact
                    target_xmu = c(x, mu);
                    target_xmu = target_xmu[!is.na(target_xmu)];
                    realized_xmu = c(qmlog(y, ml),
                        metalog:::standard_moments(ml, which(!is.na(mu0))));
                    if (mean(abs(target_xmu - realized_xmu)) > 1e-8) {
                        ml = NULL;
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







x = c(-1, 0)
y = c(0.05, 0.5)
metalog:::RF_fit_metalog(x, y, NULL, 2, c(-Inf, Inf))



mlog(c(244653.27,1134301.532), c(0.025, 0.975), method = "exact")


new_ALRI_lo = 244653.27
new_ALRI_hi = 1134301.532
new_ALRI = 489306.543

mlog(c(new_ALRI_lo, new_ALRI_hi), 0.95, mean = new_ALRI, c(0, Inf), method = "exact", refine = FALSE)
