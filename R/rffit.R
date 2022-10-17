# Function for root finding
multiroot_f = function(a, ml_n, ml_bounds,
    target_mu, target_x, target_y, testpts, Z)
{
    # Calculate current moments
    if (is.null(target_mu)) {
        current_mu = NULL;
    } else {
        current_mu = standard_moments(
            list(a = a, b = ml_bounds),
            which(!is.na(target_mu))
        );
    }

    # Calculate current quantiles
    if (is.null(Z)) {
        # no Z matrix supplied: Note target_x is supplied in terms of the
        # underlying unbounded distribution, so we use -Inf and Inf as bounds here.
        current_x = metalog:::metalog_M(target_y, a = a, n = ml_n, bl = -Inf, bu = Inf);
    } else {
        # Use supplied Z matrix
        current_x = (Z %*% a)[, 1];
    }

    # Function value
    c(target_mu[!is.na(target_mu)], target_x) - c(current_mu, current_x)
}

# Root finding
RF_fit_metalog = function(x, y, mu, n, bounds, Z = NULL)
{
    # Find root
    atol = 1e-8;
    start = rep(0, n);
    guess_a1 = mean(c(x, mlog_trans(mu[1], bounds[1], bounds[2])), na.rm = TRUE);
    if (is.finite(guess_a1)) {
        start[1] = guess_a1;
    }

    res = rootSolve::multiroot(multiroot_f, start = start, maxiter = 1000,
        rtol = 0, atol = 1e-8, ctol = 0,
        ml_n = n, ml_bounds = bounds,
        target_mu = mu, target_x = x, target_y = y,
        testpts = testpts, Z = Z);

    if (mean(abs(res$estim.precis)) > atol) {
        stop("Could not find root.")
    }

    # TODO error checking etc

    res$root
}
