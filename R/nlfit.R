# Objective function for nonlinear optimization
NL_eval_f = function(a, ml_n, ml_bounds, target_mean, target_x, target_y, testpts)
{
    # Calculate current mean
    if (is.null(target_mean)) {
        current_mean = NULL;
    } else {
        current_mean = moment(list(a = a, b = ml_bounds));
    }

    # Calculate current quantiles - target_x is supplied in terms of the
    # underlying unbounded distribution, so we use -Inf and Inf as bounds here.
    current_x = metalog_M(target_y, a = a, n = ml_n, bl = -Inf, bu = Inf);

    # SSE for mean and quantiles
    # TODO weights?
    sum( (c(target_mean, target_x) - c(current_mean, current_x))^2 )
}

# Constraint function for nonlinear optimization
NL_eval_g = function(a, ml_n, ml_bounds, target_mean, target_x, target_y, testpts)
{
    mat = LP_feas_coeffs(n = ml_n, m = length(testpts$y),
        testpts$log_odds_y, testpts$ym0.5, testpts$ycinv);
    -mat %*% a
}

# Nonlinear fitting
NL_fit_metalog = function(x, y, mean, n, bounds)
{
    # Test points for feasibility
    p = c(0.00001, 0.0001, seq(0.001, 0.999, 0.001), 0.9999, 0.99999);
    k = length(p);

    # Precalculated quantities for checking constraints
    testpts = list(
        y = p,
        log_odds_y = log(p / (1 - p)),
        ym0.5 = p - 0.5,
        ycinv = 1 / (p * (1 - p))
    );

    # Set nloptr options
    x0 = rep(0, n);
    opts = list(algorithm = "NLOPT_LN_COBYLA", stopval = 1e-24, maxeval = 1e4);

    # Run optimization
    res = nloptr::nloptr(x0 = x0,
        eval_f = NL_eval_f, eval_g_ineq = NL_eval_g,
        opts = opts,
        ml_n = n, ml_bounds = bounds,
        target_mean = mean, target_x = x, target_y = y,
        testpts = testpts);

    # TODO add nloptr to imports
    # TODO error checking etc

    res$sol
}
