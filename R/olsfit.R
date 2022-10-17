OLS_fit_metalog = function(x, y, mean, n, Y = NULL, LOOCV = FALSE)
{
    x = c(x, mean); # CDF points / moments to fit
    m = length(x);  # number of CDF points / moments to fit

    if (n > m) {
        # TODO this should say the number of unique points for y = "data" or bins for y = "bin"
        warning("n = ", n, " was requested, but this is greater than the number",
            " of points, m = ", m, ". Setting n = ", m, ".")
        n = m;
    }

    # Create matrix Y (m * n) in order to solve for metalog coefficients,
    # unless it is being supplied.
    if (is.null(Y)) {
        Y = metalog:::Y_matrix(m, y, mean, n);
    }

    # Solve for metalog coefficients: eq. 7, Keelin 2016
    # Attempts to recover from errors in solve due to the system being
    # computationally singular by decreasing the tolerance.

    # First, try solving with the default tolerance.
    tol = .Machine$double.eps;
    result = try(
        if (m == n) {
            a = solve(Y, x, tol = tol);
        } else {
            a = solve(crossprod(Y), crossprod(Y, x), tol = tol)[,1];
        },
        silent = TRUE);

    # If this hasn't worked, try lowering the tolerance. What this really
    # should do is verify that attr(try(...), "condition") is of the form
    # "system is computationally singular: reciprocal condition number = XXX".
    # But this requires different behaviour depending on the current language,
    # which may not be English, and moreover the error message could change in
    # future versions of R. Oh well. If this second attempt fails, the error
    # is propagated to the user.
    if (is(result, "try-error")) {
        warning(attr(result, "condition"),
            "Trying again with lower numerical tolerance.");

        result = try(
            if (m == n) {
                tol = rcond(Y) / 2;
                a = solve(Y, x, tol = tol);
            } else {
                cpY = crossprod(Y);
                tol = rcond(cpY) / 2;
                a = solve(cpY, crossprod(Y, x), tol = tol)[,1];
            },
            silent = TRUE);

        if (is(result, "try-error")) {
            warning("Second attempt at OLS with lower tolerance failed: ",
                attr(result, "condition"));
            return (NULL);
        }
    }

    # Calculate the leave-one-out cross-validation statistic if requested
    if (LOOCV) {
        if (m == n) {
            stop("Cannot calculate CV when m == n.")
        } else {
            cpY = crossprod(Y);
            tol = rcond(cpY) / 2;
            H = Y %*% solve(cpY, tol = t) %*% t(Y);
            xH = H %*% x; # predicted x given model
            CV = mean(((x - xH) / (1 - diag(H)))^2);

            return (list(a = a, CV = CV))
        }
    }

    # Calculate AIC if requested

    return (a)
}
