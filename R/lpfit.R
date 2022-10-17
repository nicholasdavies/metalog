pm = function(mat)
{
    pmmat = matrix(0, nrow = nrow(mat), ncol = 2 * ncol(mat));
    pmmat[1:nrow(mat), seq(1, ncol(pmmat), 2)] = mat;
    pmmat[1:nrow(mat), seq(2, ncol(pmmat), 2)] = -mat;
    pmmat
}

LP_fit_metalog = function(x, y, mean, n)
{
    x = c(x, mean); # CDF points / moments to fit
    m = length(x);  # number of CDF points / moments to fit

    # Test points for feasibility
    p = c(0.00001, 0.0001, seq(0.001, 0.999, 0.001), 0.9999, 0.99999);
    k = length(p);

    # Objective function
    cost = c(rep(1, 2*m), rep(0, 2*n));

    # Constraint coefficients
    log_odds_p = log(p / (1 - p));
    pm0.5 = p - 0.5;
    pcinv = 1 / (p * (1 - p));

    A = rbind(
        cbind(pm(diag(m)),       pm(Y_matrix(m, y, mean, n))),
        cbind(matrix(0, k, 2*m), pm(LP_feas_coeffs(n, k, log_odds_p, pm0.5, pcinv)))
    );

    # Constraint directions
    dir = c(rep("==", m), rep(">=", k));

    # Right-hand side
    diff_error = 0.001;
    rhs = c(x, rep(diff_error, k));

    # Build lpSolve model
    lp_sol = lpSolve::lp("min", cost, A, dir, rhs);
    if (lp_sol$status == 0) { # success
        coeffs = lp_sol$solution[(2 * m + 1):(length(lp_sol$solution))];
        coeffs = coeffs[seq(1, 2 * n - 1, 2)] - coeffs[seq(2, 2 * n, 2)];
    } else { # failure
        coeffs = NULL;
    }

    return (coeffs)
}

# An attempt to use an alternative version to the above LP solver which allows
# negative coefficients to be estimated directly, rather than using the trick
# that derives negative coefficients from the subtraction of two positive-bounded
# decision variables; this turned out to be slower to solve, at least for some
# test cases.
# LP_fit_metalog = function(x, y, n)
# {
#     m = length(x); # number of CDF points
#
#     # Test points for feasibility
#     k = 999;
#     p = seq(1e-3, 1 - 1e-3, length.out = k);
#
#     # Objective function
#     cost = c(rep(0, n), rep(1, m)); # ND attempt
#
#     # Constraints
#     log_odds_y = log(y / (1 - y));
#     ym0.5 = y - 0.5;
#     log_odds_p = log(p / (1 - p));
#     pm0.5 = p - 0.5;
#     pcinv = 1 / (p * (1 - p));
#
#     A = rbind( # ND attempt
#         cbind( LP_M_coeffs(n, m, log_odds_y, ym0.5), -diag(m)),
#         cbind(-LP_M_coeffs(n, m, log_odds_y, ym0.5), -diag(m)),
#         cbind(LP_feas_coeffs(n, k, log_odds_p, pm0.5, pcinv), matrix(0, k, m))
#     );
#
#     # Right-hand side
#     diff_error = 0.001;
#     rhs = c(x, -x, rep(diff_error, k)); # ND attempt
#
#     # Build lpSolveAPI model
#     lpobj = lpSolveAPI::make.lp(nrow = nrow(A), ncol = ncol(A), verbose = "normal"); # ND attempt
#
#     lpSolveAPI::set.objfn(lpobj, cost);
#     for (i in 1:nrow(A)) {
#         lpSolveAPI::set.row(lpobj, row = i, x = A[i,]);
#     }
#     dir = c(rep("<=", 2 * m), rep(">=", k));
#     lpSolveAPI::set.constr.type(lpobj, dir);
#     lpSolveAPI::set.rhs(lpobj, rhs);
#     lpSolveAPI::set.bounds(lpobj, rep(-1e6, n), rep(1e6, n), columns = 1:n);
#
#     system.time(lpSolveAPI::solve.lpExtPtr(lpobj));
#     sol = lpSolveAPI::get.variables(lpobj);
#     coeffs = sol[1:n]
#
#     return(coeffs)
# }
