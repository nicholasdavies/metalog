library(lpSolveAPI)

median_regression = function(x, y)
{
    # Ensure X and Y have same number of observations
    stopifnot(length(x) == length(y));

    n = length(x); # number of observations
    cost = c(0, 0, rep(1, n));
    A = rbind(
        cbind( x,  1, -diag(n)),
        cbind(-x, -1, -diag(n))
    );
    rhs = c(y, -y);

    # Build lpSolve model
    lpobj = make.lp(nrow = nrow(A), ncol = ncol(A));
    set.objfn(lpobj, cost);
    for (i in 1:nrow(A)) {
        set.row(lpobj, row = i, x = A[i,]);
    }
    set.constr.type(lpobj, rep("<=", nrow(A)));
    set.rhs(lpobj, rhs);
    set.bounds(lpobj, c(-Inf, -Inf), c(Inf, Inf), columns = c(1, 2));

    solve(lpobj)
    sol = get.variables(lpobj)

    return(sol[1:2])
}

npts = 100;
x = 1:npts + rnorm(npts, 0)
Y = 69 - 0.5 * x + rnorm(npts)

sol = median_regression(x, Y)
sol

plot(x, Y)
abline(sol[2], sol[1])


M_coeffs = function(n, m, log_odds_y, ym0.5)
{
    mat = matrix(0, nrow = m, ncol = n);

    # Eq. 6, Keelin 2016:
    mat[,1] = 1;
    mat[,2] = log_odds_y;
    if (n >= 3) {
        mat[,3] = ym0.5 * log_odds_y;
    }
    if (n >= 4) {
        mat[,4] = ym0.5;
    }
    if (n >= 5) {
        k = 2;
        for (j in seq(5, n, by = 2)) {
            factor = ym0.5 ^ k;
            mat[,j] = factor;
            if (j + 1 <= n) {
                mat[,j + 1] = factor * log_odds_y;
            }
            k = k + 1;
        }
    }

    mat
}

feas_coeffs = function(n, m, log_odds_y, ym0.5, ycinv)
{
    mat = matrix(0, nrow = m, ncol = n);

    # Eq. 6, Keelin 2016:
    mat[,2] = ycinv;
    if (n >= 3) {
        mat[,3] = ym0.5 * ycinv + log_odds_y;
    } else if (n >= 4) {
        mat[,4] = 1;
    }
    if (n >= 5) {
        k = 1;
        for (j in seq(5, n, by = 2)) {
            factor = ym0.5 ^ k;
            mat[,j] = (j-1)/2 * factor;
            if (j + 1 <= n) {
                mat[,j+1] = (factor * ym0.5 * ycinv + (k+1) * factor * log_odds_y);
            }
            k = k + 1;
        }
    }

    mat
}

fit_metalog = function(x, y, n)
{
    m = length(x); # number of CDF points

    # Test points for feasibility
    k = 100;
    p = seq(1e-4, 1 - 1e-4, length.out = k);

    # Objective function
    cost = c(rep(0, n), rep(1, m));

    # Constraints
    log_odds_y = log(y / (1 - y));
    ym0.5 = y - 0.5;
    log_odds_p = log(p / (1 - p));
    pm0.5 = p - 0.5;
    pcinv = 1 / (p * (1 - p));

    A = rbind(
        cbind( M_coeffs(n, m, log_odds_y, ym0.5), -diag(m)),
        cbind(-M_coeffs(n, m, log_odds_y, ym0.5), -diag(m)),
        cbind(feas_coeffs(n, k, log_odds_p, pm0.5, pcinv), matrix(0, k, n))
    );

    # Right-hand side
    rhs = c(x, -x, rep(0, k));

    # Build lpSolve model
    lpobj = make.lp(nrow = nrow(A), ncol = ncol(A), verbose = "normal");
    set.objfn(lpobj, cost);
    for (i in 1:nrow(A)) {
        set.row(lpobj, row = i, x = A[i,]);
    }
    set.constr.type(lpobj, rep("<=", 2 * m), 1:(2*m));
    set.constr.type(lpobj, rep(">=", k), (2*m+1):(2*m+k));
    set.rhs(lpobj, rhs);
    set.bounds(lpobj, rep(-Inf, n), rep(Inf, n), columns = 1:n);

    solve(lpobj)
    sol = get.variables(lpobj)

    return(sol[1:n])
}

sol = fit_metalog(c(1, 2, 5), c(0.1, 0.5, 0.9), 3)
plot(mlog(sol, "coeff"))
