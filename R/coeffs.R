LP_feas_coeffs = function(n, m, log_odds_y, ym0.5, ycinv)
{
    mat = matrix(0, nrow = m, ncol = n);

    # Eq. 6, Keelin 2016:
    mat[,2] = ycinv;
    if (n >= 3) {
        mat[,3] = ym0.5 * ycinv + log_odds_y;
    }
    if (n >= 4) {
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

    return (mat)
}

mean_coeffs = function(n)
{
    mu = rep(0, n);

    for (k in 1:n) {
        if (k <= 4) {
            mu[k] = c(1, 0, 1/2, 0)[k];
        } else if (k %% 4 == 0) {
            mu[k] = 4/k * 2^(-1-k/2) *
                (2*sum(1/(1:(k/2))) - sum(1/(1:(k/4) - 2*(1:(k/4))^2)));
        } else if (k %% 4 == 1) {
            mu[k] = 2^((3-k)/2) / (1 + k);
        } else {
            mu[k] = 0;
        }
    }

    return (mu);
}

Y_matrix = function(m, y, mean, n)
{
    # Create matrix Y which is used to solve for metalog coefficients.
    # Eq. 7, Keelin 2016:
    log_odds_y = log(y / (1 - y));
    ym0.5 = y - 0.5;
    Y = matrix(1, nrow = m, ncol = n);

    # Top rows: quantiles
    top = seq_along(y);
    Y[top, 2] = log_odds_y;
    if (n >= 3) {
        Y[top, 3] = ym0.5 * log_odds_y;
    }
    if (n >= 4) {
        Y[top, 4] = ym0.5;
    }
    if (n >= 5) {
        k = 2;
        for (j in seq(5, n, by = 2)) {
            factor = ym0.5 ^ k;
            Y[top, j] = factor;
            if (j + 1 <= n) {
                Y[top, j + 1] = factor * log_odds_y;
            }
            k = k + 1;
        }
    }

    # Bottom row: mean (if requested)
    if (!is.null(mean)) {
        Y[m, ] = mean_coeffs(n);
    }

    return (Y)
}

Yc_matrix = function(m, yc, mean, n)
{
    # Create matrix Y which is used to solve for metalog coefficients.
    # Eq. 7, Keelin 2016:
    log_odds_y = log((1 - yc) / yc);
    ym0.5 = 0.5 - yc;
    Y = matrix(1, nrow = m, ncol = n);

    # Top rows: quantiles
    top = 1:length(yc);
    Y[top, 2] = log_odds_y;
    if (n >= 3) {
        Y[top, 3] = ym0.5 * log_odds_y;
    }
    if (n >= 4) {
        Y[top, 4] = ym0.5;
    }
    if (n >= 5) {
        k = 2;
        for (j in seq(5, n, by = 2)) {
            factor = ym0.5 ^ k;
            Y[top, j] = factor;
            if (j + 1 <= n) {
                Y[top, j + 1] = factor * log_odds_y;
            }
            k = k + 1;
        }
    }

    # Bottom row: mean (if requested)
    if (!is.null(mean)) {
        Y[m, ] = mean_coeffs(n);
    }

    return (Y)
}
