#' Mean, variance, and other moments of a metalog distribution
#'
#' Calculate the moments of a metalog distribution \code{ml}. This function
#' allows calculating the mean, variance, and other quantities of potential
#' interest.
#'
#' The mean of a distribution is its first noncentral moment and the variance
#' \ifelse{html}{\eqn{\sigma}\out{<sup>2</sup>}}{\eqn{\sigma^2}} is the second
#' central moment. Standardly, the skewness is defined as the third standardized
#' moment and the kurtosis as the fourth standardized moment, where the
#' \emph{n}th standardized moment is the \emph{n}th central moment divided by
#' \ifelse{html}{\eqn{\sigma}\out{<sup><i>n</i></sup>}}{\eqn{\sigma^n}}.
#'
#' In other words:
#' \itemize{
#'  \item{Mean: \code{moment(ml)}}
#'  \item{Variance: \code{moment(ml, 2, TRUE)}}
#'  \item{Skewness: \code{moment(ml, 3, TRUE) / moment(ml, 2, TRUE)^(3/2)}}
#'  \item{Kurtosis: \code{moment(ml, 4, TRUE) / moment(ml, 2, TRUE)^2}}
#' }
#'
#' Mathematically, this function calculates the integral over \eqn{y} from 0 to
#' 1 of \ifelse{html}{\eqn{f(M(y)-c)}\out{<sup><i>k</i></sup>}}{\eqn{f(M(y)-c)^k}},
#' where: \eqn{M(y)} is the quantile function of the distribution; \eqn{k} is
#' \code{order}; \eqn{c} is equal to 0 if \code{central == FALSE} and is equal
#' to the mean of the distribution if \code{central == TRUE}; and \eqn{f(x)=x}
#' if \code{absolute == FALSE} and \eqn{f(x)=|x|} if \code{absolute == TRUE}.
#'
#' @param order Which moment to calculate.
#' @param central \code{TRUE} for moment about the mean (central moment),
#'  \code{FALSE} for moment about 0.
#' @param absolute \code{TRUE} to calculate absolute moments.
#'
#' @export
moment = function(ml, order = 1, central = FALSE, absolute = FALSE)
{
    # TODO make this work as an S3 method?? see packages e1071 and moments for
    # another definition of this function . . .
    # note e1071 has 'center'

    # Optimised cases for 1st central moment (=0) and 1st noncentral moment (=mean)
    if (!absolute && order == 1) {
        if (central) {
            return (0)
        } else if (identical(ml$b, c(-Inf, Inf))) {
            return (sum(mean_coeffs(length(ml$a)) * ml$a))
        }
    }

    # Origin for moment calculation
    if (central) {
        origin = moment(ml);
    } else {
        origin = 0;
    }

    # Calculate moment
    if (absolute) {
        integrate01(
            f = function(y) abs(metalog_M(y, a = ml$a, n = length(ml$a),
                    bl = ml$b[1], bu = ml$b[2]) - origin)^order,
            fc = function(yc) abs(metalog_Mc(yc, a = ml$a, n = length(ml$a),
                    bl = ml$b[1], bu = ml$b[2]) - origin)^order)
    } else {
        integrate01(
            f = function(y) (metalog_M(y, a = ml$a, n = length(ml$a),
                    bl = ml$b[1], bu = ml$b[2]) - origin)^order,
            fc = function(yc) (metalog_Mc(yc, a = ml$a, n = length(ml$a),
                    bl = ml$b[1], bu = ml$b[2]) - origin)^order)
    }
}

# Tanh-sinh quadrature over the interval [0,1]
integrate01 = function(f, N = 64, h = 1/16, fc = function(x) f(1 - x))
{
    # Implementation guided by DH Bailey, Tanh-Sinh High-Precision Quadrature,
    # https://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf.
    # References in comments here are to this paper.

    # Values of j for sum in eq. (1), but only from 0 to N
    j = 0:N;

    # Complement of x_j, y_j = 1 - x_j, referred to in section 3
    y = 1/(exp(0.5 * pi * sinh(j * h)) * cosh(0.5 * pi * sinh(j * h)));

    # Weights w_j
    w = (0.5 * pi * h * cosh(j * h)) / (cosh(0.5 * pi * sinh(j * h))^2);

    # Conduct integration in two halves using f and fc
    # Left half, including j = 0
    lh = sum(w * f(y * 0.5));

    # Right half, excluding j = 0
    rh = sum(w[-1] * fc(y[-1] * 0.5));

    return (0.5 * (rh + lh))
}

# Get specified standardized moments of the metalog distribution ml
# Note that for orders == 1, this returns the mean; for orders == 2, the
# standard deviation; and for all other orders, the standardized order (nth
# central order divided by the standard deviation to the nth power).
standard_moments = function(ml, orders)
{
    m = rep(NA_real_, length(orders));

    # Get standard deviation for standardizing, if needed
    if (any(orders > 1)) {
        sd = sqrt(moment(ml, 2, central = TRUE));
    }

    # Get all requested orders
    for (i in seq_along(orders))
    {
        k = orders[i];
        if (k == 1) {
            m[i] = moment(ml, 1, central = FALSE);
        } else if (k == 2) {
            m[i] = sd;
        } else {
            m[i] = moment(ml, k, central = TRUE) / sd^k;
        }
    }

    return (m);
}
