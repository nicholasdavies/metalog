// xtabulate.cpp

#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Tabulation for numeric vectors
//'
//' Base-R's \code{\link[base]{tabulate}} will not tabulate numeric values; this
//' function does.
//'
//' @param y a numeric vector.
//'
//' @return A list with two components: \code{x}, the sorted unique values in
//' \code{y}, and \code{n}, the count associated with each value in \code{x}.
// [[Rcpp::export]]
List xtabulate(NumericVector& y)
{
    // To create storage for x, we need to explicitly copy y, as Rcpp passes by
    // reference even when the argument to the C++ function is passed by value.
    NumericVector x = clone(y); // Sorted values
    IntegerVector n(x.size());  // Counts of each value in x

    if (x.size() > 0)
    {
        // Sort values in x
        std::sort(x.begin(), x.end());

        // Count each value in x
        auto counter = n.begin();
        auto item = x.begin(), last = x.end();
        auto dest = item;

        // Step through all elements of x
        while (++item != last)
        {
            // Increment counter for the current element
            ++*counter;

            // If item has reached a new value, save it and point
            // dest and counter to the next position
	        if (*dest != *item)
	        {
	            *++dest = *item;
	            ++counter;
	        }
        }
        // Increment counter for the final element
        ++*counter;

        // Resize x and n
        x.erase(dest + 1, last);
        n.erase(counter + 1, n.end());
    }

    // Return results in a list
    List L = List::create(_["x"] = x, _["n"] = n);
    return L;
}
