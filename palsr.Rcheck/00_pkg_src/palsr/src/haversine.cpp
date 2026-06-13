#include <Rcpp.h>
using namespace Rcpp;

// Great-circle (Haversine) distance between pairs of longitude/latitude points.
//
// Inputs are in decimal degrees and are recycled to a common length (like base R
// arithmetic), so any of the four coordinate arguments may be length 1. The result
// is returned in the same units as `radius` (default: kilometres, mean Earth radius
// 6371.0088 km). NA in any coordinate yields NA for that element.
//
// [[Rcpp::export]]
NumericVector haversine_cpp(NumericVector lon1, NumericVector lat1,
                            NumericVector lon2, NumericVector lat2,
                            double radius = 6371.0088) {
  R_xlen_t n1 = lon1.size(), n2 = lat1.size(), n3 = lon2.size(), n4 = lat2.size();
  R_xlen_t n = std::max(std::max(n1, n2), std::max(n3, n4));
  if (n1 == 0 || n2 == 0 || n3 == 0 || n4 == 0) return NumericVector(0);

  const double d2r = M_PI / 180.0;
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    double o1 = lon1[i % n1], a1 = lat1[i % n2];
    double o2 = lon2[i % n3], a2 = lat2[i % n4];
    if (NumericVector::is_na(o1) || NumericVector::is_na(a1) ||
        NumericVector::is_na(o2) || NumericVector::is_na(a2)) {
      out[i] = NA_REAL;
      continue;
    }
    double phi1 = a1 * d2r, phi2 = a2 * d2r;
    double dphi = (a2 - a1) * d2r;
    double dlam = (o2 - o1) * d2r;
    double s1 = std::sin(dphi / 2.0);
    double s2 = std::sin(dlam / 2.0);
    double h = s1 * s1 + std::cos(phi1) * std::cos(phi2) * s2 * s2;
    if (h > 1.0) h = 1.0;            // numerical guard for antipodal points
    out[i] = 2.0 * radius * std::asin(std::sqrt(h));
  }
  return out;
}
