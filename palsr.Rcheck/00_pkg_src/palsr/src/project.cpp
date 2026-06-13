#include <Rcpp.h>
using namespace Rcpp;

// Project a single actor's location at one prediction time from pre-extracted event
// histories. All history filtering (which events belong to the focal actor / its
// alters, and the age = predict_time - event_time in days) is done in R and passed in;
// this kernel only performs the exponential-smoothing weighting and the focal/alter
// blend. See ALGORITHM.md.
//
// Returns a length-4 numeric vector: {lon, lat, n_focal, n_alter}. If the focal actor
// has no prior events, {lon, lat} are NA (no projected location).
//
// pi_zero      : force the mixing weight pi = 0 (one-parameter / focal-only model).
// alter_legacy : reproduce the Dataverse code's scalar alter weight (un-normalized sum
//                of alter coordinates) instead of the paper's normalized weighted mean.
// eps          : numerical offset inside each age weight (default 0.01).
//
// [[Rcpp::export]]
NumericVector project_one_cpp(NumericVector focal_age,
                              NumericVector focal_lon,
                              NumericVector focal_lat,
                              NumericVector alter_age,
                              NumericVector alter_lon,
                              NumericVector alter_lat,
                              double alpha, double beta,
                              double gamma, double eta,
                              bool pi_zero = false,
                              bool alter_legacy = false,
                              double eps = 0.01) {
  R_xlen_t n_i = focal_age.size();
  R_xlen_t n_k = alter_age.size();
  NumericVector out(4);

  if (n_i == 0) {                      // no focal history -> no PAL
    out[0] = NA_REAL; out[1] = NA_REAL; out[2] = 0; out[3] = (double) n_k;
    return out;
  }

  // Focal: normalized recency weights, weighted mean of coordinates.
  double wsum = 0.0;
  std::vector<double> wi(n_i);
  for (R_xlen_t e = 0; e < n_i; ++e) {
    wi[e] = 1.0 / (std::pow(focal_age[e], alpha) + eps);
    wsum += wi[e];
  }
  double foc_lon = 0.0, foc_lat = 0.0;
  for (R_xlen_t e = 0; e < n_i; ++e) {
    double w = wi[e] / wsum;
    foc_lon += w * focal_lon[e];
    foc_lat += w * focal_lat[e];
  }

  // Mixing weight pi (weight on the alter average).
  double pi = 0.0;
  if (!pi_zero && n_k > 0) {
    double v = std::pow((double) n_i / (double) n_k, 1.0 / std::sqrt((double) n_k));
    pi = 1.0 / (1.0 + std::exp(-(gamma + eta * v)));   // plogis
  }

  double alt_lon = 0.0, alt_lat = 0.0;
  if (pi > 0.0 && n_k > 0) {
    if (alter_legacy) {
      // Dataverse behavior: scalar weight = 1 / sum(unnormalized weights), applied to
      // the (un-normalized) sum of alter coordinates.
      double ksum = 0.0;
      for (R_xlen_t e = 0; e < n_k; ++e)
        ksum += 1.0 / (std::pow(alter_age[e], beta) + eps);
      double scal = 1.0 / ksum;
      double slon = 0.0, slat = 0.0;
      for (R_xlen_t e = 0; e < n_k; ++e) { slon += alter_lon[e]; slat += alter_lat[e]; }
      alt_lon = scal * slon;
      alt_lat = scal * slat;
    } else {
      // Paper's intended normalized weighted mean.
      double ksum = 0.0;
      std::vector<double> wk(n_k);
      for (R_xlen_t e = 0; e < n_k; ++e) {
        wk[e] = 1.0 / (std::pow(alter_age[e], beta) + eps);
        ksum += wk[e];
      }
      for (R_xlen_t e = 0; e < n_k; ++e) {
        double w = wk[e] / ksum;
        alt_lon += w * alter_lon[e];
        alt_lat += w * alter_lat[e];
      }
    }
  }

  out[0] = (1.0 - pi) * foc_lon + pi * alt_lon;
  out[1] = (1.0 - pi) * foc_lat + pi * alt_lat;
  out[2] = (double) n_i;
  out[3] = (double) n_k;
  return out;
}

// Batched projection over many "units" (actor-at-time histories) sharing one parameter
// set, used inside the estimation objective so each optimizer evaluation is a single
// call. Histories are flattened: f_* are the concatenated focal events with per-unit
// offset f_off (0-based) and length f_len; k_* likewise for alter events. Returns an
// n-by-2 matrix of projected {lon, lat}; rows with no focal history are NA.
//
// [[Rcpp::export]]
NumericMatrix project_batch_cpp(NumericVector f_age, NumericVector f_lon, NumericVector f_lat,
                                IntegerVector f_off, IntegerVector f_len,
                                NumericVector k_age, NumericVector k_lon, NumericVector k_lat,
                                IntegerVector k_off, IntegerVector k_len,
                                double alpha, double beta, double gamma, double eta,
                                bool pi_zero = false, bool alter_legacy = false,
                                double eps = 0.01) {
  R_xlen_t n = f_off.size();
  NumericMatrix out(n, 2);
  for (R_xlen_t u = 0; u < n; ++u) {
    int n_i = f_len[u], n_k = k_len[u];
    int fo = f_off[u], ko = k_off[u];
    if (n_i == 0) { out(u, 0) = NA_REAL; out(u, 1) = NA_REAL; continue; }

    double wsum = 0.0;
    for (int e = 0; e < n_i; ++e) wsum += 1.0 / (std::pow(f_age[fo + e], alpha) + eps);
    double foc_lon = 0.0, foc_lat = 0.0;
    for (int e = 0; e < n_i; ++e) {
      double w = (1.0 / (std::pow(f_age[fo + e], alpha) + eps)) / wsum;
      foc_lon += w * f_lon[fo + e];
      foc_lat += w * f_lat[fo + e];
    }

    double pi = 0.0;
    if (!pi_zero && n_k > 0) {
      double v = std::pow((double) n_i / (double) n_k, 1.0 / std::sqrt((double) n_k));
      pi = 1.0 / (1.0 + std::exp(-(gamma + eta * v)));
    }
    double alt_lon = 0.0, alt_lat = 0.0;
    if (pi > 0.0 && n_k > 0) {
      if (alter_legacy) {
        double ksum = 0.0, slon = 0.0, slat = 0.0;
        for (int e = 0; e < n_k; ++e) {
          ksum += 1.0 / (std::pow(k_age[ko + e], beta) + eps);
          slon += k_lon[ko + e]; slat += k_lat[ko + e];
        }
        double scal = 1.0 / ksum;
        alt_lon = scal * slon; alt_lat = scal * slat;
      } else {
        double ksum = 0.0;
        for (int e = 0; e < n_k; ++e) ksum += 1.0 / (std::pow(k_age[ko + e], beta) + eps);
        for (int e = 0; e < n_k; ++e) {
          double w = (1.0 / (std::pow(k_age[ko + e], beta) + eps)) / ksum;
          alt_lon += w * k_lon[ko + e];
          alt_lat += w * k_lat[ko + e];
        }
      }
    }
    out(u, 0) = (1.0 - pi) * foc_lon + pi * alt_lon;
    out(u, 1) = (1.0 - pi) * foc_lat + pi * alt_lat;
  }
  return out;
}
