// Rcpp version of run_unified() in scripts/05_unified_diffusion_sweep.R
// Loaded by scripts/05_unified_engine.R (Rcpp::sourceCpp); see that file for
// the gate helpers and the R reference implementations.
//
// Same algorithm and the SAME random stream as the R loop: set.seed() is called
// from R before each simulation, and each frontier tie draws one unif_rand() in
// increasing edge order, exactly as runif(length(idx)) does in R. Results are
// therefore identical to 05 (checked by scripts/05_engine_check.R).
//
// The frontier (ties from an adopter to a non-adopter) is updated incrementally
// instead of being recomputed over all ties at every step.
//
// Extra outputs (not in 05): why the run stopped, the step at which it stopped,
// t50 (step at which non-seed adoption reaches half its final value) and the
// number of adopters at fixed horizons (the final value if it stopped earlier).
#include <Rcpp.h>
#include <algorithm>
#include <string>
using namespace Rcpp;

// src, dst : every tie stored twice (1-based node ids), as in 05's edge_arrays()
// deg      : degree, floored at 1 (the denominator of the exposure)
// w_edge   : pass probability of each directed tie (0/1 under the frozen gate)
// tau_star : max(0, (q - Gamma) / lambda), the share of contacts each node needs
// [[Rcpp::export]]
IntegerVector run_unified_cpp(IntegerVector src, IntegerVector dst, IntegerVector deg,
                              NumericVector w_edge, NumericVector tau_star,
                              NumericVector q, double Gamma,
                              IntegerVector seeds, int max_steps, int patience,
                              IntegerVector horizons) {
  const int n = tau_star.size();
  const int m = src.size();
  std::vector<char> a(n, 0), is_seed(n, 0);
  for (int s = 0; s < seeds.size(); ++s) { a[seeds[s] - 1] = 1; is_seed[seeds[s] - 1] = 1; }
  int n_adopt = 0;
  for (int i = 0; i < n; ++i) n_adopt += a[i];
  const int n_seed = n_adopt;

  std::vector<int> counts(n, 0), touched, newbies;
  touched.reserve(n); newbies.reserve(n);
  // stop reason: 1 all adopted, 2 frontier empty, 3 patience (a fixed point under
  // the frozen gate, where patience = 1), 4 step limit
  int stall = 0, steps = 0, reason = 4;
  int t_stop = max_steps;
  const int H = horizons.size();
  IntegerVector at_h(H, -1);
  int hi = 0;
  std::vector<int> ev_t, ev_n;   // adoption events: step and cumulative adopters

  // ties by source and by destination (CSR), to update the frontier
  std::vector<int> out_ptr(n + 1, 0), in_ptr(n + 1, 0), out_e(m), in_e(m);
  for (int e = 0; e < m; ++e) { ++out_ptr[src[e]]; ++in_ptr[dst[e]]; }
  for (int i = 0; i < n; ++i) { out_ptr[i + 1] += out_ptr[i]; in_ptr[i + 1] += in_ptr[i]; }
  {
    std::vector<int> fo(out_ptr.begin(), out_ptr.end() - 1), fi(in_ptr.begin(), in_ptr.end() - 1);
    for (int e = 0; e < m; ++e) { out_e[fo[src[e] - 1]++] = e; in_e[fi[dst[e] - 1]++] = e; }
  }
  // frontier kept in increasing edge order: the order in which R's which() +
  // runif() draws the ties
  std::vector<char> in_f(m, 0);
  std::vector<int> frontier;
  frontier.reserve(m);
  for (int e = 0; e < m; ++e)
    if (a[src[e] - 1] && !a[dst[e] - 1]) { in_f[e] = 1; frontier.push_back(e); }

  RNGScope scope;
  for (int t = 1; t <= max_steps; ++t) {
    if (frontier.empty()) { reason = 2; t_stop = t; break; }
    touched.clear();
    for (size_t k = 0; k < frontier.size(); ++k) {
      const int e = frontier[k];
      if (unif_rand() <= w_edge[e]) {
        const int d = dst[e] - 1;
        if (counts[d] == 0) touched.push_back(d);
        ++counts[d];
      }
    }
    newbies.clear();
    for (size_t k = 0; k < touched.size(); ++k) {
      const int i = touched[k];
      const double E = (double) counts[i] / (double) deg[i];
      if (!a[i] && E >= tau_star[i]) newbies.push_back(i);
      counts[i] = 0;
    }
    for (size_t k = 0; k < newbies.size(); ++k) a[newbies[k]] = 1;
    if (!newbies.empty()) {
      // ties into the new adopters leave the frontier, their ties to
      // non-adopters enter; then restore increasing edge order
      bool added = false;
      for (size_t k = 0; k < newbies.size(); ++k) {
        const int i = newbies[k];
        for (int p = in_ptr[i]; p < in_ptr[i + 1]; ++p) in_f[in_e[p]] = 0;
        for (int p = out_ptr[i]; p < out_ptr[i + 1]; ++p) {
          const int e = out_e[p];
          if (!a[dst[e] - 1] && !in_f[e]) { in_f[e] = 1; frontier.push_back(e); added = true; }
        }
      }
      size_t w = 0;
      for (size_t k = 0; k < frontier.size(); ++k) if (in_f[frontier[k]]) frontier[w++] = frontier[k];
      frontier.resize(w);
      if (added) std::sort(frontier.begin(), frontier.end());
      n_adopt += newbies.size(); stall = 0; steps = t;
      ev_t.push_back(t); ev_n.push_back(n_adopt);
    } else ++stall;
    while (hi < H && horizons[hi] == t) { at_h[hi] = n_adopt; ++hi; }
    if (n_adopt == n) { reason = 1; t_stop = t; break; }
    if (stall >= patience) { reason = 3; t_stop = t; break; }
  }
  for (; hi < H; ++hi) at_h[hi] = n_adopt;

  int rational = 0, social = 0;
  for (int i = 0; i < n; ++i) if (a[i] && !is_seed[i]) { if (q[i] <= Gamma) ++rational; else ++social; }

  int t50 = 0;   // 0 if no non-seed node adopted
  const double half = 0.5 * (double) (n_adopt - n_seed);
  for (size_t k = 0; k < ev_t.size(); ++k) if ((double) (ev_n[k] - n_seed) >= half) { t50 = ev_t[k]; break; }

  IntegerVector out = IntegerVector::create(
    _["num_adopters"] = n_adopt, _["num_adopted_rational"] = rational,
    _["num_adopted_social"] = social, _["num_steps"] = steps,
    _["initial_cluster_size"] = n_seed, _["stop_reason"] = reason, _["t_stop"] = t_stop,
    _["t50"] = t50);
  CharacterVector nm = out.names();
  for (int k = 0; k < H; ++k) {
    out.push_back(at_h[k]);
    nm.push_back("n_t" + std::to_string(horizons[k]));
  }
  out.names() = nm;
  return out;
}
