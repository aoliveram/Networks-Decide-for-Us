// The diffusion engine of 05, in Rcpp. Loaded by scripts/05_unified_engine.R
// (Rcpp::sourceCpp); see that file for the gate and exposure helpers and the R
// reference implementations. Two functions:
//
//   run_exposure_grid()  the engine of the sweep: frozen gate, general exposure
//                        (one formula with theta and kappa; theta = kappa = 0 is
//                        the base model), the whole (MSP, IUL) grid in one call
//   run_unified_cpp()    one simulation under the annealed or the frozen gate with
//                        the base exposure; kept as the exact reference for 05's
//                        former R loop (scripts/05_engine_check.R)
//
// ---- run_unified_cpp() ------------------------------------------------------
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
#include <vector>
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

// ---- run_exposure_grid() ------------------------------------------------------
// Rule:  adopt_i  <=>  Gamma + lambda * E_i >= q_i  and  >= 1 adopting contact whose tie
// to i is open. Frozen gate: one coin per tie, drawn once per simulation right after
// set.seed(sim_seed), shared by both directions (the first E uniforms, as frozen_weights()).
//
// Exposure of ego i ("C1"; contacts split in circles g, "pairs" = ties between two
// contacts of i, i.e. triangles through i):
//
//   E_i = [ sum_g A(a_ig) + 2 kappa_p e_a_i ] / [ sum_g A(n_ig) + 2 kappa_p e_N_i ]
//   A(x) = x (1 + kappa_c (x - 1)) / (1 + theta (x - 1)),   A(0) = 0
//
// a_ig : contacts of i in circle g that adopted AND whose tie to i is open
// n_ig : contacts of i in circle g (all of them)
// e_a_i: ties between two adopting, open contacts of i;  e_N_i: ties between any two contacts
// theta in [0, 1]: independent reinforcement (Kish discount of redundant sources in a circle)
// kappa_p >= 0   : cohesive reinforcement (a bonus per pair of adopting contacts who know
//                  each other). kappa_c (cohesion by circle) is kept for checks only.
// theta = kappa = 0 with one circle per contact gives E_i = a_i / d_i, the base model; the
// function then reproduces run_unified_cpp() with the frozen gate cell by cell.
//
// Adoption is permanent and the gate frozen, so exposure only grows when a contact adopts:
// the cascade is monotone and stops at the first step without new adopters (a fixed point;
// max_steps is only a guard). Output per grid cell (rows h-major: row = hi * nG + gi):
// num_adopters, num_adopted_rational (non-seed adopters with q <= Gamma), t50 (step at which
// non-seed adoption reaches half its final value; 0 if none), num_steps (last step with an
// adoption).
static inline double Aform(double x, double th, double ka) {
  if (x <= 0.0) return 0.0;
  return x * (1.0 + ka * (x - 1.0)) / (1.0 + th * (x - 1.0));
}

// [[Rcpp::export]]
IntegerMatrix run_exposure_grid(IntegerVector src, IntegerVector dst, NumericMatrix W,
                                IntegerVector slot, IntegerVector slot_n, IntegerVector slot_ego,
                                double theta, double kappa_c,
                                IntegerVector tri_ptr, IntegerVector tri_f, double kappa_p,
                                IntegerVector eN,
                                NumericVector q, NumericVector Gammas, double lambda,
                                IntegerVector seeds, IntegerMatrix sim_seed, int max_steps) {
  const int n = q.size();
  const int m = src.size();          // directed edges: rows 0..E-1 and E..2E-1 (same tie twice)
  const int E = m / 2;
  const int nH = W.ncol(), nG = Gammas.size();
  const bool kp_on = kappa_p > 0.0;
  const bool plain = (theta == 0.0 && kappa_c == 0.0 && kappa_p == 0.0);
  const double eps = plain ? 0.0 : 1e-12;

  // out-edges by source (CSR)
  std::vector<int> out_ptr(n + 1, 0), out_e(m);
  for (int e = 0; e < m; ++e) ++out_ptr[src[e]];
  for (int i = 0; i < n; ++i) out_ptr[i + 1] += out_ptr[i];
  {
    std::vector<int> fo(out_ptr.begin(), out_ptr.end() - 1);
    for (int e = 0; e < m; ++e) out_e[fo[src[e] - 1]++] = e;
  }
  // denominators
  std::vector<double> denom(n, 0.0);
  const int S = slot_n.size();
  std::vector<double> Aslot_n(S);
  for (int s = 0; s < S; ++s) { Aslot_n[s] = Aform(slot_n[s], theta, kappa_c); denom[slot_ego[s] - 1] += Aslot_n[s]; }
  for (int i = 0; i < n; ++i) denom[i] += 2.0 * kappa_p * eN[i];

  std::vector<char> is_seed(n, 0);
  for (int s = 0; s < seeds.size(); ++s) is_seed[seeds[s] - 1] = 1;
  const int n_seed = (int) std::count(is_seed.begin(), is_seed.end(), 1);

  IntegerMatrix out(nH * nG, 4);
  colnames(out) = CharacterVector::create("num_adopters", "num_adopted_rational", "t50", "num_steps");

  std::vector<char> a(n), seen(n), act(m), open(m);
  std::vector<int> cnt(S), touched, newbies, ev_t, ev_n;
  std::vector<double> numer(n), tau(n);
  touched.reserve(n); newbies.reserve(n);
  Function set_seed("set.seed");

  for (int hi = 0; hi < nH; ++hi) {
    for (int gi = 0; gi < nG; ++gi) {
      const double G = Gammas[gi];
      for (int i = 0; i < n; ++i) tau[i] = std::max(0.0, (q[i] - G) / lambda);
      // frozen gate: set.seed + one coin per tie (first E uniforms), same as frozen_weights()
      set_seed(sim_seed(hi, gi));
      {
        RNGScope scope;
        for (int e = 0; e < E; ++e) {
          const char o = (unif_rand() <= W(e, hi)) ? 1 : 0;
          open[e] = o; open[e + E] = o;
        }
      }
      std::fill(a.begin(), a.end(), 0); std::fill(seen.begin(), seen.end(), 0);
      std::fill(act.begin(), act.end(), 0); std::fill(cnt.begin(), cnt.end(), 0);
      std::fill(numer.begin(), numer.end(), 0.0);
      touched.clear(); newbies.clear(); ev_t.clear(); ev_n.clear();
      int n_adopt = 0;
      for (int s = 0; s < seeds.size(); ++s) { a[seeds[s] - 1] = 1; }
      n_adopt = n_seed;

      auto activate = [&](int j) {
        for (int p = out_ptr[j]; p < out_ptr[j + 1]; ++p) {
          const int e = out_e[p];
          if (!open[e]) continue;
          const int i = dst[e] - 1;
          if (a[i]) continue;
          const int s = slot[e];
          const int c = ++cnt[s];
          numer[i] += Aform(c, theta, kappa_c) - Aform(c - 1, theta, kappa_c);
          if (kp_on) {
            int pairs = 0;
            for (int r = tri_ptr[e]; r < tri_ptr[e + 1]; ++r) pairs += act[tri_f[r]];
            numer[i] += 2.0 * kappa_p * pairs;
            act[e] = 1;
          }
          if (!seen[i]) { seen[i] = 1; touched.push_back(i); }
        }
      };

      for (int s = 0; s < seeds.size(); ++s) activate(seeds[s] - 1);
      int steps = 0;
      for (int t = 1; t <= max_steps; ++t) {
        newbies.clear();
        for (size_t k = 0; k < touched.size(); ++k) {
          const int i = touched[k];
          seen[i] = 0;
          if (!a[i] && numer[i] / denom[i] >= tau[i] - eps) newbies.push_back(i);
        }
        touched.clear();
        if (newbies.empty()) break;
        for (size_t k = 0; k < newbies.size(); ++k) a[newbies[k]] = 1;
        n_adopt += (int) newbies.size();
        for (size_t k = 0; k < newbies.size(); ++k) activate(newbies[k]);
        steps = t; ev_t.push_back(t); ev_n.push_back(n_adopt);
        if (n_adopt == n) break;
      }
      int rational = 0;
      for (int i = 0; i < n; ++i) if (a[i] && !is_seed[i] && q[i] <= G) ++rational;
      int t50 = 0;
      const double half = 0.5 * (double) (n_adopt - n_seed);
      for (size_t k = 0; k < ev_t.size(); ++k) if ((double) (ev_n[k] - n_seed) >= half) { t50 = ev_t[k]; break; }
      const int row = hi * nG + gi;
      out(row, 0) = n_adopt; out(row, 1) = rational; out(row, 2) = t50; out(row, 3) = steps;
    }
  }
  return out;
}
