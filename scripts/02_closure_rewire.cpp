// Constrained rewiring toward more triadic closure (scripts/02_network_calibrated.R).
// Loaded with Rcpp::sourceCpp from that script. Two targets:
//   rewire_closure_node()  every node toward target_factor x its own triangles: level x2 (main)
//   rewire_closure()       the TOTAL triangles toward target_factor x: level x2g (validation only)
//
// Attributes (one row per node):
//   cats  integer matrix of categorical attributes (GSS/ATP/GSS14: race, sex, religion;
//         ELSOC: education group, sex, religion, ideology)
//   nums  numeric matrix of continuous attributes (age, education in years)
// Since PR #5 the attributes come as matrices (until then: race, sex, relig, age, educ as
// separate vectors). The checks consume no random numbers, so GSS and ATP networks are
// bit-identical to the earlier version.
//
// rewire_closure():
// Swap: edges (a,b),(c,d) -> (a,d),(c,b). Every node keeps its degree. The swap is allowed only if,
// for every categorical attribute separately, x(b) == x(d) or x(a) == x(c): then the multiset of
// attribute pairs over edges is unchanged, so the mixing matrix of each attribute (and its
// same-attribute count) stays exactly as it was. The total |difference| of every continuous attribute
// over edges must stay within `tol` (relative) of its starting value. Moves are accepted toward a
// target number of triangles (Metropolis on |T - T*|, temperature `temp`); neutral moves are accepted.
// After T* is reached, a mixing phase of `mix_props` more proposals keeps |T - T*| small while the
// network keeps changing. With degrees fixed, the number of connected triples is fixed, so
// transitivity scales with T.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// the swap (a,b),(c,d) -> (a,d),(c,b) keeps every categorical mixing matrix
static inline bool cats_ok(const IntegerMatrix& cats, int a, int b, int c, int d) {
  for (int k = 0; k < cats.ncol(); k++)
    if (!(cats(b, k) == cats(d, k) || cats(a, k) == cats(c, k))) return false;
  return true;
}
// new totals of |difference| of every continuous attribute; false if one leaves the tolerance
static inline bool nums_ok(const NumericMatrix& x, int a, int b, int c, int d, const std::vector<double>& s,
                           const std::vector<double>& s0, double tol, std::vector<double>& ns) {
  for (int m = 0; m < x.ncol(); m++) {
    ns[m] = s[m] - std::fabs(x(a, m) - x(b, m)) - std::fabs(x(c, m) - x(d, m)) + std::fabs(x(a, m) - x(d, m)) + std::fabs(x(c, m) - x(b, m));
    if (std::fabs(ns[m] - s0[m]) > tol * s0[m]) return false;
  }
  return true;
}

// [[Rcpp::export]]
List rewire_closure(IntegerMatrix el, int n, IntegerMatrix cats, NumericMatrix nums, double target_factor, double tol,
                    double max_props, double temp, double mix_props) {
  const int m = el.nrow(), M = nums.ncol();
  std::vector<int> A(m), B(m);
  std::vector<unsigned char> adj((size_t)n * n, 0);
  std::vector<std::vector<int> > nb(n);
  for (int e = 0; e < m; e++) {
    int a = el(e, 0) - 1, b = el(e, 1) - 1;
    A[e] = a; B[e] = b;
    adj[(size_t)a * n + b] = adj[(size_t)b * n + a] = 1;
    nb[a].push_back(b); nb[b].push_back(a);
  }
  auto cn = [&](int u, int v) {
    const bool su = nb[u].size() < nb[v].size();
    const std::vector<int>& s = su ? nb[u] : nb[v];
    const int o = su ? v : u;
    int c = 0;
    for (int x : s) if (adj[(size_t)x * n + o]) c++;
    return c;
  };
  auto rm_nb = [&](int u, int v) {
    std::vector<int>& L = nb[u];
    for (size_t i = 0; i < L.size(); i++) if (L[i] == v) { L[i] = L.back(); L.pop_back(); break; }
  };
  double T = 0;
  std::vector<double> s0(M, 0.0), ns(M, 0.0);
  for (int e = 0; e < m; e++) {
    T += cn(A[e], B[e]);
    for (int k = 0; k < M; k++) s0[k] += std::fabs(nums(A[e], k) - nums(B[e], k));
  }
  T /= 3.0;
  const double T0 = T, Tstar = target_factor * T0;
  std::vector<double> s = s0;
  double props = 0;
  long acc = 0;
  RNGScope scope;
  // phase 1: climb to T*; phase 2 (mix_props proposals after reaching T*): keep |T - T*| small while mixing
  double after = 0; bool reached = false;
  while (props < max_props && (!reached || after < mix_props)) {
    props++; if (reached) after++;
    int e1 = (int)(R::unif_rand() * m), e2 = (int)(R::unif_rand() * m);
    if (e1 == e2) continue;
    int a = A[e1], b = B[e1], c = A[e2], d = B[e2];
    if (R::unif_rand() < 0.5) { int t = c; c = d; d = t; }
    if (a == c || a == d || b == c || b == d) continue;
    if (adj[(size_t)a * n + d] || adj[(size_t)c * n + b]) continue;
    if (!cats_ok(cats, a, b, c, d)) continue;
    if (!nums_ok(nums, a, b, c, d, s, s0, tol, ns)) continue;
    int lost = cn(a, b) + cn(c, d);
    int gain = (cn(a, d) - adj[(size_t)b * n + d] - adj[(size_t)c * n + a]) +
               (cn(c, b) - adj[(size_t)d * n + b] - adj[(size_t)a * n + c]);
    double dT = gain - lost;
    double E0 = std::fabs(T - Tstar), E1 = std::fabs(T + dT - Tstar);
    if (E1 > E0 && R::unif_rand() >= std::exp(-(E1 - E0) / temp)) continue;
    adj[(size_t)a * n + b] = adj[(size_t)b * n + a] = 0; adj[(size_t)c * n + d] = adj[(size_t)d * n + c] = 0;
    adj[(size_t)a * n + d] = adj[(size_t)d * n + a] = 1; adj[(size_t)c * n + b] = adj[(size_t)b * n + c] = 1;
    rm_nb(a, b); rm_nb(b, a); rm_nb(c, d); rm_nb(d, c);
    nb[a].push_back(d); nb[d].push_back(a); nb[c].push_back(b); nb[b].push_back(c);
    A[e1] = a; B[e1] = d; A[e2] = c; B[e2] = b;
    T += dT; s = ns; acc++;
    if (T >= Tstar) reached = true;
  }
  IntegerMatrix out(m, 2);
  for (int e = 0; e < m; e++) { out(e, 0) = A[e] + 1; out(e, 1) = B[e] + 1; }
  NumericVector rel(M);
  for (int k = 0; k < M; k++) rel[k] = s[k] / s0[k] - 1;
  return List::create(_["edges"] = out, _["T0"] = T0, _["T"] = T, _["accepted"] = (double)acc,
                      _["proposals"] = props, _["num_rel_change"] = rel);
}

// Variant "node": the same swaps and constraints, but the target is per node: every node i
// moves toward target_factor x its own starting triangles t_i(0). Energy E = sum_i |t_i - t_i*|
// (Metropolis, temperature `temp`). The global target of rewire_closure() lets the triangles go
// where they are cheapest, i.e. among high-degree nodes; here the local clustering of every node
// is multiplied by the same factor, so the clustering-degree profile of the starting network is
// kept. Phase 1 ends when E <= e_tol x sum_i t_i*; then `mix_props` proposals as above.
// [[Rcpp::export]]
List rewire_closure_node(IntegerMatrix el, int n, IntegerMatrix cats, NumericMatrix nums, double target_factor, double tol,
                         double max_props, double temp, double mix_props, double e_tol) {
  const int m = el.nrow(), M = nums.ncol();
  std::vector<int> A(m), B(m);
  std::vector<unsigned char> adj((size_t)n * n, 0);
  std::vector<std::vector<int> > nb(n);
  for (int e = 0; e < m; e++) {
    int a = el(e, 0) - 1, b = el(e, 1) - 1;
    A[e] = a; B[e] = b;
    adj[(size_t)a * n + b] = adj[(size_t)b * n + a] = 1;
    nb[a].push_back(b); nb[b].push_back(a);
  }
  auto rm_nb = [&](int u, int v) {
    std::vector<int>& L = nb[u];
    for (size_t i = 0; i < L.size(); i++) if (L[i] == v) { L[i] = L.back(); L.pop_back(); break; }
  };
  // triangles at every node
  std::vector<int> t(n, 0);
  for (int e = 0; e < m; e++) {
    const int u = A[e], v = B[e];
    const std::vector<int>& s = nb[u].size() < nb[v].size() ? nb[u] : nb[v];
    const int o = nb[u].size() < nb[v].size() ? v : u;
    for (int x : s) if (adj[(size_t)x * n + o]) { t[u]++; t[v]++; t[x]++; }
  }
  for (int i = 0; i < n; i++) t[i] /= 3;            // each triangle was seen from its 3 edges
  std::vector<double> ts(n);
  double Ets = 0, E = 0, T = 0;
  for (int i = 0; i < n; i++) { ts[i] = target_factor * t[i]; Ets += ts[i]; E += std::fabs(t[i] - ts[i]); T += t[i]; }
  T /= 3.0;
  const double T0 = T, E0 = E;
  std::vector<double> s0(M, 0.0), ns(M, 0.0);
  for (int e = 0; e < m; e++) for (int k = 0; k < M; k++) s0[k] += std::fabs(nums(A[e], k) - nums(B[e], k));
  std::vector<double> sv = s0;
  double props = 0;
  long acc = 0;
  std::vector<int> dlt(n, 0), touched;
  std::vector<char> in_t(n, 0);
  touched.reserve(4096);
  auto bump = [&](int u, int s) { if (!in_t[u]) { in_t[u] = 1; touched.push_back(u); } dlt[u] += s; };
  // common neighbours of (u, v), skipping x1 and x2 (the endpoints of the removed ties)
  auto each_cn = [&](int u, int v, int x1, int x2, int s) {
    const bool su = nb[u].size() < nb[v].size();
    const std::vector<int>& L = su ? nb[u] : nb[v];
    const int o = su ? v : u;
    int c = 0;
    for (int x : L) if (x != x1 && x != x2 && adj[(size_t)x * n + o]) { bump(u, s); bump(v, s); bump(x, s); c++; }
    return c;
  };
  RNGScope scope;
  double after = 0; bool reached = false;
  while (props < max_props && (!reached || after < mix_props)) {
    props++; if (reached) after++;
    int e1 = (int)(R::unif_rand() * m), e2 = (int)(R::unif_rand() * m);
    if (e1 == e2) continue;
    int a = A[e1], b = B[e1], c = A[e2], d = B[e2];
    if (R::unif_rand() < 0.5) { int tt = c; c = d; d = tt; }
    if (a == c || a == d || b == c || b == d) continue;
    if (adj[(size_t)a * n + d] || adj[(size_t)c * n + b]) continue;
    if (!cats_ok(cats, a, b, c, d)) continue;
    if (!nums_ok(nums, a, b, c, d, sv, s0, tol, ns)) continue;
    // triangle changes per node: remove (a,b) and (c,d), then add (a,d) and (c,b)
    int dT = 0;
    dT -= each_cn(a, b, -1, -1, -1);
    dT -= each_cn(c, d, -1, -1, -1);
    dT += each_cn(a, d, b, c, +1);
    dT += each_cn(c, b, d, a, +1);
    double dE = 0;
    for (int u : touched) dE += std::fabs(t[u] + dlt[u] - ts[u]) - std::fabs(t[u] - ts[u]);
    const bool ok = !(dE > 0 && R::unif_rand() >= std::exp(-dE / temp));
    for (int u : touched) { if (ok) t[u] += dlt[u]; dlt[u] = 0; in_t[u] = 0; }
    touched.clear();
    if (!ok) continue;
    adj[(size_t)a * n + b] = adj[(size_t)b * n + a] = 0; adj[(size_t)c * n + d] = adj[(size_t)d * n + c] = 0;
    adj[(size_t)a * n + d] = adj[(size_t)d * n + a] = 1; adj[(size_t)c * n + b] = adj[(size_t)b * n + c] = 1;
    rm_nb(a, b); rm_nb(b, a); rm_nb(c, d); rm_nb(d, c);
    nb[a].push_back(d); nb[d].push_back(a); nb[c].push_back(b); nb[b].push_back(c);
    A[e1] = a; B[e1] = d; A[e2] = c; B[e2] = b;
    T += dT; E += dE; sv = ns; acc++;
    if (E <= e_tol * Ets) reached = true;
  }
  IntegerMatrix out(m, 2);
  for (int e = 0; e < m; e++) { out(e, 0) = A[e] + 1; out(e, 1) = B[e] + 1; }
  NumericVector rel(M);
  for (int k = 0; k < M; k++) rel[k] = sv[k] / s0[k] - 1;
  return List::create(_["edges"] = out, _["T0"] = T0, _["T"] = T, _["accepted"] = (double)acc,
                      _["proposals"] = props, _["num_rel_change"] = rel,
                      _["E0_rel"] = E0 / Ets, _["E_rel"] = E / Ets, _["reached"] = reached);
}
