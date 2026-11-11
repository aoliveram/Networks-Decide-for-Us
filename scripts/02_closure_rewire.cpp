// Constrained rewiring toward more triadic closure (scripts/02_network_calibrated.R).
// Loaded with Rcpp::sourceCpp from that script. Two targets:
//   rewire_closure_node()  every node toward target_factor x its own triangles: level x2 (main)
//   rewire_closure()       the TOTAL triangles toward target_factor x: level x2g (validation only)
//
// rewire_closure():
// Swap: edges (a,b),(c,d) -> (a,d),(c,b). Every node keeps its degree. The swap is allowed only if,
// for race, sex and religion separately, x(b) == x(d) or x(a) == x(c): then the multiset of attribute
// pairs over edges is unchanged, so the race mixing matrix (and same-sex, same-religion counts) stay
// exactly as they were. The total |age| and |educ| differences over edges must stay within `tol`
// (relative) of their starting values. Moves are accepted toward a target number of triangles
// (Metropolis on |T - T*|, temperature `temp`); neutral moves are accepted. After T* is reached, a
// mixing phase of `mix_props` more proposals keeps |T - T*| small while the network keeps changing.
// With degrees fixed, the number of connected triples is fixed, so transitivity scales with T.
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List rewire_closure(IntegerMatrix el, int n, IntegerVector race, IntegerVector sex, IntegerVector relig,
                    NumericVector age, NumericVector educ, double target_factor, double tol,
                    double max_props, double temp, double mix_props) {
  const int m = el.nrow();
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
  double T = 0, sa0 = 0, se0 = 0;
  for (int e = 0; e < m; e++) {
    T += cn(A[e], B[e]);
    sa0 += std::fabs(age[A[e]] - age[B[e]]); se0 += std::fabs(educ[A[e]] - educ[B[e]]);
  }
  T /= 3.0;
  const double T0 = T, Tstar = target_factor * T0;
  double sa = sa0, se = se0, props = 0;
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
    if (!(race[b] == race[d] || race[a] == race[c])) continue;
    if (!(sex[b] == sex[d] || sex[a] == sex[c])) continue;
    if (!(relig[b] == relig[d] || relig[a] == relig[c])) continue;
    double nsa = sa - std::fabs(age[a] - age[b]) - std::fabs(age[c] - age[d]) + std::fabs(age[a] - age[d]) + std::fabs(age[c] - age[b]);
    if (std::fabs(nsa - sa0) > tol * sa0) continue;
    double nse = se - std::fabs(educ[a] - educ[b]) - std::fabs(educ[c] - educ[d]) + std::fabs(educ[a] - educ[d]) + std::fabs(educ[c] - educ[b]);
    if (std::fabs(nse - se0) > tol * se0) continue;
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
    T += dT; sa = nsa; se = nse; acc++;
    if (T >= Tstar) reached = true;
  }
  IntegerMatrix out(m, 2);
  for (int e = 0; e < m; e++) { out(e, 0) = A[e] + 1; out(e, 1) = B[e] + 1; }
  return List::create(_["edges"] = out, _["T0"] = T0, _["T"] = T, _["accepted"] = (double)acc,
                      _["proposals"] = props, _["age_rel_change"] = sa / sa0 - 1, _["educ_rel_change"] = se / se0 - 1);
}

// Variant "node": the same swaps and constraints, but the target is per node: every node i
// moves toward target_factor x its own starting triangles t_i(0). Energy E = sum_i |t_i - t_i*|
// (Metropolis, temperature `temp`). The global target of rewire_closure() lets the triangles go
// where they are cheapest, i.e. among high-degree nodes; here the local clustering of every node
// is multiplied by the same factor, so the clustering-degree profile of the starting network is
// kept. Phase 1 ends when E <= e_tol x sum_i t_i*; then `mix_props` proposals as above.
// [[Rcpp::export]]
List rewire_closure_node(IntegerMatrix el, int n, IntegerVector race, IntegerVector sex, IntegerVector relig,
                         NumericVector age, NumericVector educ, double target_factor, double tol,
                         double max_props, double temp, double mix_props, double e_tol) {
  const int m = el.nrow();
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
  double Ets = 0, E = 0, T = 0, sa0 = 0, se0 = 0;
  for (int i = 0; i < n; i++) { ts[i] = target_factor * t[i]; Ets += ts[i]; E += std::fabs(t[i] - ts[i]); T += t[i]; }
  T /= 3.0;
  const double T0 = T, E0 = E;
  for (int e = 0; e < m; e++) { sa0 += std::fabs(age[A[e]] - age[B[e]]); se0 += std::fabs(educ[A[e]] - educ[B[e]]); }
  double sa = sa0, se = se0, props = 0;
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
    if (!(race[b] == race[d] || race[a] == race[c])) continue;
    if (!(sex[b] == sex[d] || sex[a] == sex[c])) continue;
    if (!(relig[b] == relig[d] || relig[a] == relig[c])) continue;
    double nsa = sa - std::fabs(age[a] - age[b]) - std::fabs(age[c] - age[d]) + std::fabs(age[a] - age[d]) + std::fabs(age[c] - age[b]);
    if (std::fabs(nsa - sa0) > tol * sa0) continue;
    double nse = se - std::fabs(educ[a] - educ[b]) - std::fabs(educ[c] - educ[d]) + std::fabs(educ[a] - educ[d]) + std::fabs(educ[c] - educ[b]);
    if (std::fabs(nse - se0) > tol * se0) continue;
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
    T += dT; E += dE; sa = nsa; se = nse; acc++;
    if (E <= e_tol * Ets) reached = true;
  }
  IntegerMatrix out(m, 2);
  for (int e = 0; e < m; e++) { out(e, 0) = A[e] + 1; out(e, 1) = B[e] + 1; }
  return List::create(_["edges"] = out, _["T0"] = T0, _["T"] = T, _["accepted"] = (double)acc,
                      _["proposals"] = props, _["age_rel_change"] = sa / sa0 - 1, _["educ_rel_change"] = se / se0 - 1,
                      _["E0_rel"] = E0 / Ets, _["E_rel"] = E / Ets, _["reached"] = reached);
}
