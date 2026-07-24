// ===============================================================
// dominance_swapper.cpp
// ===============================================================
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List allocate_dominance_flows(NumericVector pool_size, 
                              NumericVector hole_size, 
                              IntegerVector edge_p, 
                              IntegerVector edge_c) {
  
  int n_edges = edge_p.size();
  NumericVector edge_fill(n_edges, 0.0);
  
  // Create mutable copies of the pools and holes to track remainders
  NumericVector p_pool = clone(pool_size);
  NumericVector p_hole = clone(hole_size);
  
  for(int i = 0; i < n_edges; i++) {
    // R indexes are 1-based, convert to 0-based for C++ arrays
    int p_idx = edge_p[i] - 1;
    int c_idx = edge_c[i] - 1;
    
    double available_class = p_pool[c_idx];
    double available_hole = p_hole[p_idx];
    
    if (available_class > 0 && available_hole > 0) {
      // The flow is strictly the minimum of the remaining class pool and the remaining hole
      double flow = std::min(available_class, available_hole);
      
      edge_fill[i] = flow;
      
      // Drain the pool and hole
      p_pool[c_idx] -= flow;
      p_hole[p_idx] -= flow;
      
      // Fix floating point dust (if it's strictly under 1e-12, lock to 0)
      if (p_pool[c_idx] < 1e-12) p_pool[c_idx] = 0;
      if (p_hole[p_idx] < 1e-12) p_hole[p_idx] = 0;
    }
  }
  
  return List::create(
    Named("edge_fill") = edge_fill,
    Named("pool_remainder") = p_pool,
    Named("hole_remainder") = p_hole
  );
}
