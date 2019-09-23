// File mipfp/src/ipfp_core.cpp
// by Oliver Krebs
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 or 3 of the License
//  (at your option).
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  http://www.r-project.org/Licenses/
//
// -----------------------------------------------------------------------------
// This file provides the function .IpfpCoreC which implements the core 
// algorithm of the iterative proportional fitting method using Rcpp. 
// The helper functions CalcIdxC, CalcStepIdxC, CalcIdxCoreC are used to 
// calculate the indices for applying a function over a specifc margin of a
// vectorized multidimensional array, i.e. the indices of each set of variables 
// of the given margin. 
// -----------------------------------------------------------------------------

#include <Rcpp.h>

// parallel support
#ifdef _OPENMP
#include <omp.h>
#endif

void CalcIdxCoreC(const std::vector<int>& target_margin,
                  const std::vector<int>& dim_seed,
                  const int& n_dims,
                  const int current_dim,
                  int& current_result,   
                  std::vector<int>& tmp_res,
                  std::vector<int>& res) {  
  // Calculate indices of a group in a vectorized multidimensional array defined
  // by a given margin - this is a recursive function that is not supposed to be
  // called directly by the user.
  //
  // Author: O. Krebs
  //  
  // Args:
  //   target_margin: The margin for which to obtain the indices.
  //   dim_seed: Dimensions of the seed for which to obtain indices.
  //   n_dims: Number of dimensions in the given margin
  //   current_dim: The currently visited dimension in the recursive process
  //   current_result: Current result in the recursive process
  //   tmp_res: Temporary result vector
  //   res: Final result vector
  //
  // Returns: No value. Changes the passed res and tmp_res vectors. 
  
  // sum temporary result vector if end of recursion is reached for the current
  // dim of the margin
  if (current_dim > n_dims - 1) {
    res[current_result] = 1;
    for (int n = 0; n < current_dim; n++) {
      res[current_result] += tmp_res[n];  
    }  
    current_result += 1;
  } else {
    // iterate over the current dimension in the array. "Backwards" iteration 
    // ensures compatibility with apply-style margins
    int i_max = dim_seed[target_margin[n_dims - 1 - current_dim] - 1];
    for (int i = 0; i < i_max; i++) {
      tmp_res[current_dim] = i;
      // step size is defined by all dimensions below the current one
      // -1 converts R style indices to c++ style indices
      for (int j = 0; j < target_margin[n_dims - 1 - current_dim] - 1; j++) {
        tmp_res[current_dim] *= dim_seed[j];
      }
      // continue recursion for the next dimension
      CalcIdxCoreC(target_margin, dim_seed, n_dims, current_dim + 1, 
                   current_result, tmp_res, res);
    }
  }
}

std::vector<int> CalcIdxC(const std::vector<int>& target_margin,
                          const std::vector<int>& dim_seed) {
  // Calculate indices for groups - defined by a given margin - in a vectorized 
  // multidimensional array. This function sets up necessary values and calls 
  // the recursive function CalcIdxCore above.
  //
  // Author: O. Krebs
  //  
  // Args:
  //   target_margin: The margin for which to obtain the indices.
  //   dim_seed: Dimensions of the seed for which to obtain indices.
  // Returns: Vector of indicies for groups - defined by a given margin - in a 
  //   vectorized multidimensional array. 
  
  // set up for recursive call
  int n_dims = target_margin.size();
  std::vector<int> tmp_res(n_dims);
  
  // get number of groups defined by the current margin
  int res_size = 1;
  for (int i = 0; i < n_dims; i++) {
    res_size *= dim_seed[target_margin[i] - 1];
  }
  std::vector<int> res(res_size);
  
  int current_result = 0;
  int current_dim = 0;
  
  // call recursive function which changes res
  CalcIdxCoreC(target_margin, dim_seed, n_dims, current_dim, current_result, 
               tmp_res, res);
  
  return(res);
}


std::vector<int> CalcStepIdxC(const std::vector<int> target_margin,
                              const std::vector<int>& dim_seed) {
  // Calculate strides between values in one group - defined by a given margin -
  // in a vectorized multidimensional array. 
  //
  // Author: O. Krebs
  //  
  // Args:
  //   target_margin: The margin for which to obtain the index steps.
  //   dim_seed: Dimensions of the seed for which to obtain index steps.
  // Returns: Vector of strides between values of one group - defined by a given
  //   margin - in a vectorized multidimensional array. 
  
  // To get strides run CalcIdx on all dimensions not in the target_margin
  std::vector<int> neg_margin;
  for (int d = 0, d_size = dim_seed.size(); d < d_size; d++) {
    if (std::find(target_margin.begin(), target_margin.end(), d + 1) 
          == target_margin.end()) {
      // keep dimension notation as in target_margin/R, i.e. starting from 1
      neg_margin.push_back(d + 1);
    }
  }
  return(CalcIdxC(neg_margin, dim_seed));
}

// [[Rcpp::export(name = ".IpfpCoreC")]]
Rcpp::List IpfpCoreC(const Rcpp::NumericVector& seed,
                     const Rcpp::List& target_list,
                     const Rcpp::List& target_data,
                     const bool& print,
                     const int& iter,
                     const double& tol,
                     const bool& na_target,
                     const int& n_threads = 1) {
  // Update an array using the iterative proportional fitting procedure.
  //
  // Author: O. Krebs
  //         Derived from function Ipfp by J. Barthelemy
  //  
  // Args:
  //   seed: The initial multi-dimensional array in vectorized form. Each cell 
  //         must be non-negative.
  //   target_list: A list of the target margins provided in target_data. Each
  //                component of the list is a (potentially vectorized) array 
  //                whose cells indicates which dimension the corresponding 
  //                margin relates to.
  //   target_data: A list containing the data of the target margins. Each
  //                component of the list is a (potentially vectorized) array 
  //                storing a margin. The list order must follow the one defined
  //                in target_list. Note that the cells of the arrays must be 
  //                non-negative, but may contain NA values.
  //   print: Verbose parameter: if TRUE prints the current iteration number
  //          and the value of the stopping criterion.
  //   iter: The maximum number of iteration allowed; must be greater than 0.
  //   tol: If the maximum absolute difference between two iteration is lower
  //        than the value specified by tol, then ipfp has reached convergence
  //        (stopping criterion); must be greater than 0.
  
  // cast target list and target data to vector
  const int n_targets = target_list.size();
  std::vector<std::vector<int> > vtarget;
  std::vector<Rcpp::NumericVector> vtarget_data; 
  for (int t = 0; t < n_targets; t++) {
    vtarget.push_back(Rcpp::as<std::vector<int> >(target_list[t]));
    vtarget_data.push_back(Rcpp::as<Rcpp::NumericVector>(target_data[t]));
  }
  
  // transform target margins to start and stride indices for iterating the seed
  const std::vector<int> dim_seed = 
    Rcpp::as<std::vector<int> >(seed.attr("dim"));
  std::vector<std::vector<int> > start_idx;
  std::vector<std::vector<int> > step_idx;
  for (int t = 0; t < n_targets; t++) {
    start_idx.push_back(CalcIdxC(vtarget[t], dim_seed));
    step_idx.push_back(CalcStepIdxC(vtarget[t], dim_seed));
  }
  
  // remove groups, i.e. start indices, for NA values in each target
  if (na_target == true) {
    std::vector<std::vector<int> > tmp(n_targets);
    for (int t = 0; t < n_targets; t++) {
      for (int i = 0, max_i = vtarget_data[t].size(); i < max_i; i++) {
        if(!Rcpp::NumericVector::is_na(vtarget_data[t][i]))
          tmp[t].push_back(start_idx[t][i]);
      }
      vtarget_data[t] = na_omit(vtarget_data[t]);
    }
    tmp.swap(start_idx);
  }
  
  // get the number of values in each constraint and largest constraint
  std::vector<int> t_length(n_targets, 1);
  for (int t = 0; t < n_targets; t++) {
    t_length[t] = vtarget_data[t].size();
  }
  
  // get deep copies of the seed to modify
  std::vector<double> result = Rcpp::as<std::vector<double> >(seed);
  std::vector<double> result_tmp = Rcpp::as<std::vector<double> >(seed);
  std::size_t result_length = result.size();
  
  // initializations for the loop
  double result_diff = 0;
  double stp_crit = 0;
  std::vector<double> evol_stp_crit;
  
  bool converged = false;
  int i = 0;
  
  while (converged == false && i < iter) {
    if (print == true) {
      Rcpp::Rcout << "... ITER " << i + 1 << std::endl;
    }
    
    // saving previous iteration result (for testing convergence)
    result.swap(result_tmp);
    
    // reset stp.crit
    stp_crit = 0;
    
    for(std::size_t t = 0; t < n_targets; t++) {
      
      // update result with target t
      #pragma omp parallel for num_threads(n_threads)
      for (std::size_t j = 0; j < t_length[t]; j++) {
        
        long double tmp_sum = 0;
        long double update_factor = 0;
        
        // get current margin sum
        for (std::size_t s = 0, max_s = step_idx[t].size(); s < max_s; s++) {
          // pulling result from result_temp on first run, avoids having to 
          // copy result_temp to result above
          if (t  == 0) {
            tmp_sum += result_tmp[start_idx[t][j] - 1 + step_idx[t][s] - 1];
          } else {
            tmp_sum += result[start_idx[t][j] - 1 + step_idx[t][s] - 1];
          }
        }
        // check for 0's and derive update factor
        if(tmp_sum == 0) {
          update_factor = -1;
        } else {
          update_factor = (vtarget_data[t][j] - tmp_sum) / tmp_sum;
        }
        // apply update factor
        for (std::size_t s = 0, max_s = step_idx[t].size(); s < max_s; s++) {
          if (t == 0) {
            result[start_idx[t][j] - 1 +  step_idx[t][s] - 1] = 
              result_tmp[start_idx[t][j] - 1 +  step_idx[t][s] - 1] + 
              result_tmp[start_idx[t][j] - 1 +  step_idx[t][s] - 1] * 
              update_factor;
          } else {
            result[start_idx[t][j] - 1 +  step_idx[t][s] - 1] +=
              result[start_idx[t][j] - 1 +  step_idx[t][s] - 1] *
              update_factor;
          }
        }
      }
      // when run on the last margin, check for convergence
      if(t == n_targets - 1) {
        stp_crit = 0;
        #pragma omp parallel for num_threads(n_threads) reduction(max:stp_crit)
        for (std::size_t j = 0; j < result_length; j++) {
          result_diff = std::abs(result[j] - result_tmp[j]);
          if (result_diff > stp_crit) {
            stp_crit = result_diff;
          }
        }
      }
    }
    if (print == true) {
      Rcpp::Rcout <<"       stoping criterion: " << stp_crit << std::endl;
    }
    if (stp_crit < tol) {
      converged = true;
      if (print == true) {
        Rcpp::Rcout << "Convergence reached after " << i + 1 << " iterations!" 
                    << std::endl;
      }
    }
    evol_stp_crit.push_back(stp_crit);
    i++;
  }
  Rcpp::NumericVector ret_result = Rcpp::wrap(result);
  ret_result.attr("dim") = dim_seed;
  return Rcpp::List::create(Rcpp::Named("x.hat") = ret_result,
                            Rcpp::Named("evol.stp.crit") = 
                              Rcpp::wrap(evol_stp_crit),
                            Rcpp::Named("conv") = Rcpp::wrap(converged));
}
