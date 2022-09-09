////////////////////////////////////////////////////
// Description:
// Multispecies CPUE analysis (dynamics factor analysis)
// covariate effects are estimated for each species
// species interaction effect are estimated using DFA
// if model is based on "area" level, there is a need to incorporate average "targeting behavior"

#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// square
template<class Type>
Type square(Type x){
  return pow(x,2);
}

// Function to important barrier-SPDE code
template<class Type>
struct spde_barrier_t{
  vector<Type> C0;
  vector<Type> C1;
  Eigen::SparseMatrix<Type> D0;
  Eigen::SparseMatrix<Type> D1;
  Eigen::SparseMatrix<Type> I;
  spde_barrier_t(SEXP x){           // x = List passed from R
    C0 = asVector<Type>(getListElement(x,"C0"));
    C1 = asVector<Type>(getListElement(x,"C1"));
    D0 = tmbutils::asSparseMatrix<Type>(getListElement(x,"D0"));
    D1 = tmbutils::asSparseMatrix<Type>(getListElement(x,"D1"));
    I = tmbutils::asSparseMatrix<Type>(getListElement(x,"I"));
  }
};

// Function to calculate Q (precision) matrix using barrier-SPDE
template<class Type>
Eigen::SparseMatrix<Type> Q_spde1(spde_barrier_t<Type> spde_barrier, Type kappa, vector<Type> c){
  //using namespace Eigen;
  vector <Type> range(2);
  range(0) = sqrt(8.0)/kappa*c(0);
  range(1) = range(0)*c(1);
  Type pi = 3.141592;

  int dimLatent = spde_barrier.D0.row(0).size();
  vector<Type> Cdiag(dimLatent);
  Eigen::SparseMatrix<Type > Cinv(dimLatent,dimLatent);

  Cdiag = spde_barrier.C0*pow(range(0),2.0) + spde_barrier.C1*pow(range(1),2.0);
  for(int i =0; i<dimLatent; ++i){
    Cinv.coeffRef(i,i) = 1/Cdiag(i);
  }

  Eigen::SparseMatrix<Type>A = spde_barrier.I;
  A = A + (pow(range(0),2.0)/8.0) * spde_barrier.D0 + (pow(range(1),2.0)/8.0) * spde_barrier.D1;

  Eigen::SparseMatrix<Type> Q = A.transpose() * Cinv * A/pi *2 * 3;

  return Q;
}

// Repeat vector
template <class Type>
vector<Type> RepeatVector(vector<Type> x, int times)
{
  int n = x.size() * times;
  vector<Type> res(n);
  int k = 0;
  for (int i = 0; i < times; i++) {
    for (int j = 0; j < x.size(); j++) {
      res[k] = x(j);
      k++;
    }
  }
  return res;
}

// Parameter transform for the autocorrelation coefficient
// approach 1
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(0.5) * x)) - Type(1);}

// approach 2
template <class Type>
Type zerofive_to_one(Type x)
{
  return Type(1) - Type(0.5) * invlogit(x);
}


// some likelihood functions
template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0)
{
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0)
{
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template <class Type>
Type logistic_func(Type x, Type s50, Type slope) {
  // logistic function controling the density-dependent effect.
  // logistic function is similar to length or size based selectivity
  // in fisheries, parameterized by the points at which f(x) = 0.5 and the width of the slope
  Type pred = Type(1.0) / (Type(1.0) + exp(-log(Type(19.0)) * (x - s50) / slope));
  return pred;
}


// some specifications of available options
enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2,
  inverse_link  = 3
};

template <class Type>
Type InverseLink(Type eta, int link)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = exp(eta);
    break;
  case logit_link:
    out = invlogit(eta);
    break;
  case inverse_link:
    out = Type(1.0) / eta;
    break;
  default:
    error("Link not implemented.");
  }
  return out;
}

enum valid_family {
  tweedie_family  = 0
  // gaussian_family   = 1,
  // poisson_family    = 2,
  // gamma_family      = 3,
  // nb_family         = 4,
  // student_family    = 5,
  // lognormal_family  = 6
};

// Function to creating the loading matrix of Nspecies to Nfactor factors
template<class Type>
matrix<Type> create_loadings_matrix(vector<Type> L_val, int n_rows, int n_cols)
{
  matrix<Type> L_rc(n_rows, n_cols);
  int Count = 0;
  for(int r=0; r<n_rows; r++){
    for(int c=0; c<n_cols; c++){
      L_rc(r,c) = L_val(Count);
      Count++;
    }
  }
  return(L_rc);
};

// Generate loadings matrix for covariance
// zerosum_penalty -- used for EOF indices when also estimating Omega (such that EOF is zero-centered index)
// trace_sum_penalty -- used for sum of squared elements,


template<class Type>
matrix<Type> create_loadings_covariance( vector<Type> L_val, int n_rows, int n_cols){
  matrix<Type> L_rc(n_rows, n_cols);
  int Count = 0;
  for(int r=0; r<n_rows; r++){
    for(int c=0; c<n_cols; c++){
      if(r>=c){
        L_rc(r,c) = L_val(Count);
        Count++;
      }else{
        L_rc(r,c) = 0.0;
      }
    }}

  // // Zero-sum constraint
  // if( zerosum_penalty > 0 ){
  //   vector<Type> colsum( n_cols );
  //   colsum.setZero();
  //   for(int c=0; c<n_cols; c++){
  //     for(int r=0; r<n_rows; r++){
  //       colsum(c) += L_rc(r,c);
  //     }}
  //   for(int c=0; c<n_cols; c++){
  //     for(int r=0; r<n_rows; r++){
  //       L_rc(r,c) -= colsum(c) / n_rows;
  //     }}
  //   for(int c=0; c<n_cols; c++){
  //     jnll_pointer += zerosum_penalty * square(colsum(c));
  //   }
  // }
  // // Trace sum-to-one constraint
  // if( trace_sum_penalty > 0 ){
  //   Type Cov_trace = 0;
  //   for(int c=0; c<n_cols; c++){
  //     for(int r=0; r<n_rows; r++){
  //       Cov_trace += square(L_rc(r,c));
  //     }}
  //   for(int c=0; c<n_cols; c++){
  //     for(int r=0; r<n_rows; r++){
  //       L_rc(r,c) = L_rc(r,c) / sqrt(Cov_trace);  // L is the sqrt(Cov) so divide by sqrt(Cov_trace)
  //     }}
  //   for(int c=0; c<n_cols; c++){
  //     jnll_pointer += trace_sum_penalty * square(log(Cov_trace));
  //   }
  // }
  return L_rc;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;

  // Data section
  DATA_INTEGER(Nfactor);                     //The number of age groups
  DATA_MATRIX(X);                            //Design matrix for the fixed effects (same for each age groups)
  DATA_MATRIX(yobs);                          //the observed species abundance at each station

  // Random intercepts:
  DATA_IMATRIX(RE_indexes);
  DATA_IVECTOR(nobs_RE);
  DATA_IVECTOR(ln_tau_G_index);

  // Derived data
  int n_RE = nobs_RE.size();
  int Nspecies = yobs.cols();
  int Nobs = yobs.rows();

  // INLA features & model configurations (e.g. how to move from mesh to actual data point)
  // DATA_MATRIX_INDICATOR(keep, yobs);           // https://rdrr.io/cran/TMB/man/oneStepPredict.html
  // DATA_ARRAY(to_keep);                        // https://rdrr.io/cran/TMB/man/oneStepPredict.html

  // Distribution
  DATA_INTEGER(family);
  DATA_INTEGER(link);
  DATA_INTEGER(sim);
  DATA_INTEGER(incl_target);
  DATA_INTEGER(incl_sp_int);

  // Spatial model or not
  DATA_INTEGER(spatial_model);  // 0 = non spatial model, 1 = only spatial model without temporal aspect = constant spatial field
  DATA_INTEGER(include_st);     // 0 = no: only spatial model without temporal aspect = constant spatial field, 2 = IID, 3= AR1
  DATA_INTEGER(barrier);
  DATA_INTEGER(do_spatialDFA);
  DATA_STRUCT(spde, spde_t);
  // DATA_STRUCT(spde_barrier,spde_barrier_t);   		//Three matrices needed for representing the spatial field (+ 3 more for the barrier effect) 2007
  // DATA_VECTOR(barrier_scaling); // scaling of range


  // Things related to spatial model (fake info is included when using non spatial model)
  DATA_SPARSE_MATRIX(Aobs);                   //Matrix for interpolating points within triangles (for the observation) - used for the mean spatial field
  DATA_SPARSE_MATRIX(Ast);                    //Same but now divided for each year (each observation in a year) - used for the spatio-temporal field
  DATA_IVECTOR(A_spatial_index);              //Vector of stations to match up A_st output
  DATA_FACTOR(year_i);                        //Year index for the stations
  DATA_INTEGER(Nyear);                        //Number of years (needed if including the spatio-temporal effect)
  DATA_INTEGER(Nmesh);                        //Number of mesh (needed if including the spatio-temporal effect)
  DATA_INTEGER(Npred);                        //Number of prediction points

  // For the projection
  DATA_INTEGER(Do_predict);
  DATA_INTEGER(se_fit);
  DATA_MATRIX(X_proj);
  DATA_SPARSE_MATRIX(A_proj);
  DATA_IMATRIX(RE_indexes_proj);
  DATA_IVECTOR(A_spatial_index_proj);         //Vector of stations to match up A_st output

  //// Now defining all the parameters
  // Fixed effects
  PARAMETER_MATRIX(beta);                     // coefficient associated with the fixed effect: Neffect x Nspecies

  // Spatial random effect if any
  PARAMETER_MATRIX(omega);                    // The mean spatial effects: Nmesh x Nspecies
  PARAMETER_ARRAY(epsilon_st);                // The spatio-temporal effect: Nmesh x Ntime x Nspecies
  PARAMETER_VECTOR(transf_rho);               // The autocorrelation value between year
  PARAMETER_VECTOR(logKappa);                 // Spatial scale parameter in Matern covariance structures for Nspecies
  PARAMETER_VECTOR(logTauO);                  // Precision parameter for the average spatial field for Nspecies
  PARAMETER_VECTOR(logTauE);                  // Precision parameter for the spatio-temporal variation for Nspecies

	// Species interaction effect (random)
	PARAMETER_VECTOR(L_val_spatial);            // Values on the loading matrix for projecting the spatial effects from latent factors to species
	PARAMETER_VECTOR(L_val_target);             // Values on the loading matrix for the small scale targeting behavior
	PARAMETER_VECTOR(L_val_sp);             		// Values on the loading matrix needed to derive the CORR matrix
  PARAMETER_MATRIX(RE_species_latent);        // the latent targeting effect for each observation (Nobs X Nspecies)
  PARAMETER_MATRIX(RE_species_int);           // the species interaction effect for each observation (Nobs X Nspecies)
  PARAMETER_VECTOR(logsds);                   // the variance of the species interaction effect

  PARAMETER_VECTOR(thetaf);                   // tweedie only
  PARAMETER_VECTOR(ln_phi);                   // sigma / dispersion / etc.

    // Parameters for the random effect (acting on the interecept)
    PARAMETER_MATRIX(ln_tau_G);               // random intercept sigmas
    PARAMETER_MATRIX(RE);                     // random intercept deviations


  int Nsize = Nspecies;
  if (do_spatialDFA == 1) Nsize = Nfactor;
  // derived parameters
  vector<Type> rho(Nsize);
  for (int i=0; i<Nsize; i++){
    rho(i) =f(transf_rho(i));
  }

  vector<Type> sds=exp(logsds);

  // ======================== Calculate the linear predictors (link space) then back-transform ========================

  matrix<Type> eta(Nobs, Nspecies);     // this is at the link scale for each factor
  // matrix<Type> fixed(Nobs, fixed_species);   // this is at the fixed effect
  matrix<Type> mu(Nobs, Nspecies);      // this is at the real scale

  // Step 1: Add the contribution of the fixed effects to each species
  for(int j=0;j<Nspecies; j++){
    eta.col(j) = X * vector<Type>(beta.col(j));
  }

  // Step 2: Now add the contribution of the random effects when existing (on the intercept - no slope):
  if (nobs_RE(0) >0){
    for (int i = 0; i < Nobs; i++){
      for (int j = 0; j < Nspecies; j++){
        int temp = 0;
        for (int k = 0; k < n_RE; k++) {
          if (k == 0) eta(i,j) += RE(RE_indexes(i, k),j); // record it
          if (k > 0) {
            temp += nobs_RE(k-1);
            eta(i,j) += RE(RE_indexes(i, k) + temp, j); // record it
          }
        }
        // fixed(i,j) = eta(i,j);
      }
    }
  }

  matrix<Type> Lspatial = create_loadings_covariance(L_val_spatial, Nspecies, Nfactor);

  if (spatial_model == 1) {

  // Step 3: Now add the contribution of spatial effect if existing using a spatial DFA approach or NOT
    // Non spatial DFA approach = estimate for all species
      if (do_spatialDFA == 0) {

          // Step 3.1: the mean spatial field
          matrix<Type> omega_A(Nobs, Nspecies);
          for(int j=0;j<Nspecies; j++){
            omega_A.col(j) = Aobs* vector<Type>(omega.col(j));
            eta.col(j) = vector<Type>(eta.col(j)) + vector<Type>(omega_A.col(j));
          }

          // Step 3.2: the spatio-temporal random field
          if (include_st > 0) {

            array<Type> epsilon_st_A(Ast.rows(), Nyear, Nspecies);
            array<Type> epsilon_st_A_temp(Ast.rows(), Nyear);
            vector<Type> Tmp_st(Nmesh);
            for (int j = 0; j < Nspecies; j++){
              epsilon_st_A_temp.setZero();
              for (int i = 0; i < Nyear; i++){
                Tmp_st.setZero();
                for (int k=0; k < Nmesh; k++) {
                  Tmp_st(k) = epsilon_st(k,i,j);
                }
                epsilon_st_A_temp.col(i) = Ast * Tmp_st;
                for (int k=0; k < Ast.rows(); k++) {
                  epsilon_st_A(k,i,j) = epsilon_st_A_temp(k,i);
                }
              }
            }
            // now finding the YEAR the observation takes place and assign that effect value from epsilon_st_A
            matrix<Type> epsilon_st_A_mat(Nobs, Nspecies);
            for (int i = 0; i < Nobs; i++){
              for (int j = 0; j < Nspecies; j++){
                epsilon_st_A_mat(i,j) = epsilon_st_A(A_spatial_index(i), year_i(i),j);
                eta(i,j) += epsilon_st_A_mat(i,j);
              }
            }
          }


        }

    // Spatial DFA approach = estimate for all species
      if (do_spatialDFA == 1) {

          // Step 3.1: calculate contribution of mean spatial field
            matrix<Type> omega_A(Nobs, Nfactor);
            for(int j=0;j<Nfactor; j++){
              omega_A.col(j) = Aobs* vector<Type>(omega.col(j));
            }

          // now projecting omega_A to each species using the Latent factors
            for (int i = 0; i < Nobs; i++){
              for(int j=0; j<Nspecies; j++){
                for(int k=0; k<Nfactor; k++){
                  eta(i,j) += Lspatial(j,k)*omega_A(i,k);
                }
              }
            }

          // Step 3.2: the spatio-temporal random effect
            if (include_st > 0) {

              array<Type> epsilon_st_A(Ast.rows(), Nyear, Nfactor);
              array<Type> epsilon_st_A_temp(Ast.rows(), Nyear);
              vector<Type> Tmp_st(Nmesh);
              for (int j = 0; j < Nfactor; j++){
                epsilon_st_A_temp.setZero();
                for (int i = 0; i < Nyear; i++){
                  Tmp_st.setZero();
                  for (int k=0; k < Nmesh; k++) {
                    Tmp_st(k) = epsilon_st(k,i,j);
                  }
                  epsilon_st_A_temp.col(i) = Ast * Tmp_st;
                  for (int k=0; k < Ast.rows(); k++) {
                    epsilon_st_A(k,i,j) = epsilon_st_A_temp(k,i);
                  }
                }
              }
              // now finding the YEAR the observation takes place and assign that effect value from epsilon_st_A
              matrix<Type> epsilon_st_A_mat(Nobs, Nfactor);
              for (int i = 0; i < Nobs; i++){
                for (int j = 0; j < Nfactor; j++){
                  epsilon_st_A_mat(i,j) = epsilon_st_A(A_spatial_index(i), year_i(i),j);
                }
              }

              // Now projecting epsilon_st_A_mat to each species using the Latent factors
              for (int i = 0; i < Nobs; i++){
                for(int j=0; j<Nspecies == 0; j++){
                  for(int k=0; k<Nfactor; k++){
                    eta(i,j) += Lspatial(j,k)*epsilon_st_A_mat(i,k);
                  }
                }
              }

              }

        }

  } // end of spatial model type statement

  // Step 4: Now add the contribution of the fine scale (observation-level) species correlation (when catching them), if specified
  matrix<Type> Lt = create_loadings_covariance(L_val_target, Nspecies, Nfactor);
  if (incl_target == 1){
    for (int i = 0; i < Nobs; i++){
      for(int j=0; j<Nspecies; j++){
         for(int k=0; k<Nfactor; k++){
			     eta(i,j) += Lt(j,k)*RE_species_latent(i,k);
  		  }
      }
    }
  }

  // Step 5: Now add the contribution of species interaction, if specified
  if (incl_sp_int == 1){
    for (int i = 0; i < Nobs; i++){
      for(int j= 0; j < Nspecies; j++){
		     eta(i,j) += RE_species_int(i,j);
      }
    }
  }


  // Step 6: back transform to the real scale
  for (int i = 0; i < Nobs; i++){
    for (int j = 0; j < Nspecies; j++){
      mu(i,j) = InverseLink(eta(i,j), link);
    }
  }


  // ======================== The likelihood components ========================

  // Defining the NLL
  Type NLL = 0;

  // The random effects (only add the contribution if existing)
  if (nobs_RE(0) >0){
    for (int j = 0; j < Nspecies; j++){
      for (int k = 0; k < n_RE; k++){
        for (int g = 0; g < nobs_RE(k); g++) {
          if (k==0){
            NLL -= dnorm(RE(g,j), Type(0.0), exp(ln_tau_G(ln_tau_G_index(g),j)), true);
            if (sim == 1) SIMULATE{RE(g,j) = rnorm(Type(0), exp(ln_tau_G(ln_tau_G_index(g),j)));}
          }
          if (k>0){
            int temp4 =0;
            for (int bb=0; bb<(k-1); bb++) temp4 += nobs_RE(bb);
            NLL -= dnorm(RE(g+temp4,j), Type(0.0), exp(ln_tau_G(ln_tau_G_index(g+temp4),j)), true);
            if (sim == 1) SIMULATE{RE(g+temp4,j) = rnorm(Type(0), exp(ln_tau_G(ln_tau_G_index(g+temp4),j)));}
          }
        }
      }
    }
  }

  // Now add the likelihood contribution of the spatial random effects

    // The contribution of the average spatial field
    if (spatial_model == 1) {
      Eigen::SparseMatrix<Type> Q;

      if (do_spatialDFA == 0) {
        for(int i =0;i <Nspecies; i++){
          // if(barrier==1) {
            // Q = Q_spde1(spde_barrier, exp(logKappa(i)), barrier_scaling);
          // } else {
            Q = R_inla::Q_spde(spde, exp(logKappa(i)));
          // }
          NLL += SCALE(GMRF(Q,false), 1.0/exp(logTauO(i)))(omega.col(i));
        }
      }

      if (do_spatialDFA == 1) {
        for(int i =0;i <Nfactor; i++){
          // if(barrier==1) {
            // Q = Q_spde1(spde_barrier, exp(logKappa(i)), barrier_scaling);
          // } else {
            Q = R_inla::Q_spde(spde, exp(logKappa(i)));
          // }
          NLL += SCALE(GMRF(Q,false), 1.0/exp(logTauO(i)))(omega.col(i));
        }
      }

    // The contribution of the spatio-temporal random field
    if (include_st > 0) {

      array<Type> epsilon_st_temp1(Nmesh, Nyear);
      if (do_spatialDFA == 0) {
        for (int i = 0; i<Nspecies; i++){
          // if (barrier==1) {
            // Q = Q_spde1(spde_barrier, exp(logKappa(i)), barrier_scaling);
          // } else {
            Q = R_inla::Q_spde(spde, exp(logKappa(i)));
          // }
          epsilon_st_temp1.setZero();
          for(int t=0; t<Nyear; t++){
            for (int j=0; j<Nmesh; j++){
              epsilon_st_temp1(j,t) = epsilon_st(j,t,i);
            }
          }
          if (include_st == 2) NLL += SCALE(SEPARABLE(AR1(rho(i)), GMRF(Q,false)), 1.0/exp(logTauE(i)))(epsilon_st_temp1);
          if (include_st == 1) {
            for(int t=0; t<Nyear; t++){
              NLL += SCALE(GMRF(Q,false), 1.0/exp(logTauE(i)))(epsilon_st_temp1.col(t));
            }
          }

          if (sim == 1) {
            vector<Type> tmp(epsilon_st_temp1.rows());
            if (include_st == 2) {
              SIMULATE {SEPARABLE(AR1(rho(i)), GMRF(Q, false)).simulate(epsilon_st_temp1);}
              epsilon_st_temp1 *= 1./exp(logTauE(i));
            }
            if (include_st == 1){
              SIMULATE {GMRF(Q, false).simulate(tmp);}
              for(int t=0; t<Nyear; t++){
                epsilon_st_temp1.col(t) = tmp/exp(logTauE(i));
              }
            }
            for(int t=0; t<Nyear; t++){
              for (int j=0; j<Nmesh; j++){
                epsilon_st(j,t,i) = epsilon_st_temp1(j,t) ;
              }
            }
          }
        }
      }

      if (do_spatialDFA == 1) {
        for (int i = 0; i<Nfactor; i++){
          // if (barrier==1) {
            // Q = Q_spde1(spde_barrier, exp(logKappa(i)), barrier_scaling);
          // } else {
            Q = R_inla::Q_spde(spde, exp(logKappa(i)));
          // }
          epsilon_st_temp1.setZero();
          for(int t=0; t<Nyear; t++){
            for (int j=0; j<Nmesh; j++){
              epsilon_st_temp1(j,t) = epsilon_st(j,t,i);
            }
          }
          if (include_st == 2) NLL += SCALE(SEPARABLE(AR1(rho(i)), GMRF(Q,false)), 1.0/exp(logTauE(i)))(epsilon_st_temp1);
          if (include_st == 1) {
            for(int t=0; t<Nyear; t++){
              NLL += SCALE(GMRF(Q,false), 1.0/exp(logTauE(i)))(epsilon_st_temp1.col(t));
            }
          }

          if (sim == 1) {
            vector<Type> tmp(epsilon_st_temp1.rows());
            if (include_st == 2) {
              SIMULATE {SEPARABLE(AR1(rho(i)), GMRF(Q, false)).simulate(epsilon_st_temp1);}
              epsilon_st_temp1 *= 1./exp(logTauE(i));
            }
            if (include_st == 1){
              SIMULATE {GMRF(Q, false).simulate(tmp);}
              for(int t=0; t<Nyear; t++){
                epsilon_st_temp1.col(t) = tmp/exp(logTauE(i));
              }
            }
            for(int t=0; t<Nyear; t++){
              for (int j=0; j<Nmesh; j++){
                epsilon_st(j,t,i) = epsilon_st_temp1(j,t) ;
              }
            }
          }
        }
      }

    }

    }


    // The fine scale targetting behavior effect. Modeling latent variable effects to the species cpue
    if (incl_target == 1){

      for (int i=0; i<Nobs; i++) {
        for (int j=0; j<Nfactor; j++) {
          NLL -= dnorm(RE_species_latent(i,j), Type(0), Type(1), true);
          if (sim == 1) {
        		 SIMULATE{RE_species_latent(i,j) = rnorm(Type(0), Type(1));}
          }
        }
      }
    }

    // Now adding the possible effect of species interaction (co-occurence)
    matrix<Type> Cov(Nspecies,Nspecies);
    if (incl_sp_int == 1){
  		vector<Type> L_val_transf(L_val_sp.size());
      for (int ii =0; ii<L_val_sp.size(); ii++) {
        L_val_transf(ii) = f(L_val_sp(ii));
      }
      density::UNSTRUCTURED_CORR_t<Type> nldens(L_val_transf);
      density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sds);
      Cov = nldens.cov();

  		for (int i=0; i<Nobs; i++) {
  			NLL += scnldens(RE_species_int.row(i));
  			if (sim == 1) {
  			  RE_species_int.row(i) = sds * nldens.simulate();
  			}
  		}
    }

   // The observation likelihood: Only tweedie at the moment
    vector<Type> phi = exp(ln_phi);
    for (int i=0; i<Nobs; i++) {
      for (int j = 0; j < Nspecies; j++){
        // case tweedie_family:
        Type s1 = invlogit(thetaf(j)) + Type(1.0);
        // if (!isNA(yobs(i,j))) NLL -= to_keep(i,j) * keep(i,j) * dtweedie(yobs(i,j), mu(i,j), phi(j), s1, true);
        if (!isNA(yobs(i,j))) NLL -= dtweedie(yobs(i,j), mu(i,j), phi(j), s1, true);
        if (sim == 1) {
          SIMULATE { yobs(i,j) = rtweedie(mu(i,j), phi(j), s1);}
        }   // break;
      }
    }

    // Adding the penalty term related to the latent factor approach (based on the newest paper from)


  // ================= Creating the indices of abundance ====================

  if (Do_predict == 1) {
		int Npred_all = Npred*Nyear;
		matrix<Type> eta_proj(Npred_all, Nspecies);     // this is at the link scale for each factor
		matrix<Type> mu_proj(Npred_all, Nspecies);      // this is at the real scale

		// Step 1: Add the contribution of the fixed effects to each species
		for(int j=0;j<Nspecies; j++){
			eta_proj.col(j) = X_proj * vector<Type>(beta.col(j));
		}

		// Step 2: Add the contribution of the spatial effects
		if (spatial_model == 1) {

		  if (do_spatialDFA == 0) {
		    array<Type> epsilon_st_A_proj(Npred, Nyear, Nspecies);
    		array<Type> epsilon_st_A_proj_temp(Npred, Nyear);
    		vector<Type> Tmp_st_proj(Nmesh);
    		matrix<Type> omega_A_proj(Npred, Nspecies);

    		// the average spatial field
    		for(int j=0;j<Nspecies; j++){
    		  omega_A_proj.col(j) = A_proj * vector<Type>(omega.col(j));
    		}
   	 		for (int i = 0; i < Npred; i++){
    		  for(int j=0;j<Nspecies; j++){
    		    for(int k=0;k<Nyear; k++){
    		      eta_proj(i+Npred*k,j) += omega_A_proj(i,j);   // same value of omega each year
    		    }
    		  }
    		}

   	 		// the spatio-temporal field
   	 		if (include_st > 0) {
   	 		  for (int j = 0; j < Nspecies; j++){
   	 		    epsilon_st_A_proj_temp.setZero();
   	 		    for (int i = 0; i < Nyear; i++){
   	 		      Tmp_st_proj.setZero();
   	 		      for (int k=0; k < Nmesh; k++) {
   	 		        Tmp_st_proj(k) = epsilon_st(k,i,j);
   	 		      }
   	 		      epsilon_st_A_proj_temp.col(i) = A_proj * Tmp_st_proj;
   	 		      for (int k=0; k < Npred; k++) {
   	 		        epsilon_st_A_proj(k,i,j) = epsilon_st_A_proj_temp(k,i);
   	 		      }
   	 		    }
   	 		  }

   	 		  // allocating the annual projection to each row of prediction data.frame
   	 		  Type temp2 = Type(0);
   	 		  for (int ii = 0; ii < Npred; ii++){
   	 		    for (int jj = 0; jj < Nspecies; jj++){
   	 		      for (int kk=0; kk < Nyear; kk++){
   	 		        int index = ii+Npred*kk;
   	 		        temp2 = epsilon_st_A_proj(A_spatial_index_proj(index), kk, jj);
   	 		        eta_proj(index,jj) += temp2;
   	 		      }
   	 		    }
   	 		  }
   	 		}

		  }

		  if (do_spatialDFA == 1) {
		    array<Type> epsilon_st_A_proj(Npred, Nyear, Nfactor);
    		array<Type> epsilon_st_A_proj_temp(Npred, Nyear);
    		vector<Type> Tmp_st_proj(Nmesh);
    		matrix<Type> omega_A_proj(Npred, Nfactor);

    		// the average spatial field
    		for(int j=0;j<Nfactor; j++){
    		  omega_A_proj.col(j) = A_proj * vector<Type>(omega.col(j));
    		}
   	 		for (int i = 0; i < Npred; i++){
    		  for(int j=0;j<Nspecies; j++){
    		    for(int k=0;k<Nyear; k++){
    		      for(int l=0;l<Nfactor; l++){
    		        eta_proj(i+Npred*k,j) += Lspatial(j,l)*omega_A_proj(i,l);   // same value of omega each year
    		      }
    		    }
    		  }
    		}

   	 		// the spatio-temporal field
   	 		if (include_st > 0) {
   	 		  for (int j = 0; j < Nfactor; j++){
   	 		    epsilon_st_A_proj_temp.setZero();
   	 		    for (int i = 0; i < Nyear; i++){
   	 		      Tmp_st_proj.setZero();
   	 		      for (int k=0; k < Nmesh; k++) {
   	 		        Tmp_st_proj(k) = epsilon_st(k,i,j);
   	 		      }
   	 		      epsilon_st_A_proj_temp.col(i) = A_proj * Tmp_st_proj;
   	 		      for (int k=0; k < Npred; k++) {
   	 		        epsilon_st_A_proj(k,i,j) = epsilon_st_A_proj_temp(k,i);
   	 		      }
   	 		    }
   	 		  }

   	 		  // allocating the annual projection to each row of prediction data.frame
   	 		  Type temp2 = Type(0);
   	 		  for (int ii = 0; ii < Npred; ii++){
   	 		    for (int jj = 0; jj < Nspecies; jj++){
   	 		      for (int kk=0; kk < Nyear; kk++){
   	 		        for(int ll=0;ll<Nfactor; ll++){
   	 		          int index = ii+Npred*kk;
   	 		          eta_proj(index,jj) += Lspatial(jj,ll)*epsilon_st_A_proj(A_spatial_index_proj(index), kk, ll);
   	 		        }
   	 		      }
   	 		    }
   	 		  }
   	 		}
		  }


		}


		// Step 3: Now add the contribution of the random effects when existing (on the intercept - no slope):
		// if the random effect if not present, its takes automatically a 0 value
		if (nobs_RE(0) >0){
			for (int i = 0; i < Npred_all; i++){
				for (int j = 0; j < Nspecies; j++){
					int temp = 0;
					for (int k = 0; k < n_RE; k++) {
						if (k == 0) {
							if(RE_indexes_proj(i, k) != 999) eta_proj(i,j) += RE(RE_indexes_proj(i, k),j); // record it
						}
						if (k > 0) {
							temp += nobs_RE(k-1);
							if(RE_indexes_proj(i, k) != 999) eta_proj(i,j) += RE(RE_indexes_proj(i, k) + temp, j); // record it
						}
					}
					// fixed(i,j) = eta(i,j);
				}
			}
		}


//   if (incl_target == 1){
//  		for (int i = 0; i < Npred; i++){
// 		  for(int j=0; j<Nspecies; j++){
// 		    for(int k=0; k<Nfactor; k++){
// 		      eta_proj(i,j) += Lt(j,k)*RE_species_latent(target_index,k);
// 		    }
// 		  }
// 		}
// 	}

		// Step 5: back transform to the real scale
		for (int i = 0; i < Npred_all; i++){
			for (int j = 0; j < Nspecies; j++){
				mu_proj(i,j) = InverseLink(eta_proj(i,j), link);
			}
		}

		if (se_fit == 0) REPORT(mu_proj);
		if (se_fit == 1) ADREPORT(mu_proj);

  }

  // ======================== Reporting ========================

  // Model parameters
  REPORT(beta);

  if (sim == 1) {
    SIMULATE {
    REPORT(yobs);
    REPORT(RE);
    REPORT(RE_species_latent);
    }
  }
  // Some outputs for the diagnostics
  REPORT(mu);
  REPORT(Lt);
  REPORT(RE_species_latent);
  if (incl_sp_int == 1) {
    REPORT(Cov);
    REPORT(sds);
  }

  return NLL;

}

