#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace density;



template<class Type>
Type dcrwSSM(objective_function<Type> * obj) {

  // Input data
  DATA_ARRAY(y); //track coordinates
  DATA_IVECTOR(b); //vector of behavioural states, treated as known
  DATA_IVECTOR(idx); // vector of the intervals of y in xhat
  DATA_VECTOR(jidx); // proportion of the time interval of y in xhat
  DATA_MATRIX(ae); //matrix nx4, all of the parameters for the t-distributions for lon and lat

  // Input parameters - i.e. parameters to estimate.
  PARAMETER_VECTOR(working_theta);
  PARAMETER_VECTOR(working_gamma); // First autocorrelation - logit because 0 < gamma < 1

  PARAMETER(working_sigma_lon); // Process standard deviation in lon - log because sdlon > 0
  PARAMETER(working_sigma_lat); // Process standard deviation in lat - log because sdlat > 0

  PARAMETER_MATRIX(working_A); //matrix of switching probabilities

  PARAMETER(working_psi); //working value for measurement error

  PARAMETER_MATRIX(x);

  // Transformation of the input parameters to model format
  vector<Type> theta(2);
  theta(0) = 2*M_PI/(1.0+exp(-working_theta(0))) - M_PI; // -pi to pi
  theta(1) = 2*M_PI/(1.0+exp(-working_theta(1))); // 0 to 2pi
  //vector<Type> theta = 2*M_PI/(1.0+exp(-logitTheta1)) - M_PI; // -pi to pi

  // vector<Type> gamma = working_gamma; //1.0/(1.0+exp(-working_gamma));
//  vector<Type> gamma = 1.0/(1.0+exp(-working_gamma));
  vector<Type> gamma(2);
  gamma(0) = 1.0/(1.0+exp(-working_gamma(0)));
  gamma(1) = gamma(0)/(1.0+exp(-working_gamma(1)));

  Type sigma_lon = exp(working_sigma_lon);
  Type sigma_lat = exp(working_sigma_lat);
  Type psi = exp(working_psi);

  //setting up the matrix of switching probabilities
  int m = working_A.rows(); // number of states
  matrix<Type> A(m, m);
  for(int i = 0; i < m; ++i){
      for(int j = 0; j < (m-1); ++j){
          A(i,j) = exp(working_A(i,j));  //ensures the probabilities we estimate are >0
      }
      A(i,m-1) = 1.0;
      A.row(i) *= 1.0/A.row(i).sum();
  }

  // System of equations for finding stationary distribution:    'I - Gamma + 1'
  matrix<Type> system(m, m);
  for(int i=0; i<m; i++) {
      for(int j=0; j<m; j++) {
          if( i ==j ) {
              system(i,j) = 1 - A(j,i) + 1;
          } else {
              system(i,j) = 0 - A(j,i) + 1;
          }
      }
  }

  // vector of ones for finding stationary distribution
  vector<Type> ones(m);
  for(int i=0; i<m; i++) {
      ones(i) = 1.0;
  }
  // take the inverse of the system matrix and then multiply by the vector of ones
  matrix<Type> sys_inverse = system.inverse();
  vector<Type> delta = sys_inverse*ones;

  // Report the parameters and their standard errors in their model format
  REPORT(theta);
  REPORT(gamma);
  REPORT(sigma_lon);
  REPORT(sigma_lat);
  REPORT(psi);
  REPORT(A);
  ADREPORT(theta);
  ADREPORT(gamma);
  ADREPORT(sigma_lon);
  ADREPORT(sigma_lat);
  ADREPORT(psi);


 int n_obs = y.cols();
 int n_x = x.cols();




    // Variance-covariance matrix for the process equation.
    matrix<Type> Sigma(2,2);
    Sigma << sigma_lon*sigma_lon, 0.0,
    0.0, sigma_lat*sigma_lat;
    MVNORM_t<Type> nll_proc(Sigma);

    //Calculate logLikelihood
    Type nll = 0.0;

    vector<Type> tmp(2);
    vector<Type> tmp2(2);
    matrix<Type> T(2,2); //rotational matrix
    vector<Type> xhat(2);


    // behavioural state markov chain
    // parameters (A) should always be fixed, but the constant is important to have
    // in the likelihood for determining the best model

    Type markovchain = log(delta(b(0)));
    for(int i=1; i < n_x; ++i){
        markovchain += log(A(b(i-1), b(i)));
        // the A matrix is the previous state in rows, the next state in columns
        // so the likelihood of this observed markov chain is just the right entries of the matrix added together
    }
    REPORT(markovchain);
    nll -= markovchain;

    //  report the likelihood values to troubleshoot
    vector<Type> pe_vec(n_x);
    vector<Type> me_vec_lon(n_obs);
    vector<Type> me_vec_lat(n_obs);


    // Initial process contribution

    tmp = x.col(1)-x.col(0);
    nll += nll_proc(tmp);
    pe_vec(1) = nll_proc(tmp);

    //movement (process) equation
    for(int i = 2; i < n_x; ++i){

        T << cos(Type(theta(b(i-1)))), -sin(Type(theta(b(i-1)))),
        sin(Type(theta(b(i-1)))), cos(Type(theta(b(i-1)))); //might need to put types on this?
        tmp = (x.col(i)-x.col(i-1)) - T*(x.col(i-1) - x.col(i-2))*Type(gamma(b(i)));

        nll += nll_proc(tmp);
        pe_vec(i) = nll_proc(tmp);

    }

    //measurement equation
    for(int i = 0; i < n_obs; ++i){

        xhat = (1.0-jidx(i)) * x.col(idx(i)-1) + jidx(i)*x.col(idx(i));
        tmp = y.col(i);
        tmp2 = tmp-xhat;
        // nll += nll_meas(tmp-xhat); //MVNORM_t returns the negative log likelihood
        // taken from Marie's DCRW_Argos.cpp file
        // based on this pretty sure dt returns the log likelihood
        nll -= ( 0.5*log(psi) - log(ae(i,0)) + dt(sqrt(psi)*tmp2(0)/ae(i,0),ae(i,1),true) ); // Longitude
        me_vec_lon(i) = -(0.5*log(psi) - log(ae(i,0)) + dt(sqrt(psi)*tmp2(0)/ae(i,0),ae(i,1),true)); // Longitude
        nll -= ( 0.5*log(psi) - log(ae(i,2)) + dt(sqrt(psi)*tmp2(1)/ae(i,2),ae(i,3),true) ); // Latitude
        me_vec_lat(i) = -(0.5*log(psi) - log(ae(i,2)) + dt(sqrt(psi)*tmp2(1)/ae(i,2),ae(i,3),true)); // Latitude
    }

    Type nll2 = pe_vec.sum() + me_vec_lon.sum() + me_vec_lat.sum();
    REPORT(nll2);
    REPORT(nll);
    REPORT(pe_vec);
    REPORT(me_vec_lon);
    REPORT(me_vec_lat);

  return nll2;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
