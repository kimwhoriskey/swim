#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace density;

template<class Type>
struct track {
  array<Type> y;
  vector<int> b;
  vector<int> idx;
  vector<Type> jidx;
  matrix<Type> ae;
  track(SEXP x) {
    y = tmbutils::asArray<Type>(getListElement(x, "y"));
    b = asVector<int>(getListElement(x, "b"));
    idx = asVector<int>(getListElement(x, "idx"));
    jidx = asVector<Type>(getListElement(x, "jidx"));
    ae = asMatrix<Type>(getListElement(x, "ae"));
  }
};


template<class Type>
Type dcrwSSM(objective_function<Type> * obj) {

  // Input data
  DATA_STRUCT(dat, track);

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


  // Variance-covariance matrix for the process equation.
  matrix<Type> Sigma(2,2);
  Sigma << sigma_lon*sigma_lon, 0.0,
  0.0, sigma_lat*sigma_lat;




   //////// likelihood

  //Calculate logLikelihood
  Type nll = 0.0;

  // behaviour equation
  Type behav = behaviour(dat.b, delta, A);
  REPORT(behav);
  nll -= behav;

  // movement equation
  vector<Type> movevec = movement(x, dat.b, gamma, theta, Sigma);
  nll += movevec.sum();
  REPORT(movevec);

  // measurement equation
  vector<Type> measvec = measurement(dat.y, x, dat.idx, dat.jidx, dat.ae, psi);
  REPORT(measvec);
  nll += measvec.sum();

  REPORT(nll);
  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
