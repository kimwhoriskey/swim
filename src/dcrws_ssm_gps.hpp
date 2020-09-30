#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace density;


template<class Type>
Type dcrwSSMgps(objective_function<Type> * obj) {

  // data

  // names of each track
  vector<std::string> tracknames = StringList("tracknames",obj);
  int ntracks = tracknames.size();
  vector<track<Type> > tracks(ntracks); // vector< kind of element> name(size)

  // tracks
  for(int i = 0; i < ntracks; i++){
    track<Type> indtrack(tracknames[i],obj);
    tracks[i] = indtrack;
  }



  REPORT(ntracks);


  // parameters

  // Input parameters - i.e. parameters to estimate.
  PARAMETER_VECTOR(working_theta);
  vector<Type> theta(2);
  theta(0) = 2*M_PI/(1.0+exp(-working_theta(0))) - M_PI; // -pi to pi
  theta(1) = 2*M_PI/(1.0+exp(-working_theta(1))); // 0 to 2pi
  ADREPORT(theta);

  PARAMETER_VECTOR(working_gamma);
  vector<Type> gamma(2);
  gamma(0) = 1.0/(1.0+exp(-working_gamma(0)));
  gamma(1) = gamma(0)/(1.0+exp(-working_gamma(1)));
  ADREPORT(gamma);

  PARAMETER(working_sigma_lon);
  Type sigma_lon = exp(working_sigma_lon);
  ADREPORT(sigma_lon);
  PARAMETER(working_sigma_lat);
  Type sigma_lat = exp(working_sigma_lat);
  ADREPORT(sigma_lat);
  // Variance-covariance matrix for the process equation.
  matrix<Type> Sigma(2,2);
  Sigma << sigma_lon*sigma_lon, 0.0,
  0.0, sigma_lat*sigma_lat;

  PARAMETER_MATRIX(working_A); //matrix of switching probabilities
  int m = working_A.rows(); // number of states
  matrix<Type> A(m, m);
  for(int i = 0; i < m; ++i){
      for(int j = 0; j < (m-1); ++j){
          A(i,j) = exp(working_A(i,j));  //ensures the probabilities we estimate are >0
      }
      A(i,m-1) = 1.0;
      A.row(i) *= 1.0/A.row(i).sum();
  }
  ADREPORT(A);
  // system of equations stat dist
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
  // vector of ones for stat dist
  vector<Type> ones(m);
  for(int i=0; i<m; i++) {
      ones(i) = 1.0;
  }
  // take the inverse of the system matrix and then multiply by the vector of ones
  matrix<Type> sys_inverse = system.inverse();
  vector<Type> delta = sys_inverse*ones;


  PARAMETER_VECTOR(working_gps_err); // working value for gps measurement error
  vector<Type> gps_err = exp(working_gps_err);
  matrix<Type> gpsSigma(2,2);
  gpsSigma << gps_err(0)*gps_err(0), 0.0,
  0.0, gps_err(1)*gps_err(1);





   // likelihood calculations

  //Calculate logLikelihood
  Type nll = 0.0;

  // behaviour equation
  vector<Type> behav(ntracks);
  Type tmpb = 0.0;
  for(int i=0; i < ntracks; ++i){
    tmpb = behaviour(tracks[i].b, delta, A);
    behav(i) = tmpb;
  }
  REPORT(behav);
  nll -= behav.sum();

  // movement equation
  Type tmpp = 0.0;
  vector<Type> move(ntracks);
  for(int i=0; i < ntracks; ++i){
    tmpp = movement(tracks[i].x, tracks[i].b, gamma, theta, Sigma);
    move(i) = tmpp;
  }
  REPORT(move);
  nll += move.sum();

  // measurement equation
  Type tmpm = 0.0;
  vector<Type> meas(ntracks);
  for(int i=0; i < ntracks; ++i){
    tmpm = measurement_gps(tracks[i].y, tracks[i].x, tracks[i].idx, tracks[i].jidx, tracks[i].ae, gpsSigma);
    meas(i) = tmpm;
  }

  REPORT(meas);
  nll += meas.sum();

  for(int i=0; i<ntracks; ++i) tracks[i].Rep();

  REPORT(nll);
  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
