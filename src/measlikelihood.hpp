
// measurement process
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type measurement(objective_function<Type> * obj, array<Type> y, matrix<Type> x, vector<int> kind, vector<int> idx, vector<Type> jidx, matrix<Type> ae, int datatype){

  vector<Type> measnll(y.cols());
  vector<Type> xhat(2);
  vector<Type> tmp(2);
  vector<Type> tmp2(2);
  // density::MVNORM_t<Type> gps_nll(gpsSigma);

  if(datatype==0){

    PARAMETER(working_psi); //working value for measurement error
    Type psi = exp(working_psi);
    ADREPORT(psi);

    for(int i = 0; i < y.cols(); ++i){
        xhat = (1.0-jidx(i)) * x.col(idx(i)-1) + jidx(i)*x.col(idx(i));
        tmp = y.col(i);
        tmp2 = tmp-xhat;

        // taken from Marie's DCRW_Argos.cpp file
        // dt returns the log likelihood
        measnll(i) = -(0.5*log(psi) - log(ae(i,0)) + dt(sqrt(psi)*tmp2(0)/ae(i,0),ae(i,1),true) ); // Longitude
        measnll(i) -= (0.5*log(psi) - log(ae(i,2)) + dt(sqrt(psi)*tmp2(1)/ae(i,2),ae(i,3),true) ); // Latitude

    }

  } else if (datatype==1){

    PARAMETER_VECTOR(working_gps_err); // working value for gps measurement error
    vector<Type> gps_err = exp(working_gps_err);
    ADREPORT(gps_err);
    matrix<Type> gpsSigma(2,2);
    gpsSigma << gps_err(0)*gps_err(0), 0.0,
    0.0, gps_err(1)*gps_err(1);
    density::MVNORM_t<Type> gps_nll(gpsSigma);

    for(int i = 0; i < y.cols(); ++i){

        xhat = (1.0-jidx(i)) * x.col(idx(i)-1) + jidx(i)*x.col(idx(i));
        tmp = y.col(i);
        tmp2 = tmp-xhat;

        measnll(i) = gps_nll(tmp2); //MVNORM_t returns the negative log likelihood

    }


  } else if (datatype==2) {

    PARAMETER(working_psi); //working value for measurement error
    Type psi = exp(working_psi);
    ADREPORT(psi);
    PARAMETER_VECTOR(working_gps_err); // working value for gps measurement error
    vector<Type> gps_err = exp(working_gps_err);
    ADREPORT(gps_err);
    matrix<Type> gpsSigma(2,2);
    gpsSigma << gps_err(0)*gps_err(0), 0.0,
    0.0, gps_err(1)*gps_err(1);
    density::MVNORM_t<Type> gps_nll(gpsSigma);

    for(int i = 0; i < y.cols(); ++i){
        xhat = (1.0-jidx(i)) * x.col(idx(i)-1) + jidx(i)*x.col(idx(i));
        tmp = y.col(i);
        tmp2 = tmp-xhat;

        // different measurement error dists
        // 0 for gps, 1 for argos LS, could do 2 for Argos ellipse
        if(kind(i) == 0){
          measnll(i) = gps_nll(tmp2); //MVNORM_t returns the negative log likelihood
        } else if(kind(i) == 1) {
          // taken from Marie's DCRW_Argos.cpp file
          // dt returns the log likelihood
          measnll(i) = -(0.5*log(psi) - log(ae(i,0)) + dt(sqrt(psi)*tmp2(0)/ae(i,0),ae(i,1),true) ); // Longitude
          measnll(i) -= (0.5*log(psi) - log(ae(i,2)) + dt(sqrt(psi)*tmp2(1)/ae(i,2),ae(i,3),true) ); // Latitude
        }

    }
  }




  Type rslt = measnll.sum();
  return rslt; // negative log-lik
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
