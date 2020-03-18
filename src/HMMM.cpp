 #include <TMB.hpp>
 
using namespace density;


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Input data
  DATA_MATRIX(x); //track coordinates
  DATA_FACTOR(stateSpace);  //values in state-space
  DATA_VECTOR(initDist);  //specify initial probability
  
  // Input parameters - i.e. parameters to estimate.
  PARAMETER(logitTheta1); //first turning angle
  PARAMETER(logitTheta2); //second turning angle, add to the first and attempt to make it larger
  PARAMETER(logitGamma1); // First autocorrelation - logit because 0 < gamma < 1
  PARAMETER(logitGamma2); // First autocorrelation - logit because 0 < gamma < 1
  PARAMETER(logSdlon); // Process standard deviation in lon - log because sdlon > 0
  PARAMETER(logSdlat); // Process standard deviation in lat - log because sdlat > 0

  PARAMETER_MATRIX(logA); //matrix of switching probabilities
  
  // Transformation of the input parameters to model format
  Type theta1 = 2*M_PI/(1.0+exp(-logitTheta1)) - M_PI; //between -pi and pi
  Type theta2 = 2*M_PI/(1.0+exp(-logitTheta2)); //between 0 and 2*pi
  vector<Type> theta(2);
  theta(0) = theta1;
  theta(1) = theta2;
  Type gamma1 = 1.0/(1.0+exp(-logitGamma1));
  //Type gamma2 = gamma1/(1.0+exp(-logitGamma2)); //autocorrelation gamma2 < gamma1
  Type gamma2 = 1.0/(1.0+exp(-logitGamma2));  //gamma's are unrelated
  vector<Type> gamma(2);
  gamma(0) = gamma1;
  gamma(1) = gamma2;
  Type sdLon = exp(logSdlon);
  Type sdLat = exp(logSdlat);

  // Report the parameters and their standard errors in their model format
  REPORT(theta1);
  REPORT(theta2);
  REPORT(gamma1);
  REPORT(gamma2);
  REPORT(sdLon);
  REPORT(sdLat);
  ADREPORT(theta1);
  ADREPORT(theta2);
  ADREPORT(gamma1);
  ADREPORT(gamma2);
  ADREPORT(sdLon);
  ADREPORT(sdLat);


 //setting up the matrix of switching probabilities
 matrix<Type> A(logA.rows(),logA.cols()+1); /*adding one more column to A than was 
 present in logA */

 for(int i = 0; i < A.rows(); ++i){
   for(int j = 0; j < logA.cols(); ++j){
     A(i,j) = exp(logA(i,j));  //ensures the probabilities we estimate are >0
   }
   A(i,A.cols()-1) = 1.0; 
   A.row(i) *= 1.0/A.row(i).sum(); 
 }

  REPORT(A);
  ADREPORT(A);


  // Variance-covariance matrix for the process equation.
  matrix<Type> Sigma(2,2);
  Sigma << sdLon*sdLon, 0.0,
    0.0, sdLat*sdLat;

  //Setting up process error
  MVNORM_t<Type> nll_dens(Sigma);


  // Calculate log-likelihood
  Type ll = 0.0; 
  matrix<Type> phi(1,stateSpace.size()); /* Here we are setting up the matrix of partial 
  or forward probabilities. The size is the number of states*/
  matrix<Type> T(2,2); //rotational matrix
  vector<Type> TMP(2); //temporary variable
  
  //Initial contribution
  for(int i = 0; i < stateSpace.size(); ++i){
    phi(i) = initDist(i);  /*this is data input; an initial state probability vector.  */
  }

  TMP = x.col(1)-x.col(0);
  for(int j=0; j < stateSpace.size(); ++j){
    phi(j) *= exp(-nll_dens(TMP));
  }
  ll += log(phi.sum());  //add contribution to the likelihood
  phi *= 1.0/phi.sum();  //reset the phi matrix

  //the rest of the data
  for(int i = 2; i < x.cols(); ++i){
      phi = phi*A;  //multiply phi by the switching probabilities
      
    for(int j = 0; j < stateSpace.size(); ++j){
      
      //process equation
      T << cos(theta(j)), -sin(theta(j)),
           sin(theta(j)), cos(theta(j));
      TMP = (x.col(i)-x.col(i-1)) - T*(x.col(i-1) - x.col(i-2))*Type(gamma(j));
      
      phi(j) *= exp(-nll_dens(TMP));
    }
    
    ll += log(phi.sum());  /*Summing up each of the probabilities and then adding the log of this
    to the log-likelihood. */
    phi *= 1.0/phi.sum();  /*resets the phi matrix to ensure the probabilities of being in each state
    sum to 1. Then the calculation kicks in again. */
  }



  //Viterbi Algorithm to estimate states
  matrix<Type> v(x.cols(),stateSpace.size()); //T x s 
  /* These are our delta-values, also the partial probabilities. Throughout
  we deal with log-probabilities.*/
  matrix<int> ba(x.cols(),stateSpace.size()); //T x s, 
  /*This gives the back pointer, it saves the identity of the most likely state at time 
  t-1 that led to a particular state of interest at time t.*/
  vector<int> states(x.cols()); //n, 
  /*This is where we will save the estimated states, i.e. bhat. */


  
  matrix<Type> T2(2,2);
  vector<Type> tmp(2);
  
  tmp = x.col(1)-x.col(0);   
  for(int j = 0; j < stateSpace.size(); ++j){
    v(0,j) = log(initDist(j))-nll_dens(tmp); 
  }
  
  //Now we calculate the deltas for the rest of the data set
  for(int k = 1; k < (x.cols()-1); ++k){
    for(int j = 0; j < stateSpace.size(); ++j){
      Type tmp = v(k-1,0) + log(A(0,j));  
      ba(k,j) = 0; //back pointer
      for(int i = 1; i < stateSpace.size(); ++i){
        if(v(k-1,i) + log(A(i,j)) > tmp){
          tmp = v(k-1,i) + log(A(i,j)); 
          ba(k,j) = i;
        }
      }
      T2 << cos(theta(j)), -sin(theta(j)),
            sin(theta(j)), cos(theta(j));
      TMP = (x.col(k+1)-x.col(k)) - T2*(x.col(k) - x.col(k-1))*Type(gamma(j));
      v(k,j) = tmp-nll_dens(TMP);
    }
  }
  
  
//the last step  
  for(int j = 0; j < stateSpace.size(); ++j){
    Type tmp = v(x.cols()-2,0) + log(A(0,j));  
    
    ba(x.cols()-1,j) = 0; //back pointer
    
    for(int i = 1; i < stateSpace.size(); ++i){
      if(v(x.cols()-2,i) + log(A(i,j)) > tmp){
        tmp = v(x.cols()-2,i) + log(A(i,j)); 
        ba(x.cols()-1,j) = i;
      }
    }

    v(x.cols()-1,j) = tmp; 
  }

  
  REPORT(v); //In case we want to look at these to troubleshoot. 
  REPORT(ba);
  
  //Decoding the states vector
  states(x.cols()-1) = 0;
  for(int i = 1; i < stateSpace.size(); ++i){
    if(v(x.cols()-1,i)>v(x.cols()-1,i-1)){
      states(x.cols()-1) = i;
    }
  }
  for(int k = 2; k < x.cols()+1; ++k){
    states(x.cols()-k) = ba(x.cols()-k+1,states(x.cols()-k+1)); 
  }
  
  REPORT(states);

 


  return -ll;

}

