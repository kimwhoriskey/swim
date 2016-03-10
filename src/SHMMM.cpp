 #include <TMB.hpp>
 
using namespace density;


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Input data
  DATA_MATRIX(x); //track coordinates
  DATA_IVECTOR(b);  //set of possible states
  DATA_FACTOR(stateSpace);  //values in state-space
  DATA_VECTOR(initDist);  //specify initial probability
  
  // Input parameters - i.e. parameters to estimate.
  PARAMETER(logitTheta1); //first turning angle
  PARAMETER(logitTheta2); //second turning angle, add to the first and attempt to make it larger
  PARAMETER(logitGamma1); // First autocorrelation - logit because 0 < gamma < 1
  PARAMETER(logitGamma2); // First autocorrelation - logit because 0 < gamma < 1
  PARAMETER(logSdlat); // Process standard deviation in lat - log because sdlat > 0
  PARAMETER(logSdlon); // Process standard deviation in lon - log because sdlon > 0
  
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
  Type sdLat = exp(logSdlat);
  Type sdLon = exp(logSdlon);
  
  // Report the parameters and their standard errors in their model format
  REPORT(theta1);
  REPORT(theta2);
  REPORT(gamma1);
  REPORT(gamma2);
  REPORT(sdLat);
  REPORT(sdLon);
  ADREPORT(theta1);
  ADREPORT(theta2);
  ADREPORT(gamma1);
  ADREPORT(gamma2);
  ADREPORT(sdLat);
  ADREPORT(sdLon);
  //ADREPORT(rho);


  //setting up the matrix of switching probabilities
//   matrix<Type> A(logA.rows(),logA.cols()+1); /*adding one more column to A than was 
//   present in logA */
  
//   A(0,0) = (1.0/2.0)*(1.0/(1.0+exp(-logA(0,0)))) + (1.0/2.0);
//   A(0,1) = 1.0 - A(0,0);
//   A(1,0) = (1.0/2.0)*(1.0/(1.0+exp(-logA(1,0))));
//   A(1,1) = 1.0 - A(1,0);

 //setting up the matrix of switching probabilities
 matrix<Type> A(logA.rows(),logA.cols()+1); /*adding one more column to A than was 
 present in logA */

 for(int i = 0; i < A.rows(); ++i){
   for(int j = 0; j < logA.cols(); ++j){
     A(i,j) = exp(logA(i,j));  //ensures the probabilities we estimate are >0
   }
   A(i,A.cols()-1) = 1.0; /*Remember that c++ is base 0, so the first position in a vector is 
   actually called position 0. So this is setting the last column (which is not present in logA
   equal to 1. We should be able to switch this 1 to any value we want, because the next line
   ensures that all of the probabilities in each row are intrinsically linked, and will add to 1.*/
   A.row(i) *= 1.0/A.row(i).sum(); /*keeps the sum across rows =1, specifically each element is 
   divided by its row sum*/
   //Note x *= 2 means x = 2*x
 }

  REPORT(A);
  ADREPORT(A);


  // Variance-covariance matrix for the process equation.
  matrix<Type> covp(2,2);
  covp << sdLon*sdLon, 0.0,
    0.0, sdLat*sdLat;

  /*Setting up the likelihood variable for the process equation. nll_dens will return the 
  negative log-likelihood value for the multivariate normal distribution with covariance
  matrix covp. */
  MVNORM_t<Type> nll_dens(covp);


  //Forward the Viterbi algorithm to calculate logLikelihood'
  Type ll = 0.0; /*Going to calculate the log-likelihood, and then multiply by -1 at the 
  end to get the negative log-likelihood*/
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
      T << cos(theta(b(j))), -sin(theta(b(j))),
           sin(theta(b(j))), cos(theta(b(j)));
      TMP = (x.col(i)-x.col(i-1)) - T*(x.col(i-1) - x.col(i-2))*Type(gamma(b(j)));
      
      phi(j) *= exp(-nll_dens(TMP));
    }
    
    
    ll += log(phi.sum());  /*Summing up each of the probabilities and then adding the log of this
    to the log-likelihood. */
    phi *= 1.0/phi.sum();  /*resets the phi matrix to ensure the probabilities of being in each state
    sum to 1. Then the calculation kicks in again. */
  }



  //Viterbi Algorithm to estimate states
  matrix<Type> v(x.cols(),stateSpace.size()); //T x s 
  /* These are our delta-values, also the forward or partial probabilities. Throughout
  we deal with log-probabilities.*/
  matrix<int> ba(x.cols(),stateSpace.size()); //T x s, 
  /*This gives the back pointer, it saves the identity of the most likely state at time 
  t-1 that led to a particular state of interest at time t.*/
  vector<int> states(x.cols()); //n, 
  /*This is where we will save the estimated states, i.e. bhat. */

  for(int j = 0; j < stateSpace.size(); ++j){
    v(0,j) = log(initDist(j));
    /*The first deltas are equal to the initial probabilities of being in 
    each state because we have no difference data for the first time step. */
  } 
   
//   TMP = x.col(1)-x.col(0);   
//   for(int j = 0; j < stateSpace.size(); ++j){
//     v(1,j) = log(initDist(j))-nll_dens(TMP); 
//     /*The second deltas are equal to multiplying the first delta by the density or likelihood
//     value of the first differences. Since we are using logs for v, this corresponds to adding
//     the first v, which are just equal to the log-initial probabilities, to the log-density
//     of the first differences. Because nll_dens returns the negative log-likelihood, we have
//     to multiply by -1.*/
//   }
  
  matrix<Type> T2(2,2);
  vector<Type> tmp(2);

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
      /*Now we're adding in the information from the observations*/
      T2 << cos(theta(b(j))), -sin(theta(b(j))),
            sin(theta(b(j))), cos(theta(b(j)));
      TMP = (x.col(k+1)-x.col(k)) - T2*(x.col(k) - x.col(k-1))*Type(gamma(b(j)));
      
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
    /*Now we're adding in the information from the observations*/
    
    v(x.cols()-1,j) = tmp;
  }
  

  
  
//   
//   
//   //Now we calculate the deltas for the rest of the data set
//   for(int k = 1; k < (x.cols()-1); ++k){
//     
//     for(int j = 0; j < stateSpace.size(); ++j){
//       Type tmp = v(k-1,0) + log(A(j,0));  
//       
//       ba(k,j) = 0; //back pointer
//       
//       for(int i = 1; i < stateSpace.size(); ++i){
//         if(v(k-1,i) + log(A(j,i)) > tmp){
//           tmp = v(k-1,i) + log(A(j,i)); 
//           ba(k,j) = i;
//         }
//       }
//       /*Now we're adding in the information from the observations*/
//       T2 << cos(theta(b(j))), -sin(theta(b(j))),
//             sin(theta(b(j))), cos(theta(b(j)));
//       TMP = (x.col(k+1)-x.col(k)) - T2*(x.col(k) - x.col(k-1))*Type(gamma(b(j)));
//       
//       v(k,j) = v(k-1,ba(k,j))-nll_dens(TMP);
//     }
//   }
//   
//   
//   //the last step  
//   for(int j = 0; j < stateSpace.size(); ++j){
//     Type tmp = v(x.cols()-2,0) + log(A(j,0));  
//     
//     ba(x.cols()-1,j) = 0; //back pointer
//     
//     for(int i = 1; i < stateSpace.size(); ++i){
//       if(v(x.cols()-2,i) + log(A(j,i)) > tmp){
//         tmp = v(x.cols()-2,i) + log(A(j,i)); 
//         ba(x.cols()-1,j) = i;
//       }
//     }
//     /*Now we're adding in the information from the observations*/
//     
//     v(x.cols()-1,j) = v(x.cols()-2,ba(x.cols()-1,j));
//   }
//   
//   
//   
  
  
  REPORT(v); //In case we want to look at these to troubleshoot. 
  REPORT(ba);
  
  //Writing out the states vector
  //setting the last value of the states vector to zero
  states(x.cols()-1) = 0;
  for(int i = 1; i < stateSpace.size(); ++i){
    if(v(x.cols()-1,i)>v(x.cols()-1,i-1)){
      states(x.cols()-1) = i;
      
    }
  }
  for(int k = 2; k < x.cols()+1; ++k){
    states(x.cols()-k) = ba(x.cols()-k+1,states(x.cols()-k+1)); /*works backwards from the final
    state estimate to get the sequence of previous states */
  }
  
  REPORT(states);

 


  return -ll;

}

