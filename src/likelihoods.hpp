//behaviour process
template<class Type>
Type behaviour(vector<int> b, vector<Type> delta, matrix<Type> A){

  Type res = log(delta(b(0)));

  for(int i=1; i < b.size(); ++i){
      res += log(A(b(i-1), b(i)));
      // the A matrix is the previous state in rows, the next state in columns
      // so the likelihood of this observed markov chain is just the right entries of the matrix added together
  }

  return res; // log likelihood
}



//movement process
template<class Type>
Type movement(matrix<Type> x, vector<int> b, vector<Type> gamma, vector<Type> theta, matrix<Type> Sigma){

  density::MVNORM_t<Type> nllproc(Sigma);
  vector<Type> tmp(2);
  vector<Type> movenll(x.cols());
  matrix<Type> T(2,2);


  //first loc
  tmp = x.col(1)-x.col(0);
  movenll(1) = nllproc(tmp);

  // other locs
  for(int i = 2; i < x.cols(); ++i){

      T << cos(Type(theta(b(i-1)))), -sin(Type(theta(b(i-1)))),
      sin(Type(theta(b(i-1)))), cos(Type(theta(b(i-1))));
      tmp = (x.col(i)-x.col(i-1)) - T*(x.col(i-1) - x.col(i-2))*Type(gamma(b(i)));
      movenll(i) = nllproc(tmp);

  }

  Type rslt = movenll.sum();
  return rslt;
}

// measurement process
template<class Type>
Type measurement(array<Type> y, matrix<Type> x, vector<int> idx, vector<Type> jidx, matrix<Type> ae, Type psi){

  vector<Type> measnll(y.cols());
  vector<Type> xhat(2);
  vector<Type> tmp(2);
  vector<Type> tmp2(2);

  for(int i = 0; i < y.cols(); ++i){
      xhat = (1.0-jidx(i)) * x.col(idx(i)-1) + jidx(i)*x.col(idx(i));
      tmp = y.col(i);
      tmp2 = tmp-xhat;
      // nll += nll_meas(tmp-xhat); //MVNORM_t returns the negative log likelihood
      // taken from Marie's DCRW_Argos.cpp file
      // dt returns the log likelihood
      measnll(i) = -(0.5*log(psi) - log(ae(i,0)) + dt(sqrt(psi)*tmp2(0)/ae(i,0),ae(i,1),true) ); // Longitude
      measnll(i) -= (0.5*log(psi) - log(ae(i,2)) + dt(sqrt(psi)*tmp2(1)/ae(i,2),ae(i,3),true) ); // Latitude
  }

  Type rslt = measnll.sum();
  return rslt; // negative log-lik
}






//hmm likelihood calc






// probability density arrays
template<class Type>
array<Type> fdens(int m, matrix<Type> x, vector<Type> theta, vector<Type> gamma, matrix<Type> Tau, bool cdens, bool lon){

  int n_obs = x.cols();
  array<Type> densarray(m, m, n_obs);
  matrix<Type> T(2,2); //rotational matrix
  vector<Type> tmp(2); //temporary variable
  density::MVNORM_t<Type> nll_dens(Tau);

  for(int k=0; k<n_obs; k++) {    // k indexes time
      for(int i=0; i<m; i++) {
          for(int j=0; j<m; j++) { // i and j both index behavioural state
              if( i==j ) { // diagonal entries
                  //process equation
                  if(k == 0){
                      densarray(i,j,k) = 0.0;
                  } else if (k == 1){
                      tmp = x.col(1)-x.col(0);
                      if(cdens){
                        if(lon){
                          densarray(i,j,k) = pnorm(tmp(0)/Type(sqrt(Tau(0,0))));
                        } else {
                          densarray(i,j,k) = pnorm(tmp(1)/Type(sqrt(Tau(1,1))));
                        }
                      } else {
                        densarray(i,j,k) = exp(-nll_dens(tmp));
                      }
                  } else {
                      T << cos(theta(i)), -sin(theta(i)),
                      sin(theta(i)), cos(theta(i));
                      tmp = (x.col(k)-x.col(k-1)) - T*(x.col(k-1) - x.col(k-2))*Type(gamma(i));
                      if(cdens){
                        if(lon){
                          densarray(i,j,k) = pnorm(tmp(0)/Type(sqrt(Tau(0,0))));
                        } else {
                          densarray(i,j,k) = pnorm(tmp(1)/Type(sqrt(Tau(1,1))));
                        }
                      } else {
                        densarray(i,j,k) = exp(-nll_dens(tmp));
                      }
                  }
              } else {
                  densarray(i,j,k) = 0.0; //off-diagonals are zero
              }
          }
      }
  } // matrix form for ease of forward algorithm

  return densarray;

}



// viterbi algorithm to get the states

template<class Type>
vector<int> viterbi(int m, matrix<Type> x, vector<Type> delta, matrix<Type> A, array<Type> P_array){

  array<Type> forwards(m, x.cols());    // Columns index time, rows index state
  array<Type> forward_max(m);    // Vector of most likely states
  vector<int> states(x.cols());


  //starting state likelihoods
  for(int i=0; i<m; i++) {
      forwards(i,0) = log(delta(i));
  }

  for(int k=1; k<x.cols(); k++) {
      for(int i=0; i<m; i++) {
          for(int j=0; j<m; j++) {
              forward_max(j) = forwards(j,k-1) + log(A(i,j)); // transition probabilities
          }
          forwards(i,k) = forward_max.maxCoeff() + log(P_array(i,i,k));    // Choose most likely transition
      }
  }

  Eigen::Index max_row; // looks like this is the index of a vector
  Eigen::Index max_col;
  Type foo;

  // Get the last column of viterbi matrix for backwards probabilities
  matrix<Type> column = forwards.col(x.cols()-1);  // column matrix, the last column, should be mx1
  foo = column.maxCoeff(&max_row, &max_col);    // Sends index of largest coefficient to max_row (don't care about max_col)
  //okay so foo contains the maximum value. max_row contains the row where foo occurs.
  //max_col contains the column at which foo occurs (doesn't matter cause just one column)
  states(x.cols()-1) = max_row + 1;    // Eigen::Index starts at 0, want it to start at 1

  // Backwards recursion
  for(int k=x.cols()-2; k>=0; k--) {
      for(int i=0; i<m; i++) { // index state
          column(i,0) = forwards(i,k) + log(A(i,max_row));    // backwards Viterbi coefficient
      }
      foo = column.maxCoeff(&max_row, &max_col);    // sends index of largest coefficient to max_row
      // figures out which of the previous states led to the next one
      states(k) = max_row + 1;
  }

  return states;

}



// likelihood calculation



// just the pseudoresiduals
template<class Type>
matrix<Type> pseudoresids(matrix<Type> x, vector<Type> delta, matrix<Type> A, array<Type> P_array, array<Type> cdfslon, array<Type> cdfslat){

  matrix<Type> pseudo(x.cols(), 2);
  matrix<Type> alphas(x.cols(), A.rows());
  matrix<Type> alpha(1,A.rows());
  matrix<Type> delta_row(1,A.rows());
    for(int i=0; i<A.rows(); i++) {
        delta_row(0,i) = delta(i); //creating a new vector identical to start_probs
    }


  for(int k=0; k<x.cols(); k++) {
      if( k == 0 ) {
          pseudo.row(k) = delta_row;
          alpha = delta_row;
          alphas.row(k) = alpha;
      } else {
          pseudo(k,0) = ( alpha*A * ( cdfslon.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth longitude
          pseudo(k,1) = ( alpha*A * ( cdfslat.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth latitude
          alpha = alpha* ( A ) * ( P_array.col(k).matrix() ); // Add k'th observation to forward probability
          alphas.row(k) = alpha;
      }
      alpha = alpha/alpha.sum();    // rescale forward probabilities to sum to 1
  }
  return pseudo;
}

// just the likelihood
template<class Type>
Type loglik(matrix<Type> x, vector<Type> delta, matrix<Type> A, array<Type> P_array){

  matrix<Type> alphas(x.cols(), A.rows());
  matrix<Type> alpha(1,A.rows());
  matrix<Type> delta_row(1,A.rows());
    for(int i=0; i<A.rows(); i++) {
        delta_row(0,i) = delta(i); //creating a new vector identical to start_probs
    }
  Type ll = 0.0;

  for(int k=0; k<x.cols(); k++) {
      if( k == 0 ) {
          alpha = delta_row;
          alphas.row(k) = alpha;
      } else {
          alpha = alpha* ( A ) * ( P_array.col(k).matrix() ); // Add k'th observation to forward probability
          alphas.row(k) = alpha;
      }
      ll += log(alpha.sum());  // add log of forward probability to recursive likelihood
      alpha = alpha/alpha.sum();    // rescale forward probabilities to sum to 1
  }
  return ll;
}


// forward algorithm
template<class Type>
struct forwardalgrslt{
  matrix<Type> pseudo;
  Type ll;
};
// forward alg
template<class Type>
forwardalgrslt<Type> forwardalg(matrix<Type> x, vector<Type> delta, matrix<Type> A, array<Type> P_array, array<Type> cdfslon, array<Type> cdfslat){

  forwardalgrslt<Type> rslt;
  matrix<Type> _pseudo(x.cols(), 2);
  matrix<Type> alphas(x.cols(), A.rows());
  matrix<Type> alpha(1,A.rows());
  matrix<Type> delta_row(1,A.rows());
    for(int i=0; i<A.rows(); i++) {
        delta_row(0,i) = delta(i); //creating a new vector identical to start_probs
    }
  Type _ll = 0.0;


  for(int k=0; k<x.cols(); k++) {
      if( k == 0 ) {
          _pseudo.row(k) = delta_row;
          alpha = delta_row;
          alphas.row(k) = alpha;
      } else {
          _pseudo(k,0) = ( alpha*A * ( cdfslon.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth longitude
          _pseudo(k,1) = ( alpha*A * ( cdfslat.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth latitude
          alpha = alpha* ( A ) * ( P_array.col(k).matrix() ); // Add k'th observation to forward probability
          alphas.row(k) = alpha;
      }
      _ll += log(alpha.sum());  // add log of forward probability to recursive likelihood
      alpha = alpha/alpha.sum();    // rescale forward probabilities to sum to 1
  }

  rslt.pseudo = _pseudo;
  rslt.ll = _ll;
  return rslt;
}
