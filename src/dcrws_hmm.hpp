#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace density;


template<class Type>
Type dcrwHMM(objective_function<Type> * obj) {
    // Input data
    DATA_MATRIX(x); //track coordinates

    // Input parameters - i.e. parameters to estimate.
    PARAMETER_VECTOR(working_theta);
    PARAMETER_VECTOR(working_gamma); // First autocorrelation - logit because 0 < gamma < 1

    PARAMETER(working_tau_lon); // Process standard deviation in lon - log because sdlon > 0
    PARAMETER(working_tau_lat); // Process standard deviation in lat - log because sdlat > 0

    PARAMETER_MATRIX(working_A); //matrix of switching probabilities

    // Transformation of the input parameters to model format
    vector<Type> theta(2);
    theta(0) = 2*M_PI/(1.0+exp(-working_theta(0))) - M_PI; // -pi to pi
    theta(1) = 2*M_PI/(1.0+exp(-working_theta(1))); // 0 to 2pi
    // vector<Type> theta = 2*M_PI/(1.0+exp(-logitTheta1)) - M_PI; //not sure if needed

    // vector<Type> gamma = working_gamma; //1.0/(1.0+exp(-working_gamma));
//    vector<Type> gamma = 1.0/(1.0+exp(-working_gamma));
    vector<Type> gamma(2);
    gamma(0) = 1.0/(1.0+exp(-working_gamma(0)));
    gamma(1) = gamma(0)/(1.0+exp(-working_gamma(1)));

    Type tau_lon = exp(working_tau_lon);
    Type tau_lat = exp(working_tau_lat);

    // Report the parameters and their standard errors in their model format
    REPORT(theta);
    REPORT(gamma);
    REPORT(tau_lon);
    REPORT(tau_lat);
    ADREPORT(theta);
    ADREPORT(gamma);
    ADREPORT(tau_lon);
    ADREPORT(tau_lat);


    //setting up the matrix of switching probabilities
    int m = working_A.rows(); // number of states
    int n_obs = x.cols();

    matrix<Type> A(m, m); /*adding one more column to A than was
                           present in logA */

    for(int i = 0; i < m; ++i){
        for(int j = 0; j < (m-1); ++j){
            A(i,j) = exp(working_A(i,j));  //ensures the probabilities we estimate are >0
        }
        A(i,m-1) = 1.0;
        A.row(i) *= 1.0/A.row(i).sum();
    }

    REPORT(A);
    ADREPORT(A);
    REPORT(m);


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

    matrix<Type> delta_row(1,m);
    for(int i=0; i<m; i++) {
        delta_row(0,i) = delta(i); //creating a new vector identical to start_probs
    }

    //
    REPORT(delta);
    REPORT(delta_row);
    //



    // Variance-covariance matrix for the process equation.
    matrix<Type> Tau(2,2);
    Tau << tau_lon*tau_lon, 0.0,
    0.0, tau_lat*tau_lat;

    //Setting up process error
    MVNORM_t<Type> nll_dens(Tau);


    //// probability density matrices
    matrix<Type> T(2,2); //rotational matrix
    vector<Type> TMP(2); //temporary variable

    array<Type> P_array(m, m, n_obs);
    array<Type> C_lon_array(m, m, n_obs);
    array<Type> C_lat_array(m, m, n_obs);
    for(int k=0; k<n_obs; k++) {    // k indexes time
        for(int i=0; i<m; i++) {
            for(int j=0; j<m; j++) { // i and j both index behavioural state
                if( i==j ) { // diagonal entries
                    //process equation
                    if(k == 0){
                        P_array(i,j,k) = 0.0;
                        C_lon_array(i,j,k) = 0.0;
                        C_lat_array(i,j,k) = 0.0;
                    } else if (k == 1){
                        TMP = x.col(1)-x.col(0);
                        P_array(i,j,k) = exp(-nll_dens(TMP));
                        C_lon_array(i,j,k) = pnorm(TMP(0)/tau_lon);
                        C_lat_array(i,j,k) = pnorm(TMP(1)/tau_lat);
                    } else {
                        T << cos(theta(i)), -sin(theta(i)),
                        sin(theta(i)), cos(theta(i));
                        TMP = (x.col(k)-x.col(k-1)) - T*(x.col(k-1) - x.col(k-2))*Type(gamma(i));
                        P_array(i,j,k) = exp(-nll_dens(TMP));
                        C_lon_array(i,j,k) = pnorm(TMP(0)/tau_lon);
                        C_lat_array(i,j,k) = pnorm(TMP(1)/tau_lat);
                    }
                } else {
                    P_array(i,j,k) = 0.0; //off-diagonals are zero
                    C_lon_array(i,j,k) = 0.0;
                    C_lat_array(i,j,k) = 0.0;
                }
            }
        }
    } // matrix form for ease of forward algorithm

    REPORT(P_array);

    array<Type> viterbi(m, n_obs);    // Columns index time, rows index state
    array<Type> forward_max(m);    // Vector of most likely states
    // array<Type> P_array = ta_array * sl_array; // the matrix of pdf's
    vector<int> states(n_obs);


    int min = 0;
    int max = n_obs;

    //starting state likelihoods
    for(int i=0; i<m; i++) {
        viterbi(i,min) = log(delta(i));
    }

    for(int k=(min+1); k<max; k++) {
        for(int i=0; i<m; i++) {
            for(int j=0; j<m; j++) {
                forward_max(j) = viterbi(j,k-1) + log(A(i,j)); // transition probabilities
            }
            viterbi(i,k) = forward_max.maxCoeff() + log(P_array(i,i,k));    // Choose most likely transition
        }
    }

    REPORT(viterbi);

    Eigen::Index max_row; // looks like this is the index of a vector
    Eigen::Index max_col;
    Type foo;

    // Get the last column of viterbi matrix for backwards probabilities
    matrix<Type> column = viterbi.col(max-1);  // column matrix, the last column, should be mx1
    foo = column.maxCoeff(&max_row, &max_col);    // Sends index of largest coefficient to max_row (don't care about max_col)
    //okay so foo contains the maximum value. max_row contains the row where foo occurs.
    //max_col contains the column at which foo occurs (doesn't matter cause just one column)
    states(max-1) = max_row + 1;    // Eigen::Index starts at 0, want it to start at 1

    // Backwards recursion
    for(int k=n_obs-2; k>=0; k--) {
        for(int i=0; i<m; i++) { // index state
            column(i,0) = viterbi(i,k) + log(A(i,max_row));    // backwards Viterbi coefficient
        }
        foo = column.maxCoeff(&max_row, &max_col);    // sends index of largest coefficient to max_row
        // figures out which of the previous states led to the next one
        states(k) = max_row + 1;    // Eigen::Index starts at 0, want it to start at 1
    }

    REPORT(states);
    //



    Type ll = 0.0;

    // Forward Probabilities as row vector
    matrix<Type> alpha(1,m);
    // let's store them all , then do the rest of the pseudoresids in R
    matrix<Type> alphas(n_obs, m);
    matrix<Type> pseudo(n_obs, 2);

    for(int k=0; k<n_obs; k++) {
        if( k == 0 ) {
            pseudo(k) = delta.sum();
            alpha = delta_row;
            alphas.row(k) = alpha;
        } else {
            pseudo(k,0) = ( alpha*A * ( C_lon_array.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth longitude
            pseudo(k,1) = ( alpha*A * ( C_lat_array.col(k).matrix())).sum() / alpha.sum();   // pseudoresidual for the kth latitude
            alpha = alpha* ( A ) * ( P_array.col(k).matrix() ); // Add k'th observation to forward probability
            alphas.row(k) = alpha;
        }
        ll += log(alpha.sum());  // add log of forward probability to recursive likelihood
        alpha = alpha/alpha.sum();    // rescale forward probabilities to sum to 1
    }

    REPORT(alphas);
    REPORT(pseudo);
    return -ll; //+ dummy*dummy;




}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
