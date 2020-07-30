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

    REPORT(delta);




    // Variance-covariance matrix for the process equation.
    matrix<Type> Tau(2,2);
    Tau << tau_lon*tau_lon, 0.0,
    0.0, tau_lat*tau_lat;



    //// probability density matrices
    matrix<Type> T(2,2); //rotational matrix
    vector<Type> TMP(2); //temporary variable


    array<Type> pdfs(m, m, n_obs);
    pdfs = fdens(m, x, theta, gamma, Tau, false, false);
    REPORT(pdfs);


    // viterbi algorithm
    vector<int> states(n_obs);
    states = viterbi(m, x, delta, A, pdfs);
    REPORT(states);


    array<Type> cdfslon(m, m, n_obs);
    cdfslon = fdens(m, x, theta, gamma, Tau, true, true);
    REPORT(cdfslon);

    array<Type> cdfslat(m, m, n_obs);
    cdfslat = fdens(m, x, theta, gamma, Tau, true, false);
    REPORT(cdfslat);


    forwardalgrslt<Type> forward;
    forward  = forwardalg(x, delta, A, pdfs, cdfslon, cdfslat);
    matrix<Type> pseudo = forward.pseudo;
    Type ll = forward.ll;
    REPORT(pseudo);
    REPORT(ll);

    return -ll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
