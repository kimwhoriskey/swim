#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace density;


template<class Type>
Type dcrwHMM(objective_function<Type> * obj) {

    // data

    vector<std::string> tracknames = StringList("tracknames",obj);
    int ntracks = tracknames.size();

    vector<justthelocs<Type> > tracks(ntracks); // vector< kind of element> name(size)
    for(int i = 0; i < ntracks; i++){
      justthelocs<Type> indtrack(tracknames[i],obj);
      tracks[i] = indtrack;
    }



    // parameters

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

    PARAMETER(working_tau_lon);
    Type tau_lon = exp(working_tau_lon);
    ADREPORT(tau_lon);
    PARAMETER(working_tau_lat);
    Type tau_lat = exp(working_tau_lat);
    ADREPORT(tau_lat);
    // Variance-covariance matrix for the process equation.
    matrix<Type> Tau(2,2);
    Tau << tau_lon*tau_lon, 0.0,
    0.0, tau_lat*tau_lat;

    PARAMETER_MATRIX(working_A); //matrix of switching probabilities
    int m = working_A.rows(); // number of states
    // int n_obs = x.cols();
    matrix<Type> A(m, m);
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < (m-1); ++j){
            A(i,j) = exp(working_A(i,j));
        }
        A(i,m-1) = 1.0;
        A.row(i) *= 1.0/A.row(i).sum();
    }
    ADREPORT(A);


    // system of eqns for stat dist
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
    REPORT(delta);





    // // hmm results struct
    // template<class Type>
    // struct hmmrslts {
    //   vector<int> states;
    //   matrix<Type> pseudo;
    // };

    vector<hmmrslts<Type> > rslts(ntracks);
    vector<Type> lls(ntracks);

    // for(int i=0; i < ntracks; ++i){
    //
    //   // get pdfs
    //   array<Type> pdfs(m, m, tracks[i].x.cols());
    //   pdfs = fdens(m, tracks[i].x, theta, gamma, Tau, false, false);
    //
    //   // run viterbi algorithm
    //   // vector<int> states(n_obs);
    //   rslts[i].states = viterbi(m, tracks[i].x, delta, A, pdfs);
    //
    //   // get cdfs of lon and lat
    //   array<Type> cdfslon(m, m, tracks[i].x.cols());
    //   cdfslon = fdens(m, tracks[i].x, theta, gamma, Tau, true, true);
    //   array<Type> cdfslat(m, m, tracks[i].x.cols());
    //   cdfslat = fdens(m, tracks[i].x, theta, gamma, Tau, true, false);
    //
    //   // run the forward algorithm, get pseudos and ll
    //   forwardalgrslt<Type> forward;
    //   forward  = forwardalg(tracks[i].x, delta, A, pdfs, cdfslon, cdfslat);
    //   rslts[i].pseudo = forward.pseudo;
    //   lls[i] = forward.ll;
    //
    //
    // }


    for(int i=0; i < ntracks; ++i){

      hmmrslts<Type> singlerslt(tracknames[i],obj);

      // get pdfs
      array<Type> pdfs(m, m, tracks[i].x.cols());
      pdfs = fdens(m, tracks[i].x, theta, gamma, Tau, false, false);

      // run viterbi algorithm
      // vector<int> states(n_obs);
      singlerslt.states = viterbi(m, tracks[i].x, delta, A, pdfs);

      // get cdfs of lon and lat
      array<Type> cdfslon(m, m, tracks[i].x.cols());
      cdfslon = fdens(m, tracks[i].x, theta, gamma, Tau, true, true);
      array<Type> cdfslat(m, m, tracks[i].x.cols());
      cdfslat = fdens(m, tracks[i].x, theta, gamma, Tau, true, false);

      // run the forward algorithm, get pseudos and ll
      forwardalgrslt<Type> forward;
      forward  = forwardalg(tracks[i].x, delta, A, pdfs, cdfslon, cdfslat);
      singlerslt.pseudo = forward.pseudo;
      lls[i] = forward.ll;

      rslts[i] = singlerslt; 
      rslts[i].Rep();

    }

    REPORT(ntracks);




    return -lls.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
