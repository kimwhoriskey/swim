// track struct for ssm
// jonathan babyn is the man

// declaring data and methods in struct
template<class Type>
struct track {
  track();
  track(std::string name,objective_function<Type> *obj);
  array<Type> y;
  vector<int> b;
  vector<int> kind;
  vector<int> datatype;
  vector<int> idx;
  vector<Type> jidx;
  matrix<Type> ae;
  matrix<Type> x;
  char* namm;
  objective_function<Type> *object;
  void Rep();
};

// constructor 1
template<class Type>
track<Type>::track(std::string name, objective_function<Type> *obj) { // objective function for parameter
  char* nam = new char[name.size()+1];
  namm = nam;
  object = obj;
  std::strcpy(nam,name.c_str());
  SEXP tracking_data = getListElement(obj -> data,nam);
  y = tmbutils::asArray<Type>(getListElement(tracking_data, "y"));
  b = asVector<int>(getListElement(tracking_data, "b"));
  kind = asVector<int>(getListElement(tracking_data, "kind"));
  datatype = asVector<int>(getListElement(tracking_data, "datatype")) ;
  idx = asVector<int>(getListElement(tracking_data, "idx"));
  jidx = asVector<Type>(getListElement(tracking_data, "jidx"));
  ae = asMatrix<Type>(getListElement(tracking_data, "ae"));
  x = fetchParameterMatrix(".x", obj, name);
}

// right now just have one, can add more
template<class Type>
void track<Type>::Rep(){
  const char *names[] = {"x"};
  if(isDouble<Type>::value &&
     object -> current_parallel_region <0)
    {
      SEXP _TMB_temporary_sexp_ = PROTECT(allocVector(VECSXP,1));
      SEXP temp = asSEXP(x);
      SET_VECTOR_ELT(_TMB_temporary_sexp_,0,temp);
      Rf_defineVar(Rf_install(namm),_TMB_temporary_sexp_, object -> report);

  //Name components
      SEXP nms = PROTECT(allocVector(STRSXP,1));
      for(int i = 0; i < 1; i++){
	SET_STRING_ELT(nms, i, mkChar(names[i]));
      }
      setAttrib(_TMB_temporary_sexp_, R_NamesSymbol, nms);
      UNPROTECT(2);
    }
}


 // empty constructor, some kind of hack needed
template<class Type>
track<Type>::track(){
}






// track struct for hmm

template<class Type>
struct justthelocs {
  justthelocs();
  justthelocs(std::string name, objective_function<Type> *obj);
  matrix<Type> x;
  char* namm;
  objective_function<Type> *object;
};

template<class Type>
justthelocs<Type>::justthelocs(std::string name, objective_function<Type> *obj) { // objective function for parameter
  char* nam = new char[name.size()+1];
  namm = nam;
  object = obj;
  std::strcpy(nam,name.c_str());
  SEXP tracking_data = getListElement(obj -> data,nam);
  x = asMatrix<Type>(getListElement(tracking_data, "x"));
}

// empty constructor
template<class Type>
justthelocs<Type>::justthelocs(){
}



//
// // hmm results struct
// template<class Type>
// struct hmmrslts {
//   vector<int> states;
//   matrix<Type> pseudo;
// };
//




// hmm results struct
template<class Type>
struct hmmrslts {
  hmmrslts();
  hmmrslts(std::string name, objective_function<Type> *obj);
  vector<int> states;
  matrix<Type> pseudo;
  char* namm;
  objective_function<Type> *object;
  void Rep();
};

template<class Type>
hmmrslts<Type>::hmmrslts(std::string name, objective_function<Type> *obj) { // objective function for parameter
  char* nam = new char[name.size()+1];
  namm = nam;
  object = obj;
  std::strcpy(nam,name.c_str());
}

template<class Type>
void hmmrslts<Type>::Rep(){
  const char *names[] = {"states","pseudo"};
  if(isDouble<Type>::value &&
     object -> current_parallel_region <0)
    {
      SEXP _TMB_temporary_sexp_ = PROTECT(allocVector(VECSXP,2));
      SEXP temp = asSEXP(states);
      SET_VECTOR_ELT(_TMB_temporary_sexp_,0,temp);
      temp = asSEXP(pseudo);
      SET_VECTOR_ELT(_TMB_temporary_sexp_,1,temp);
      Rf_defineVar(Rf_install(namm),_TMB_temporary_sexp_, object -> report);

  //Name components
      SEXP nms = PROTECT(allocVector(STRSXP,2));
      for(int i = 0; i < 2; i++){
	SET_STRING_ELT(nms, i, mkChar(names[i]));
      }
      setAttrib(_TMB_temporary_sexp_, R_NamesSymbol, nms);
      UNPROTECT(2);
    }
}

template<class Type>
hmmrslts<Type>::hmmrslts(){
}
