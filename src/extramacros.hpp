// extra macros needed to read in a list of lists of data and a list of matrices of parameters
// jonathan babyn is the man

/**
 * Read in a list of strings from the R TMB data list
 * @param name the name of the list of strings
 * @param obj the pointer to the TMB objective function
 */
template<class Type>
vector<std::string> StringList(const char *name,objective_function<Type> *obj){
  SEXP strings = getListElement(obj -> data,name);
  vector<std::string> strs(LENGTH(strings));
  for(int i = 0; i < LENGTH(strings);i++){
    SEXP v_element = VECTOR_ELT(strings,i);
    strs[i] = CHAR(STRING_ELT(v_element,0));
  }
  return strs;
}


/**
 * Get named parameter matrix from R TMB parameter list with given prefix
 * @param name name of parameter matrix to retrieve from list
 * @param obj Pointer to TMB objective function
 * @param prefix the prefix to attach to the name for finding from the R parameter list
 */
template<class Type>
inline matrix<Type> fetchParameterMatrix(std::string name,objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  tmbutils::matrix<Type> ret(obj -> fillShape(asMatrix<Type> ( obj -> getShape(nam, &Rf_isMatrix) ),nam));
  return ret;
}
