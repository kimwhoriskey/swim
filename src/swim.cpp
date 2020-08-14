#include <TMB.hpp>
#include <math.h>
#include "extramacros.hpp"
#include "trackstruct.hpp"
#include "likelihoods.hpp"
#include "dcrws_hmm.hpp"
#include "dcrws_ssm.hpp"


// map the data_integer to a character
enum model_choices {
  dcrwhmm = 0,
  dcrwssm = 1,
};

template<class Type>
Type objective_function<Type>::operator() ()
{

DATA_INTEGER(model);

switch(model){
  case dcrwhmm:
    return dcrwHMM(this);
    break;
  case dcrwssm:
    return dcrwSSM(this);
    break;
  }
}
