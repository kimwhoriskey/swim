#include <TMB.hpp>
#include <math.h>
#include "extramacros.hpp"
#include "trackstruct.hpp"
#include "likelihoods.hpp"
#include "measlikelihood.hpp"
#include "dcrws_hmm.hpp"
#include "dcrws_ssm_combo.hpp"
#include "dcrws_ssm_argos.hpp"
#include "dcrws_ssm_gps.hpp"


// map the data_integer to a character
enum model_choices {
  dcrwhmm = 0,
  dcrwssmargos = 1,
  dcrwssmgps = 2,
  dcrwssmcombo = 3,
};

template<class Type>
Type objective_function<Type>::operator() ()
{

DATA_INTEGER(model);

switch(model){
  case dcrwhmm:
    return dcrwHMM(this);
    break;
  case dcrwssmargos:
    return dcrwSSMargos(this);
    break;
  case dcrwssmgps:
    return dcrwSSMgps(this);
    break;
  case dcrwssmcombo:
    return dcrwSSMcombo(this);
    break;
  }
}
