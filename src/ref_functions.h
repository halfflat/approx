#pragma once

#include <cmath>
#include <functional>

#include "unary_fn.h"

// Wrappers, metadata for reference functions used for
// evaluating approximations.
//
// Will wrap stdlib and use libarb if available.

extern unary_fn<double> std_exp_double;
extern unary_fn<double> std_expm1_double;
extern unary_fn<double> std_log_double;
extern unary_fn<double> std_log1p_double;

#ifdef USE_ARB

extern unary_fn<double> arb_exp_double;
extern unary_fn<double> arb_expm1_double;
extern unary_fn<double> arb_log_double;
extern unary_fn<double> arb_log1p_double;
extern unary_fn<double> arb_exprel_double;
extern unary_fn<double> arb_exprelr_double;

#endif
