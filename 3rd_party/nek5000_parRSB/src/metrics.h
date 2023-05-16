#ifndef _PARRSB_METRICS_H_
#define _PARRSB_METRICS_H_

#include "gslib.h"

//------------------------------------------------------------------------------
// Metrics
//
typedef enum {
  RSB_COMPONENTS = 0,
  RSB_FIEDLER,
  RSB_FIEDLER_SETUP,
  RSB_FIEDLER_CALC,
  RSB_FIEDLER_CALC_NITER,
  RSB_LANCZOS_SETUP,
  RSB_LANCZOS,
  RSB_LANCZOS_TQLI,
  RSB_INVERSE_SETUP,
  RSB_PROJECT_AX,
  RSB_PROJECT_MG,
  RSB_INVERSE,
  RSB_SORT,
  RSB_PRE,
  RSB_REPAIR,
  RSB_BALANCE,
  SCHUR_PROJECT_NITER,
  SCHUR_PROJECT_OPERATOR,
  SCHUR_PROJECT_OPERATOR_FXI,
  SCHUR_PROJECT_OPERATOR_CHOL,
  SCHUR_PROJECT_OPERATOR_EZL,
  SCHUR_PROJECT_OPERATOR_MATVEC,
  SCHUR_PROJECT_PRECOND,
  SCHUR_SOLVE_CHOL1,
  SCHUR_SOLVE_CHOL2,
  SCHUR_SOLVE_PROJECT,
  SCHUR_SOLVE_SETRHS1,
  SCHUR_SOLVE_SETRHS2,
  TOL_FNL,
  TOL_TGT,
  TOL_INIT
} metric;

void metric_init();
void metric_acc(metric m, double val);
void metric_set(metric m, double val);
void metric_tic(struct comm *c, metric m);
void metric_toc(struct comm *c, metric m);
double metric_get_value(int level, metric m);
void metric_push_level();
uint metric_get_levels();
void metric_rsb_print(struct comm *c, int profile_level);
void metric_crs_print(struct comm *c, int profile_level);
void metric_finalize();

#endif
