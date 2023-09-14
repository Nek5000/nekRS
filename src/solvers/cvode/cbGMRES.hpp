#ifndef cbGMRES_HPP
#define cbGMRES_HPP

#ifdef ENABLE_CVODE

#include "sundials/sundials_types.h"
#include "nvector/nvector_serial.h"
#include "nvector/nvector_mpiplusx.h"

#ifdef ENABLE_CUDA
#include "nvector/nvector_cuda.h"
#endif
#ifdef ENABLE_HIP
#include "nvector/nvector_hip.h"
#endif

#endif

int cbGMRESSolve(SUNLinearSolver S, N_Vector x, N_Vector b, realtype delta);
void cbGMRESSetup(SUNLinearSolver S);

#endif
