#include "sundialsUtils.h"
#include "debugUtils.h"
#include <iostream>
#include "mathUtils.h" // debug

SundialsCvode::SundialsCvode(unsigned int n)
    : sunContext(new SundialsContext())
    , y(n, *sunContext)
    , abstol(n, *sunContext)
    , findRoots(false)
    , maxNumSteps(500)
    , minStep(0)
    , tInt(0)
    , nRoots(0)
    , errorCount(0)
    , errorStopCount(0)
    , theODE(NULL)
    , sundialsMem(NULL)
    , sundialsLinsol(NULL)
    , sundialsLinsolMatrix(NULL)
    , nEq(n)
    , bandwidth_upper(-1)
    , bandwidth_lower(-1)
    , _initialized(false)
{
}

SundialsCvode::~SundialsCvode()
{
    CVodeFree(&sundialsMem);
}

void SundialsCvode::initialize()
{
    tInt = t0;
    int flag = 0;
    if (_initialized) {
        // Starting over with a new IC, but the same ODE
        flag = CVodeReInit(sundialsMem, t0, y.forSundials());
        if (check_flag(&flag, "CVodeReInit", 1)) {
            throw DebugException("SundialsCvode::reInitialize: error in CVodeReInit");
        }

        CVodeSetMaxNumSteps(sundialsMem, maxNumSteps);
        CVodeSetMinStep(sundialsMem, minStep);
        #if EMBER_SUNDIALS_VERSION >= 66
            CVodeSetInterpolateStopTime(sundialsMem, true);
        #endif
        return;
    }

    // On the first call to initialize, need to allocate and set up
    // the Sundials solver, set tolerances and link the ODE functions
    sundialsMem = CVodeCreate(linearMultistepMethod, sunContext->get());
    if (check_flag((void *)sundialsMem, "CVodeCreate", 0)) {
        throw DebugException("SundialsCvode::initialize: error in CVodeCreate");
    }

    flag = CVodeInit(sundialsMem, f, t0, y.forSundials());
    if (check_flag(&flag, "CVodeMalloc", 1)) {
        throw DebugException("SundialsCvode::initialize: error in CVodeMalloc");
    }

    CVodeSVtolerances(sundialsMem, reltol, abstol.forSundials());
    CVodeSetUserData(sundialsMem, theODE);
    CVodeSetMaxNumSteps(sundialsMem, maxNumSteps);
    CVodeSetMinStep(sundialsMem, minStep);
    #if EMBER_SUNDIALS_VERSION >= 66
        CVodeSetInterpolateStopTime(sundialsMem, true);
    #endif

    if (findRoots) {
        rootsFound.resize(nRoots);
        // Call CVodeRootInit to specify the root function g with nRoots components
        flag = CVodeRootInit(sundialsMem, nRoots, g);
        if (check_flag(&flag, "CVodeRootInit", 1)) {
            throw DebugException("SundialsCvode::initialize: error in CVodeRootInit");
        }
    }

    if (bandwidth_upper == -1 && bandwidth_lower == -1) {
        SUNLinSolFree((SUNLinearSolver) sundialsLinsol);
        SUNMatDestroy((SUNMatrix) sundialsLinsolMatrix);
        sundialsLinsolMatrix = SUNDenseMatrix(y.size(), y.size(), sunContext->get());
        if (sundialsLinsolMatrix == nullptr) {
            throw DebugException((format("SundialsCvode::initialize: "
                "Unable to create SUNDenseMatrix of size %i x %i")
                % y.size() % y.size()).str()
            );
        }
        #if CT_SUNDIALS_USE_LAPACK
            sundialsLinsol = SUNLinSol_LapackDense(y.forSundials(),
                (SUNMatrix) sundialsLinsolMatrix, sunContext->get());
        #else
            sundialsLinsol = SUNLinSol_Dense(y.forSundials(), (SUNMatrix) sundialsLinsolMatrix,
                                        sunContext->get());
        #endif
        flag = CVodeSetLinearSolver(sundialsMem, (SUNLinearSolver) sundialsLinsolMatrix,
                                    (SUNMatrix) sundialsLinsolMatrix);

        flag = CVodeSetJacFn(sundialsMem, denseJac);
        if (check_flag(&flag, "CVodeSetJacFn", 1)) {
            throw DebugException("SundialsCvode::initialize: error in CVDlsSetDenseJacFn");
        }
    } else {
        // Specify the banded linear solver
        SUNLinSolFree((SUNLinearSolver) sundialsLinsol);
        SUNMatDestroy((SUNMatrix) sundialsLinsolMatrix);
        sundialsLinsolMatrix = SUNBandMatrix(y.size(), bandwidth_upper,
                                             bandwidth_lower, sunContext->get());
        if (sundialsLinsolMatrix == nullptr) {
            throw DebugException((format("SundialsCvode::initialize: "
                "Unable to create SUNBandMatrix of size %i with bandwidths "
                "%i and %i") % y.size() % bandwidth_upper % bandwidth_lower).str()
            );
        }
        #if CT_SUNDIALS_USE_LAPACK
            sundialsLinsol = SUNLinSol_LapackBand(y.forSundials(),
                (SUNMatrix) sundialsLinsolMatrix, sunContext->get());
        #else
            sundialsLinsol = SUNLinSol_Band(y.forSundials(),
                (SUNMatrix) sundialsLinsolMatrix, sunContext->get());
        #endif
        if (sundialsLinsol == nullptr) {
            throw DebugException("SundialsCvode::initialize: "
                "Unable to create Sundials Band matrix solver"
            );
        }
        CVodeSetLinearSolver(sundialsMem, (SUNLinearSolver) sundialsLinsol,
                             (SUNMatrix) sundialsLinsolMatrix);
    }

    _initialized = true;
}

bool SundialsCvode::initialized() const
{
    return _initialized;
}

void SundialsCvode::setBandwidth(int upper, int lower)
{
    bandwidth_upper = upper;
    bandwidth_lower = lower;
}

int SundialsCvode::integrateToTime(realtype t)
{
    assert(mathUtils::notnan(y));
    int flag = CVode(sundialsMem, t, y.forSundials(), &tInt, CV_NORMAL);
    if (flag != CV_SUCCESS) {
        errorCount += 1;
        if (errorCount > errorStopCount) {
            throw DebugException("CVODE Integrator had too many errors");
        }
    }
    return flag;
}

int SundialsCvode::integrateOneStep(realtype tf)
{
    assert(mathUtils::notnan(y));
    CVodeSetStopTime(sundialsMem, tf);
    int flag = CVode(sundialsMem, tf, y.forSundials(), &tInt, CV_ONE_STEP);
    if (flag != CV_SUCCESS && flag != CV_TSTOP_RETURN) {
        errorCount += 1;
        if (errorCount > errorStopCount) {
            throw DebugException("CVODE Integrator had too many errors");
        }
    }
    return flag;
}

int SundialsCvode::getRootInfo()
{
    int flag = CVodeGetRootInfo(sundialsMem, &rootsFound[0]);
    if (check_flag(&flag, "CVodeGetRootInfo", 1)) {
        throw DebugException("SundialsCvode::getRootInfo: error in CVodeGetRootInfo");
    }
    return flag;
}

void SundialsCvode::printStats()
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(sundialsMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(sundialsMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(sundialsMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(sundialsMem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(sundialsMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(sundialsMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  #if SUNDIALS_VERSION_MAJOR >= 6
      flag = CVodeGetNumJacEvals(sundialsMem, &nje);
  #else
      flag = CVDlsGetNumJacEvals(sundialsMem, &nje);
  #endif
  check_flag(&flag, "CVodeGetNumJacEvals", 1);
  #if SUNDIALS_VERSION_MAJOR >= 6
      flag = CVodeGetNumRhsEvals(sundialsMem, &nfeLS);
  #else
      flag = CVDlsGetNumRhsEvals(sundialsMem, &nfeLS);
  #endif
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(sundialsMem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
     nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
     nni, ncfn, netf, nge);
}

long int SundialsCvode::getNumSteps()
{
    long int n;
    CVodeGetNumSteps(sundialsMem, &n);
    return n;
}

int SundialsCvode::getLastOrder()
{
    int n;
    CVodeGetLastOrder(sundialsMem, &n);
    return n;
}

realtype SundialsCvode::getLastStep()
{
    realtype h;
    CVodeGetLastStep(sundialsMem, &h);
    return h;
}

int SundialsCvode::check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    return(1); }

  return(0);
}

// f routine. Compute function f(t,y).
int SundialsCvode::f(realtype t, N_Vector yIn, N_Vector ydotIn, void *f_data)
{
    sdODE* ode = (sdODE*) f_data;
    sdVector y(yIn);
    sdVector ydot(ydotIn);
    // f_data contains a pointer to the "theODE" object
    return ode->f(t, y, ydot);
}

// g routine. Compute functions g_i(t,y) for i = 0,1.
int SundialsCvode::g(realtype t, N_Vector yIn, realtype *gout, void *g_data)
{
    sdODE* ode = (sdODE*) g_data;
    sdVector y(yIn);
    // g_data contains a pointer to the "theODE" object
    return ode->g(t, y, gout);
}

// Jacobian routine. Compute J(t,y) = df/dy. *
int SundialsCvode::denseJac(
    realtype t, N_Vector yIn,
    N_Vector fyIn, SUNMatrix JIn, void *user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdODE* ode = (sdODE*) user_data;
    sdVector y(yIn);
    sdVector fy(fyIn);
    sdMatrix J(JIn);
    // user_data contains a pointer to the "theODE" object
    return ode->denseJacobian(t, y, fy, J);
}

int SundialsCvode::bandJac(
	realtype t, N_Vector yIn, N_Vector fyIn, SUNMatrix JIn, void* user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdODE* ode = (sdODE*) user_data;
    sdVector y(yIn);
    sdVector fy(fyIn);
    sdBandMatrix J(JIn);
    // user_data contains a pointer to the "theODE" object
    return ode->bandedJacobian(t, y, fy, J);
}

void SundialsCvode::setODE(sdODE* newODE)
{
    theODE = newODE;
}

sdVector::sdVector()
{
    v = nullptr;
    alloc = false;
    n = 0;
}

sdVector::sdVector(unsigned int N, SundialsContext& context)
{
    alloc = true;
    v = N_VNew_Serial(N, context.get());
    n = N;
    if (SundialsCvode::check_flag((void *)v, "N_VNew_Serial", 0)) {
        throw DebugException("sdVector: error allocating vector");
    }
}

sdVector::sdVector(N_Vector other)
{
    alloc = false;
    v = other;
    n = NV_LENGTH_S(v);
}

sdVector::sdVector(const sdVector& other)
{
    alloc = false; // This is bad if the copy persists longer than the original...
    v = other.v;
    n = other.n;
}

sdVector::~sdVector()
{
    if (alloc) {
        N_VDestroy_Serial(v);
    }
}

std::ostream& operator<<(std::ostream& os, const sdVector& v)
{
    for (unsigned int i=0; i<(v.length()-1); i++) {
        os << v[i] << ", ";
    }
    os << v[v.length()-1];
    return os;
}

sdMatrix::sdMatrix(unsigned int n, unsigned int m, SundialsContext& context)
{
    alloc = true;
    M = SUNDenseMatrix(n,m, context.get());
}

sdMatrix::sdMatrix(SUNMatrix other)
{
    alloc = false;
    M = other;
}

sdMatrix::sdMatrix()
{
    alloc = false;
    M = NULL;
}

sdMatrix::~sdMatrix() {
    if (alloc) {
        SUNMatDestroy(M);
    }
}

// Band Matrix

sdBandMatrix::sdBandMatrix(long int N, long int bwUpper, long int bwLower,
                           SundialsContext& context)
{
    alloc = true;
    M = SUNBandMatrix(N, bwUpper, bwLower, context.get());
}

sdBandMatrix::sdBandMatrix(SUNMatrix other)
{
    alloc = false;
    M = other;
}

sdBandMatrix::sdBandMatrix()
{
    alloc = false;
    M = NULL;
}

sdBandMatrix::~sdBandMatrix() {
    if (alloc) {
        SUNMatDestroy(M);
    }
}

void sdBandMatrix::print(const std::string& name) const
{
    std::cout << SM_ROWS_B(M) << std::endl;
    std::cout << SM_LBAND_B(M) << std::endl;
    std::cout << SM_UBAND_B(M) << std::endl;
    for (int i = 0; i < SM_ROWS_B(M); i++) {
        for (int j = i - SM_LBAND_B(M); j <= i + SM_UBAND_B(M); j++) {
            if (j < 0 || j >= SM_COLUMNS_B(M)) {
                continue;
            }
            std::cout << boost::format("%s[%i,%i] = %e") % name % i % j % SM_ELEMENT_B(M,i,j) << std::endl;
        }
    }
    std::cout.flush();
}

// Sundials IDA Solver

#if EMBER_ENABLE_IDA
SundialsIda::SundialsIda(unsigned int n)
    : abstol(n)
    , findRoots(false)
    , calcIC(false)
    , imposeConstraints(false)
    , y(n)
    , ydot(n)
    , y0(n)
    , ydot0(n)
    , componentId(n)
    , constraints(n)
    , nRoots(0)
    , theDAE(NULL)
    , sundialsMem(NULL)
    , nEq(n)
{
}

SundialsIda::~SundialsIda()
{
    IDAFree(&sundialsMem);
}

void SundialsIda::initialize()
{
    tInt = t0;
    sundialsMem = IDACreate();
    if (check_flag((void *)sundialsMem, "IDACreate", 0)) {
        throw DebugException("SundialsIda::initialize: error in IDACreate");
    }

    IDASetUserData(sundialsMem, theDAE);

    int flag;
    if (calcIC)    {
        // Pick an appropriate initial condition for ydot and algebraic components of y
        flag = IDASetId(sundialsMem, componentId.forSundials());
    }

    flag = IDAInit(sundialsMem, f, t0, y.forSundials(), ydot.forSundials());
    if (check_flag(&flag, "IDAMalloc", 1)) {
        throw DebugException("SundialsIda::initialize: error in IDAInit");
    }

    IDASVtolerances(sundialsMem, reltol, abstol.forSundials());

    if (findRoots) {
        rootsFound.resize(nRoots);
        // Call IDARootInit to specify the root function g with nRoots components
        flag = IDARootInit(sundialsMem, nRoots, g);
        if (check_flag(&flag, "IDARootInit", 1)) {
            throw DebugException("SundialsIda::initialize: error in IDARootInit");
        }
    }

    // Call IDASpbcg to specify the IDASpbcg dense linear solver
    flag = IDASpbcg(sundialsMem, 0);
    if (check_flag(&flag, "IDASpbcg", 1)) {
        throw DebugException("SundialsIda::initialize: error in IDASpbcg");
    }

    if (imposeConstraints) {
        flag = IDASetConstraints(sundialsMem, constraints.forSundials());
        if (check_flag(&flag, "IDASetConstraints", 1)) {
            throw DebugException("SundialsIda::initialize: error in IDASetConstraints");
        }
    }

    // this seems to work better using the default J-v function rather than specifying our own...
    //flag = IDASpilsSetJacTimesVecFn(sundialsMem, JvProd, theDAE);
    //if (check_flag(&flag, "IDASpilsSetJacTimesVecFn", 1)) {
    //    throw myException("SundialsIda::initialize: error in IDASpilsSetJacTimesVecFn");
    //}

    flag = IDASpilsSetPreconditioner(sundialsMem, preconditionerSetup, preconditionerSolve);
    if (check_flag(&flag, "IDASpilsSetPreconditioner", 1)) {
        throw DebugException("SundialsIda::initialize: error in IDASpilsSetPreconditioner");
    }

    if (calcIC) {
        flag = IDACalcIC(sundialsMem, IDA_YA_YDP_INIT, t0+1e-4);
        if (check_flag(&flag, "IDACalcIC", 1)) {
            logFile.write("IDACalcIC Error");
            throw DebugException("SundialsIda::initialize: error in IDACalcIC");
        }

        flag = IDAGetConsistentIC(sundialsMem, y0.forSundials(), ydot0.forSundials());
        if (check_flag(&flag, "IDAGetConsistentIC", 1)) {
            throw DebugException("SundialsIda::initialize: error in IDAGetConsistentIC");
        }
    }

}

int SundialsIda::integrateToTime(realtype t)
{
    int flag = IDASolve(sundialsMem, t, &tInt, y.forSundials(),
        ydot.forSundials(), IDA_NORMAL);
    return flag;
}

int SundialsIda::integrateOneStep()
{
    int flag = IDASolve(sundialsMem, tInt+1, &tInt, y.forSundials(),
        ydot.forSundials(), IDA_ONE_STEP);
    return flag;
}

int SundialsIda::getRootInfo()
{
    int flag = IDAGetRootInfo(sundialsMem, &rootsFound[0]);
    if (check_flag(&flag, "IDAGetRootInfo", 1))
        throw DebugException("SundialsIda::getRootInfo: error in IDAGetRootInfo.");
    return flag;
}

void SundialsIda::printStats(clock_t dt)
{
    long int nst, nje, nre, netf, ncfn, nge, npe, nps;
    int retval;

    retval = IDAGetNumSteps(sundialsMem, &nst);
    check_flag(&retval, "IDAGetNumSteps", 1);
    retval = IDAGetNumResEvals(sundialsMem, &nre);
    check_flag(&retval, "IDAGetNumResEvals", 1);
    retval = IDAGetNumErrTestFails(sundialsMem, &netf);
    check_flag(&retval, "IDAGetNumErrTestFails", 1);
    retval = IDAGetNumNonlinSolvConvFails(sundialsMem, &ncfn);
    check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
    if (findRoots) {
        retval = IDAGetNumGEvals(sundialsMem, &nge);
        check_flag(&retval, "IDAGetNumGEvals", 1);
    }
    retval = IDASpilsGetNumJtimesEvals(sundialsMem, &nje);
    check_flag(&retval, "IDASpilsGetNumJtimesEvals", 1);
    retval = IDASpilsGetNumPrecEvals(sundialsMem, &npe);
    check_flag(&retval, "IDASpilsGetPrecEvals", 1);
    retval = IDASpilsGetNumPrecSolves(sundialsMem, &nps);
    check_flag(&retval, "IDASpilsGetNumPrecSolves", 1);

    printf("\nIDA Solver Statistics: \n\n");
    printf("Number of steps                    = %ld\n", nst);
    printf("Number of residual evaluations     = %ld\n", nre);
    printf("Number of error test failures      = %ld\n", netf);
    printf("Number of nonlinear conv. failures = %ld\n", ncfn);
    if (findRoots) {
        printf("Number of root fn. evaluations     = %ld\n", nge);
    }
    printf("Number of J-v Evaluations          = %ld\n", nje);
    printf("Number of preconditioner evals.    = %ld\n", npe);
    printf("Number of preconditioner solves    = %ld\n", nps);

    retval = IDAGetNumSteps(sundialsMem, &nst);
    if (dt == 0) {
        logFile.write(format("IDA solver took %i steps.\n") % nst);
    } else {
        logFile.write(format("IDA solver took %i steps in %f seconds.\n") %
                nst % (((double) dt)/CLOCKS_PER_SEC));
    }
}

int SundialsIda::check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}

double SundialsIda::getStepSize()
{
    double dt;
    IDAGetCurrentStep(sundialsMem, &dt);
    return dt;
}

void SundialsIda::setInitialStepSize(double dt)
{
    IDASetInitStep(sundialsMem, dt);
}

void SundialsIda::setMaxStepSize(double dt)
{
    IDASetMaxStep(sundialsMem, dt);
}

int SundialsIda::getLastOrder()
{
    int order;
    IDAGetLastOrder(sundialsMem, &order);
    return order;
}

void SundialsIda::disableErrorOutput()
{
    IDASetErrFile(sundialsMem, NULL);
}

// f routine. Compute function f(t,y,y') = res
int SundialsIda::f(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn, void *f_data)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn);
    // f_data is a pointer to the "theDAE" object
    return ((sdDAE*) f_data)->f(t, y, ydot, res);
}

// g routine. Compute functions g_i(t,y)
int SundialsIda::g(realtype t, N_Vector yIn, N_Vector ydotIn, realtype *gout, void *g_data)
{
    sdVector y(yIn), ydot(ydotIn);
    // g_data is a pointer to the "theDAE" object
    return ((sdDAE*) g_data)->g(t, y, ydot, gout);
}

int SundialsIda::Jac(long int N, realtype t, N_Vector yIn, N_Vector ydotIn,
                     N_Vector resIn, realtype c_j, void *jac_data, DenseMat Jin,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn);
    sdMatrix J(Jin);
    // jac_data is a pointer to the "theDAE" object
    return ((sdDAE*) jac_data)->Jac(t, y, ydot, res, c_j, J);
}

int SundialsIda::JvProd(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn,
                        N_Vector vIn, N_Vector JvIn, realtype c_j, void* jac_data,
                        N_Vector tmp1, N_Vector tmp2)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn), v(vIn), Jv(JvIn);
    return ((sdDAE*) jac_data)->JvProd(t,y , ydot, res, v, Jv , c_j);
}

int SundialsIda::preconditionerSetup(realtype t, N_Vector yIn, N_Vector ydotIn,
                                     N_Vector resIn, realtype c_j, void* p_data,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn);
    return ((sdDAE*) p_data)->preconditionerSetup(t, y, ydot, res, c_j);
}

int SundialsIda::preconditionerSolve(realtype t, N_Vector yIn, N_Vector ydotIn,
                                     N_Vector resIn, N_Vector rhsIn, N_Vector vIn,
                                     realtype c_j, realtype delta, void* p_data,
                                     N_Vector tmp)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn), rhs(rhsIn), v(vIn);
    return ((sdDAE*) p_data)->preconditionerSolve(t, y, ydot, res, rhs, v, c_j, delta);
}

void SundialsIda::setDAE(sdDAE* newDAE)
{
    theDAE = newDAE;
}
#endif // EMBER_ENABLE_IDA
