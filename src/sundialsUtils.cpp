#include "sundialsUtils.h"
#include "debugUtils.h"
#include <iostream>
#include "mathUtils.h" // debug

sundialsCVODE::sundialsCVODE(unsigned int n)
    : y(n)
    , abstol(n)
    , tInt(0)
    , bandwidth_upper(-1)
    , bandwidth_lower(-1)
    , _initialized(false)
{
    nEq = n;
    sundialsMem = NULL;
    findRoots = false;
    nRoots = 0;
    maxNumSteps = 500;
    minStep = 0;
    errorCount = 0;
    errorStopCount = 0;
}

sundialsCVODE::~sundialsCVODE(void)
{
    CVodeFree(&sundialsMem);
}

void sundialsCVODE::initialize()
{
    theODE->initialize();

    if (_initialized) {
        // Starting over with a new IC, but the same ODE
        flag = CVodeReInit(sundialsMem, t0, y.forSundials());
        if (check_flag(&flag, "CVodeReInit", 1)) {
            throw debugException("sundialsCVODE::reInitialize: error in CVodeReInit");
        }

        CVodeSetMaxNumSteps(sundialsMem, maxNumSteps);
        CVodeSetMinStep(sundialsMem, minStep);
        return;
    }

    // On the first call to initialize, need to allocate and set up
    // the Sundials solver, set tolerances and link the ODE functions
    sundialsMem = CVodeCreate(linearMultistepMethod, nonlinearSolverMethod);
    if (check_flag((void *)sundialsMem, "CVodeCreate", 0)) {
        throw debugException("sundialsCVODE::initialize: error in CVodeCreate");
    }

    flag = CVodeInit(sundialsMem, f, t0, y.forSundials());
    if (check_flag(&flag, "CVodeMalloc", 1)) {
        throw debugException("sundialsCVODE::initialize: error in CVodeMalloc");
    }

    CVodeSVtolerances(sundialsMem, reltol, abstol.forSundials());
    CVodeSetUserData(sundialsMem, theODE);
    CVodeSetMaxNumSteps(sundialsMem, maxNumSteps);
    CVodeSetMinStep(sundialsMem, minStep);

    if (findRoots) {
        rootsFound.resize(nRoots);
        // Call CVodeRootInit to specify the root function g with nRoots components
        flag = CVodeRootInit(sundialsMem, nRoots, g);
        if (check_flag(&flag, "CVodeRootInit", 1)) {
            throw debugException("sundialsCVODE::initialize: error in CVodeRootInit");
        }
    }

    if (bandwidth_upper == -1 && bandwidth_lower == -1) {
        // Call CVDense to specify the CVDENSE dense linear solver
        flag = CVDense(sundialsMem, nEq);
        if (check_flag(&flag, "CVDense", 1)) {
            throw debugException("sundialsCVODE::initialize: error in CVDense");
        }

        // Set the Jacobian routine to denseJac (user-supplied)
        flag = CVDlsSetDenseJacFn(sundialsMem, denseJac);
        if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) {
            throw debugException("sundialsCVODE::initialize: error in CVDlsSetDenseJacFn");
        }
    } else {
        // Call CVDense to specify the CVBAND Banded linear solver
        flag = CVBand(sundialsMem, nEq, bandwidth_upper, bandwidth_lower);
        if (check_flag(&flag, "CVBand", 1)) {
            throw debugException("sundialsCVODE::initialize: error in CVBand");
        }

        // Set the Jacobian routine to bandJac (user-supplied)
        flag = CVDlsSetBandJacFn(sundialsMem, bandJac);
        if (check_flag(&flag, "CVDlsSetBandJacFn", 1)) {
            throw debugException("sundialsCVODE::initialize: error in CVDlsSetBandJacFn");
        }
    }

    _initialized = true;
}

bool sundialsCVODE::initialized() const
{
    return _initialized;
}

void sundialsCVODE::setBandwidth(int upper, int lower)
{
    bandwidth_upper = upper;
    bandwidth_lower = lower;
}

int sundialsCVODE::integrateToTime(realtype t)
{
    flag = CVode(sundialsMem, t, y.forSundials(), &tInt, CV_NORMAL);
    if (flag != CV_SUCCESS) {
        errorCount += 1;
        if (errorCount > errorStopCount) {
            throw debugException("CVODE Integrator had too many errors");
        }
    }
    return flag;
}

int sundialsCVODE::integrateOneStep(realtype tf)
{
    CVodeSetStopTime(sundialsMem, tf);
    flag = CVode(sundialsMem, tf, y.forSundials(), &tInt, CV_ONE_STEP);
    if (flag != CV_SUCCESS && flag != CV_TSTOP_RETURN) {
        errorCount += 1;
        if (errorCount > errorStopCount) {
            throw debugException("CVODE Integrator had too many errors");
        }
    }
    return flag;
}

int sundialsCVODE::getRootInfo(void)
{
    flagr = CVodeGetRootInfo(sundialsMem, &rootsFound[0]);
    if (check_flag(&flagr, "CVodeGetRootInfo", 1)) {
        throw debugException("sundialsCVODE::getRootInfo: error in CVodeGetRootInfo");
    }
    return flagr;
}

void sundialsCVODE::printStats(void)
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

  flag = CVDlsGetNumJacEvals(sundialsMem, &nje);
  check_flag(&flag, "CVDenseGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(sundialsMem, &nfeLS);
  check_flag(&flag, "CVDenseGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(sundialsMem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
     nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
     nni, ncfn, netf, nge);
}

long int sundialsCVODE::getNumSteps()
{
    long int n;
    CVodeGetNumSteps(sundialsMem, &n);
    return n;
}

int sundialsCVODE::getLastOrder()
{
    int n;
    CVodeGetLastOrder(sundialsMem, &n);
    return n;
}

realtype sundialsCVODE::getLastStep()
{
    realtype h;
    CVodeGetLastStep(sundialsMem, &h);
    return h;
}

int sundialsCVODE::check_flag(void *flagvalue, const char *funcname, int opt)
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
int sundialsCVODE::f(realtype t, N_Vector yIn, N_Vector ydotIn, void *f_data)
{
    sdVector y(yIn);
    sdVector ydot(ydotIn);
    // f_data contains a pointer to the "theODE" object
    return ((sdODE*) f_data)->f(t, y, ydot);
}

// g routine. Compute functions g_i(t,y) for i = 0,1.
int sundialsCVODE::g(realtype t, N_Vector yIn, realtype *gout, void *g_data)
{
    sdVector y(yIn);
    // g_data contains a pointer to the "theODE" object
    return ((sdODE*) g_data)->g(t, y, gout);
}

// Jacobian routine. Compute J(t,y) = df/dy. *
int sundialsCVODE::denseJac(int N, realtype t, N_Vector yIn,
                            N_Vector fyIn, DenseMat JIn, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdVector y(yIn);
    sdVector fy(fyIn);
    sdMatrix J(JIn);
    // user_data contains a pointer to the "theODE" object
    return ((sdODE*) user_data)->denseJacobian(t, y, fy, J);
}

int sundialsCVODE::bandJac(int N, int mupper, int mLower, realtype t,
        N_Vector yIn, N_Vector fyIn, DlsMat JIn, void* user_data,
        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdVector y(yIn);
    sdVector fy(fyIn);
    sdBandMatrix J(JIn);
    // user_data contains a pointer to the "theODE" object
    return ((sdODE*) user_data)->bandedJacobian(t, y, fy, J);
}

void sundialsCVODE::setODE(sdODE* newODE)
{
    theODE = newODE;
}

sdVector::sdVector(unsigned int N)
{
    alloc = true;
    v = N_VNew_Serial(N);
    n = N;
    if (sundialsCVODE::check_flag((void *)v, "N_VNew_Serial", 0)) {
        throw debugException("sdVector: error allocating vector");
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

sdVector::~sdVector(void)
{
    if (alloc) {
        N_VDestroy_Serial(v);
    }
}

realtype& sdVector::operator[] (unsigned int i)
{
    return NV_Ith_S(v,i);
}

realtype& sdVector::operator[] (unsigned int i) const
{
    return NV_Ith_S(v,i);
}

std::ostream& operator<<(std::ostream& os, const sdVector& v)
{
    for (unsigned int i=0; i<(v.length()-1); i++) {
        os << v[i] << ", ";
    }
    os << v[v.length()-1];
    return os;
}

sdMatrix::sdMatrix(unsigned int n, unsigned int m)
{
    alloc = true;
    M = NewDenseMat(n,m);
}

sdMatrix::sdMatrix(DenseMat other)
{
    alloc = false;
    M = other;
}

sdMatrix::sdMatrix(void)
{
    alloc = false;
    M = NULL;
}

sdMatrix::~sdMatrix(void) {
    if (alloc) {
    	DestroyMat(M);
    }
}

realtype& sdMatrix::operator() (unsigned int i, unsigned int j)
{
    return DENSE_ELEM(M,i,j);
}

realtype& sdMatrix::operator() (unsigned int i, unsigned int j) const
{
    return DENSE_ELEM(M,i,j);
}

// Band Matrix

sdBandMatrix::sdBandMatrix(long int N, long int bwUpper, long int bwLower)
{
    alloc = true;
    M = NewBandMat(N,bwUpper,bwLower, std::min(N-1, bwUpper+bwLower));
}

sdBandMatrix::sdBandMatrix(BandMat other)
{
    alloc = false;
    M = other;
}

sdBandMatrix::sdBandMatrix(void)
{
    alloc = false;
    M = NULL;
}

sdBandMatrix::~sdBandMatrix(void) {
    if (alloc) {
        DestroyMat(M);
    }
}

realtype& sdBandMatrix::operator() (long int i, long int j)
{
    return BAND_ELEM(M,i,j);
}

realtype& sdBandMatrix::operator() (long int i, long int j) const
{
    return BAND_ELEM(M,i,j);
}

// Sundials IDA Solver

sundialsIDA::sundialsIDA(unsigned int n)
    : abstol(n)
    , y(n)
    , ydot(n)
    , y0(n)
    , ydot0(n)
    , componentId(n)
    , constraints(n)
{
    nEq = n;
    sundialsMem = NULL;
    findRoots = false;
    nRoots = 0;
    imposeConstraints = false;
    calcIC = false;
}

sundialsIDA::~sundialsIDA(void)
{
    IDAFree(&sundialsMem);
}

void sundialsIDA::initialize(void)
{
    sundialsMem = IDACreate();
    if (check_flag((void *)sundialsMem, "IDACreate", 0)) {
        throw debugException("sundialsIDA::initialize: error in IDACreate");
    }

    IDASetUserData(sundialsMem, theDAE);

    if (calcIC)    {
        // Pick an appropriate initial condition for ydot and algebraic components of y
        flag = IDASetId(sundialsMem, componentId.forSundials());
    }

    flag = IDAInit(sundialsMem, f, t0, y.forSundials(), ydot.forSundials());
    if (check_flag(&flag, "IDAMalloc", 1)) {
        throw debugException("sundialsIDA::initialize: error in IDAInit");
    }

    IDASVtolerances(sundialsMem, reltol, abstol.forSundials());

    if (findRoots) {
        rootsFound.resize(nRoots);
        // Call IDARootInit to specify the root function g with nRoots components
        flag = IDARootInit(sundialsMem, nRoots, g);
        if (check_flag(&flag, "IDARootInit", 1)) {
            throw debugException("sundialsIDA::initialize: error in IDARootInit");
        }
    }

    // Call IDASpbcg to specify the IDASpbcg dense linear solver
    flag = IDASpbcg(sundialsMem, 0);
    if (check_flag(&flag, "IDASpbcg", 1)) {
        throw debugException("sundialsIDA::initialize: error in IDASpbcg");
    }

    if (imposeConstraints) {
        flag = IDASetConstraints(sundialsMem, constraints.forSundials());
        if (check_flag(&flag, "IDASetConstraints", 1)) {
            throw debugException("sundialsIDA::initialize: error in IDASetConstraints");
        }
    }

    // this seems to work better using the default J-v function rather than specifying our own...
    //flag = IDASpilsSetJacTimesVecFn(sundialsMem, JvProd, theDAE);
    //if (check_flag(&flag, "IDASpilsSetJacTimesVecFn", 1)) {
    //    throw myException("sundialsIDA::initialize: error in IDASpilsSetJacTimesVecFn");
    //}

    flag = IDASpilsSetPreconditioner(sundialsMem, preconditionerSetup, preconditionerSolve);
    if (check_flag(&flag, "IDASpilsSetPreconditioner", 1)) {
        throw debugException("sundialsIDA::initialize: error in IDASpilsSetPreconditioner");
    }

    if (calcIC) {
        flag = IDACalcIC(sundialsMem, IDA_YA_YDP_INIT, t0+1e-4);
        if (check_flag(&flag, "IDACalcIC", 1)) {
            logFile.write("IDACalcIC Error");
            throw debugException("sundialsIDA::initialize: error in IDACalcIC");
        }

        flag = IDAGetConsistentIC(sundialsMem, y0.forSundials(), ydot0.forSundials());
        if (check_flag(&flag, "IDAGetConsistentIC", 1)) {
            throw debugException("sundialsIDA::initialize: error in IDAGetConsistentIC");
        }
    }

}

int sundialsIDA::integrateToTime(realtype t)
{
    flag = IDASolve(sundialsMem, t, &tInt, y.forSundials(),
        ydot.forSundials(), IDA_NORMAL);
    return flag;
}

int sundialsIDA::integrateOneStep(void)
{
    flag = IDASolve(sundialsMem, tInt+1, &tInt, y.forSundials(),
        ydot.forSundials(), IDA_ONE_STEP);
    return flag;
}

int sundialsIDA::getRootInfo(void)
{
    flagr = IDAGetRootInfo(sundialsMem, &rootsFound[0]);
    if (check_flag(&flagr, "IDAGetRootInfo", 1))
        throw debugException("sundialsIDA::getRootInfo: error in IDAGetRootInfo.");
    return flagr;
}

void sundialsIDA::printStats(clock_t dt)
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

int sundialsIDA::check_flag(void *flagvalue, const char *funcname, int opt)
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

double sundialsIDA::getStepSize(void)
{
    double dt;
    IDAGetCurrentStep(sundialsMem, &dt);
    return dt;
}

void sundialsIDA::setInitialStepSize(double dt)
{
    IDASetInitStep(sundialsMem, dt);
}

void sundialsIDA::setMaxStepSize(double dt)
{
    IDASetMaxStep(sundialsMem, dt);
}

int sundialsIDA::getLastOrder(void)
{
    int order;
    IDAGetLastOrder(sundialsMem, &order);
    return order;
}

void sundialsIDA::disableErrorOutput(void)
{
    IDASetErrFile(sundialsMem, NULL);
}

// f routine. Compute function f(t,y,y') = res
int sundialsIDA::f(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn, void *f_data)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn);
    // f_data is a pointer to the "theDAE" object
    return ((sdDAE*) f_data)->f(t, y, ydot, res);
}

// g routine. Compute functions g_i(t,y)
int sundialsIDA::g(realtype t, N_Vector yIn, N_Vector ydotIn, realtype *gout, void *g_data)
{
    sdVector y(yIn), ydot(ydotIn);
    // g_data is a pointer to the "theDAE" object
    return ((sdDAE*) g_data)->g(t, y, ydot, gout);
}

// Jacobian routine. Compute J(t,y) = df/dy.
int sundialsIDA::Jac(long int N, realtype t, N_Vector yIn, N_Vector ydotIn,
                     N_Vector resIn, realtype c_j, void *jac_data, DenseMat Jin,
                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn);
    sdMatrix J(Jin);
    // jac_data is a pointer to the "theDAE" object
    return ((sdDAE*) jac_data)->Jac(t, y, ydot, res, c_j, J);
}

int sundialsIDA::JvProd(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn,
                        N_Vector vIn, N_Vector JvIn, realtype c_j, void* jac_data,
                        N_Vector tmp1, N_Vector tmp2)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn), v(vIn), Jv(JvIn);
    return ((sdDAE*) jac_data)->JvProd(t,y , ydot, res, v, Jv , c_j);
}

int sundialsIDA::preconditionerSetup(realtype t, N_Vector yIn, N_Vector ydotIn,
                                     N_Vector resIn, realtype c_j, void* p_data,
                                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn);
    return ((sdDAE*) p_data)->preconditionerSetup(t, y, ydot, res, c_j);
}

int sundialsIDA::preconditionerSolve(realtype t, N_Vector yIn, N_Vector ydotIn,
                                     N_Vector resIn, N_Vector rhsIn, N_Vector vIn,
                                     realtype c_j, realtype delta, void* p_data,
                                     N_Vector tmp)
{
    sdVector y(yIn), ydot(ydotIn), res(resIn), rhs(rhsIn), v(vIn);
    return ((sdDAE*) p_data)->preconditionerSolve(t, y, ydot, res, rhs, v, c_j, delta);
}

void sundialsIDA::setDAE(sdDAE* newDAE)
{
    theDAE = newDAE;
}
