#include "sundialsUtils.h"
#include <iostream>

using std::cout;
using std::endl;

sundialsCVODE::sundialsCVODE(unsigned int n)
	: abstol(n)
	, y0(n)
	, y(n)
{
	nEq = n;
	sundialsMem = NULL;
	findRoots = false;
	nRoots = 0;
}

sundialsCVODE::~sundialsCVODE(void)
{
	CVodeFree(&sundialsMem);
}

void sundialsCVODE::initialize(void)
{
	
	
	sundialsMem = CVodeCreate(linearMultistepMethod, nonlinearSolverMethod);
	if (check_flag((void *)sundialsMem, "CVodeCreate", 0)) {
		throw;
	}

	flag = CVodeMalloc(sundialsMem, f, t0, y.forSundials(), CV_SV, reltol,
			abstol.forSundials());
	if (check_flag(&flag, "CVodeMalloc", 1)) {
		throw;
	}

	if (findRoots) {
		rootsFound.resize(nRoots);
		// Call CVodeRootInit to specify the root function g with nRoots components
		flag = CVodeRootInit(sundialsMem, nRoots, g, theODE);
		if (check_flag(&flag, "CVodeRootInit", 1)) {
			throw;
		}
	}

	// Call CVDense to specify the CVDENSE dense linear solver
	flag = CVDense(sundialsMem, nEq);
	if (check_flag(&flag, "CVDense", 1)) {
		throw;
	}

	// Set the Jacobian routine to Jac (user-supplied)
	flag = CVDenseSetJacFn(sundialsMem, Jac, theODE);
	if (check_flag(&flag, "CVDenseSetJacFn", 1)) {
		throw;
	}

	CVodeSetFdata(sundialsMem, theODE);
}

int sundialsCVODE::integrateToTime(realtype t)
{
	flag = CVode(sundialsMem, t, y.forSundials(), &tInt, CV_NORMAL);
	return flag;
}

int sundialsCVODE::getRootInfo(void)
{
    flagr = CVodeGetRootInfo(sundialsMem, &rootsFound[0]);
    if (check_flag(&flagr, "CVodeGetRootInfo", 1)) throw;
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

  flag = CVDenseGetNumJacEvals(sundialsMem, &nje);
  check_flag(&flag, "CVDenseGetNumJacEvals", 1);
  flag = CVDenseGetNumRhsEvals(sundialsMem, &nfeLS);
  check_flag(&flag, "CVDenseGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(sundialsMem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

int sundialsCVODE::check_flag(void *flagvalue, char *funcname, int opt)
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
int sundialsCVODE::f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	// f_data contains a pointer to the "theODE" object
	return ((sdODE*) f_data)->f(t, sdVector(y), sdVector(ydot));
}

// g routine. Compute functions g_i(t,y) for i = 0,1. 
int sundialsCVODE::g(realtype t, N_Vector y, realtype *gout, void *g_data)
{
	// g_data contains a pointer to the "theODE" object
	return ((sdODE*) g_data)->g(t, sdVector(y), gout);
}

// Jacobian routine. Compute J(t,y) = df/dy. *
int sundialsCVODE::Jac(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	// jac_data contains a pointer to the "theODE" object
	return ((sdODE*) jac_data)->Jac(t, sdVector(y), sdVector(fy), sdMatrix(J));
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
	if (sundialsCVODE::check_flag((void *)v, "N_VNew_Serial", 0)) throw;
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

realtype& sdVector::operator() (unsigned int i) 
{
	return NV_Ith_S(v,i);
}

realtype sdVector::operator() (unsigned int i) const
{
	return NV_Ith_S(v,i);
}

std::ostream& operator<<(std::ostream& os, const sdVector& v)
{
	for (unsigned int i=0; i<(v.length()-1); i++) {
		os << v(i) << ", ";
	}
	os << v(v.length()-1);
	return os;
}

sdMatrix::sdMatrix(unsigned int n, unsigned int m) 
{
	alloc = true;
	M = DenseAllocMat(n,m);
}

sdMatrix::sdMatrix(DenseMat other) 
{
	alloc = false;
	M = other;
}

sdMatrix::~sdMatrix(void) {
	if (alloc) {
		DenseFreeMat(M);
	}
}

realtype& sdMatrix::operator() (unsigned int i, unsigned int j)
{
	return DENSE_ELEM(M,i,j);
}

realtype sdMatrix::operator() (unsigned int i, unsigned int j) const
{
	return DENSE_ELEM(M,i,j);
}


// Sundials IDA Solver

sundialsIDA::sundialsIDA(unsigned int n)
	: abstol(n)
	, y0(n)
	, y(n)
	, ydot(n)
	, ydot0(n)
{
	nEq = n;
	sundialsMem = NULL;
	findRoots = false;
	nRoots = 0;
}

sundialsIDA::~sundialsIDA(void)
{
	IDAFree(&sundialsMem);
}

void sundialsIDA::initialize(void)
{
	
	
	sundialsMem = IDACreate();
	if (check_flag((void *)sundialsMem, "IDACreate", 0)) {
		throw;
	}

	flag = IDAMalloc(sundialsMem, f, t0, y.forSundials(), ydot.forSundials(),
		IDA_SV, reltol, abstol.forSundials());
	if (check_flag(&flag, "IDAMalloc", 1)) {
		throw;
	}
	IDASetRdata(sundialsMem, theDAE);

	if (findRoots) {
		rootsFound.resize(nRoots);
		// Call IDARootInit to specify the root function g with nRoots components
		flag = IDARootInit(sundialsMem, nRoots, g, theDAE);
		if (check_flag(&flag, "IDARootInit", 1)) {
			throw;
		}
	}

	// Call IDADense to specify the IDADENSE dense linear solver
	flag = IDADense(sundialsMem, nEq);
	if (check_flag(&flag, "IDADense", 1)) {
		throw;
	}

	// Set the Jacobian routine to Jac (user-supplied)
	flag = IDADenseSetJacFn(sundialsMem, Jac, theDAE);
	if (check_flag(&flag, "IDADenseSetJacFn", 1)) {
		throw;
	}

	// Pick an appropriate initial condition for ydot
	flag = IDAGetConsistentIC(sundialsMem, NULL, ydot.forSundials());
	if (check_flag(&flag, "IDAGetConsistentIC", 1)) {
		throw;
	}
	std::cout << ydot << std::endl;
}

int sundialsIDA::integrateToTime(realtype t)
{
	flag = IDASolve(sundialsMem, t, &tInt, y.forSundials(),
		ydot.forSundials(), IDA_NORMAL);
	return flag;
}

int sundialsIDA::getRootInfo(void)
{
    flagr = IDAGetRootInfo(sundialsMem, &rootsFound[0]);
    if (check_flag(&flagr, "IDAGetRootInfo", 1)) throw;
    return flagr;
}

void sundialsIDA::printStats(void)
{
  long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;
  int retval;

  retval = IDAGetNumSteps(sundialsMem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(sundialsMem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDADenseGetNumJacEvals(sundialsMem, &nje);
  check_flag(&retval, "IDADenseGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(sundialsMem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(sundialsMem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(sundialsMem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDADenseGetNumResEvals(sundialsMem, &nreLS);
  check_flag(&retval, "IDADenseGetNumResEvals", 1);
  retval = IDAGetNumGEvals(sundialsMem, &nge);
  check_flag(&retval, "IDAGetNumGEvals", 1);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  printf("Number of root fn. evaluations     = %ld\n", nge);
}

int sundialsIDA::check_flag(void *flagvalue, char *funcname, int opt)
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

// f routine. Compute function f(t,y,y') = res
int sundialsIDA::f(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn, void *f_data)
{
	// f_data contains a pointer to the "theODE" object
	return ((sdDAE*) f_data)->f(t, sdVector(yIn), sdVector(ydotIn), sdVector(resIn));
}

// g routine. Compute functions g_i(t,y)
int sundialsIDA::g(realtype t, N_Vector yIn, N_Vector ydotIn, realtype *gout, void *g_data)
{
	// g_data contains a pointer to the "theODE" object
	return ((sdDAE*) g_data)->g(t, sdVector(yIn), sdVector(ydotIn), gout);
}

// Jacobian routine. Compute J(t,y) = df/dy. *
int sundialsIDA::Jac(long int N, realtype t, N_Vector yIn, N_Vector ydotIn, 
		             N_Vector res, realtype c_j, void *jac_data, DenseMat Jin, 
				     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	// jac_data contains a pointer to the "theODE" object
	return ((sdDAE*) jac_data)->Jac(t, sdVector(yIn), sdVector(ydotIn),
						        sdVector(res), c_j, sdMatrix(Jin));
}

void sundialsIDA::setDAE(sdDAE* newDAE)
{
	theDAE = newDAE;
}
