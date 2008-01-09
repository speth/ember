#include <iostream>
#include "sundialsUtils.h"

using std::cout;
using std::endl;

sundialsSolver::sundialsSolver(unsigned int n)
	: abstol(n)
	, y0(n)
	, y(n)
{
	nEq = n;
	cvode_mem = NULL;
	findRoots = false;
	nRoots = 0;
}

sundialsSolver::~sundialsSolver(void)
{
	CVodeFree(&cvode_mem);
}

void sundialsSolver::initialize(void)
{
	
	
	cvode_mem = CVodeCreate(linearMultistepMethod, nonlinearSolverMethod);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) {
		throw;
	}

	flag = CVodeMalloc(cvode_mem, f, t0, y.forSundials(), CV_SV, reltol,
			abstol.forSundials());
	if (check_flag(&flag, "CVodeMalloc", 1)) {
		throw;
	}

	flag = CVodeMalloc(cvode_mem, f, t0, y0.forSundials(), CV_SV, reltol,
			abstol.forSundials());
	if (check_flag(&flag, "CVodeMalloc", 1)) {
		throw;
	}

	if (findRoots) {
		rootsFound.resize(nRoots);
		// Call CVodeRootInit to specify the root function g with nRoots components
		flag = CVodeRootInit(cvode_mem, nRoots, g, theODE);
		if (check_flag(&flag, "CVodeRootInit", 1)) {
			throw;
		}
	}

	// Call CVDense to specify the CVDENSE dense linear solver
	flag = CVDense(cvode_mem, nEq);
	if (check_flag(&flag, "CVDense", 1)) {
		throw;
	}

	// Set the Jacobian routine to Jac (user-supplied)
	flag = CVDenseSetJacFn(cvode_mem, Jac, theODE);
	if (check_flag(&flag, "CVDenseSetJacFn", 1)) {
		throw;
	}

	CVodeSetFdata(cvode_mem, theODE);
}

int sundialsSolver::integrateToTime(realtype t)
{
	flag = CVode(cvode_mem, t, y.forSundials(), &tInt, CV_NORMAL);
	return flag;
}

int sundialsSolver::getRootInfo(void)
{
    flagr = CVodeGetRootInfo(cvode_mem, &rootsFound[0]);
    if (check_flag(&flagr, "CVodeGetRootInfo", 1)) throw;
    return flagr;
}

void sundialsSolver::printStats(void)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDenseGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDenseGetNumJacEvals", 1);
  flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDenseGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

sdVector::sdVector(unsigned int N) 
{
	alloc = true;
	v = N_VNew_Serial(N);
	n = N;
	if (sundialsSolver::check_flag((void *)v, "N_VNew_Serial", 0)) throw;
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


int sundialsSolver::check_flag(void *flagvalue, char *funcname, int opt)
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
int sundialsSolver::f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	// f_data contains a pointer to the "theODE" object
	return ((sdODE*) f_data)->f(t, sdVector(y), sdVector(ydot));
}

// g routine. Compute functions g_i(t,y) for i = 0,1. 
int sundialsSolver::g(realtype t, N_Vector y, realtype *gout, void *g_data)
{
	// g_data contains a pointer to the "theODE" object
	return ((sdODE*) g_data)->g(t, sdVector(y), gout);
}

// Jacobian routine. Compute J(t,y) = df/dy. *
int sundialsSolver::Jac(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	// jac_data contains a pointer to the "theODE" object
	return ((sdODE*) jac_data)->Jac(t, sdVector(y), sdVector(fy), sdMatrix(J));
}

void sundialsSolver::setODE(sdODE* newODE)
{
	theODE = newODE;
}
