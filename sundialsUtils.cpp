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
int sundialsCVODE::Jac(long int N, DenseMat JIn, realtype t,
               N_Vector yIn, N_Vector fyIn, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	sdVector y(yIn);
	sdVector fy(fyIn);
	sdMatrix J(JIn);
	// jac_data contains a pointer to the "theODE" object
	return ((sdODE*) jac_data)->Jac(t, y, fy, J);
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

sdMatrix::sdMatrix(void)
{
	alloc = false;
	M = NULL;
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

realtype& sdMatrix::operator() (unsigned int i, unsigned int j) const
{
	return DENSE_ELEM(M,i,j);
}

// Band Matrix

sdBandMatrix::sdBandMatrix(long int N, long int bwUpper, long int bwLower, long int storeUpper)
{
	alloc = true;
	M = BandAllocMat(N,bwUpper,bwLower,storeUpper);
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
		BandFreeMat(M);
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

	IDASetRdata(sundialsMem, theDAE);

	
	// Pick an appropriate initial condition for ydot
	flag = IDASetId(sundialsMem, componentId.forSundials());

	flag = IDAMalloc(sundialsMem, f, t0, y.forSundials(), ydot.forSundials(),
		IDA_SV, reltol, abstol.forSundials());
	if (check_flag(&flag, "IDAMalloc", 1)) {
		throw;
	}

	if (findRoots) {
		rootsFound.resize(nRoots);
		// Call IDARootInit to specify the root function g with nRoots components
		flag = IDARootInit(sundialsMem, nRoots, g, theDAE);
		if (check_flag(&flag, "IDARootInit", 1)) {
			throw;
			std::cout << "IDARootInit Error" << std::endl;
		}
	}

	// Call IDASpgmr to specify the IDASpgmr dense linear solver
	flag = IDASpbcg(sundialsMem, 0);
	if (check_flag(&flag, "IDASpbcg", 1)) {
		throw;
	}

	// this seems to work better using the default J-v function rather than specifying our own...
	//flag = IDASpilsSetJacTimesVecFn(sundialsMem, JvProd, theDAE);
	//if (check_flag(&flag, "IDASpilsSetJacTimesVecFn", 1)) {
	//	throw;
	//}

	flag = IDASpilsSetPreconditioner(sundialsMem, preconditionerSetup, preconditionerSolve, theDAE);
	if (check_flag(&flag, "IDASpilsSetPreconditioner", 1)) {
		throw;
	}

	flag = IDACalcIC(sundialsMem, IDA_YA_YDP_INIT, t0+1e-4);
	if (check_flag(&flag, "IDACalcIC", 1)) {
		std::cout << "IDACalcIC Error" << std::endl;
		throw;
	}

	flag = IDAGetConsistentIC(sundialsMem, y0.forSundials(), ydot0.forSundials());
	if (check_flag(&flag, "IDAGetConsistentIC", 1)) {
		throw;
	}

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
  long int nst, nni, nje, nre, nreLS, netf, ncfn, nge, nli, npe, nps;
  int retval;

  retval = IDAGetNumSteps(sundialsMem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(sundialsMem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDAGetNumNonlinSolvIters(sundialsMem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(sundialsMem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(sundialsMem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDAGetNumGEvals(sundialsMem, &nge);
  check_flag(&retval, "IDAGetNumGEvals", 1);
  retval = IDASpilsGetNumJtimesEvals(sundialsMem, &nje);
  check_flag(&retval, "IDASpilsGetNumJtimesEvals", 1);
  retval = IDASpilsGetNumLinIters(sundialsMem, &nli);
  check_flag(&retval, "IDASpilsGetNumLinIters", 1);
  retval = IDASpilsGetNumResEvals(sundialsMem, &nreLS);
  check_flag(&retval, "IDASpilsGetNumResEvals", 1);
  retval = IDASpilsGetNumPrecEvals(sundialsMem, &npe);
  check_flag(&retval, "IDASpilsGetPrecEvals", 1);
  retval = IDASpilsGetNumPrecSolves(sundialsMem, &nps);
  check_flag(&retval, "IDASpilsGetNumPrecSolves", 1);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  printf("Number of root fn. evaluations     = %ld\n", nge);
  printf("Number of J-v Evaluations          = %ld\n", nje);
  printf("Number of linear iterations        = %ld\n", nli);
  printf("Number of preconditioner evals.    = %ld\n", npe);
  printf("Number of preconditioner solves    = %ld\n", nps);
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
	sdVector y(yIn), ydot(ydotIn), res(resIn);
	// f_data contains a pointer to the "theODE" object
	return ((sdDAE*) f_data)->f(t, y, ydot, res);
}

// g routine. Compute functions g_i(t,y)
int sundialsIDA::g(realtype t, N_Vector yIn, N_Vector ydotIn, realtype *gout, void *g_data)
{
	sdVector y(yIn), ydot(ydotIn);
	// g_data contains a pointer to the "theODE" object
	return ((sdDAE*) g_data)->g(t, y, ydot, gout);
}

// Jacobian routine. Compute J(t,y) = df/dy. *
int sundialsIDA::Jac(long int N, realtype t, N_Vector yIn, N_Vector ydotIn, 
		             N_Vector resIn, realtype c_j, void *jac_data, DenseMat Jin, 
				     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	sdVector y(yIn), ydot(ydotIn), res(resIn);
	sdMatrix J(Jin);
	// jac_data contains a pointer to the "theODE" object
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
