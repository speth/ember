#pragma once

#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_band.h>

#include <sundials/sundials_types.h>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <ida/ida.h>
//#include <ida/ida_dense.h>
#include <ida/ida_spbcgs.h>

#include <iostream>
#include <vector>

// wrapper class for Sundials "N_Vector"
class sdVector {
public:
	sdVector(unsigned int n);
	sdVector(N_Vector other);
	sdVector(const sdVector& other);
	~sdVector(void);

	realtype& operator()(unsigned int i);
	realtype operator()(unsigned int i) const;
	N_Vector& forSundials(void) {return v;}
	unsigned int length(void) const {return n;}

private:
	N_Vector v;
	bool alloc;
	unsigned int n;
};

std::ostream& operator<<(std::ostream& os, const sdVector& v);

// wrapper class for Sundials DenseMat
class sdMatrix 
{
public:
	sdMatrix(unsigned int n, unsigned int m);
	sdMatrix(DenseMat other);
	sdMatrix(void);
	~sdMatrix(void);
	realtype& operator()(unsigned int i, unsigned int j);
	realtype& operator()(unsigned int i, unsigned int j) const;
	realtype** forSundials(void) {return M->data;}

private:
	DenseMat M;
	bool alloc;
};

class sdBandMatrix
{
public:
	sdBandMatrix(long int N, long int bwUpper, long int bwLower, long int storeUpper);
	sdBandMatrix(BandMat other);
	sdBandMatrix(void);
	~sdBandMatrix(void);
	realtype& operator()(long int i, long int j);
	realtype& operator()(long int i, long int j) const;
	BandMat& forSundials(void) {return M;}
	
private:
	BandMat M;
	bool alloc;
};

class sdODE {
public:
	virtual int f(realtype t, sdVector& y, sdVector& ydot)=0;
	virtual int g(realtype t, sdVector& y, realtype* gOut) {return 0;}
	virtual int Jac(realtype t, sdVector& y, sdVector& fy, sdMatrix& J)=0;
	virtual ~sdODE(void) {}
};

// wrapper class for the Sundials CVODE solver
class sundialsCVODE
{
public:
	sundialsCVODE(unsigned int n);
	~sundialsCVODE(void);

	void initialize(void);
	int integrateToTime(realtype t);
	void setODE(sdODE* newODE);
	int getRootInfo(); 
	void printStats();
	
	// Check function return value...
	//   opt == 0 means SUNDIALS function allocates memory so check if
	//            returned NULL pointer
	//   opt == 1 means SUNDIALS function returns a flag so check if
	//            flag >= 0
	//   opt == 2 means function allocates memory so check if returned
	//            NULL pointer 
	static int check_flag(void *flagvalue, char *funcname, int opt);
	
	realtype reltol;
	sdVector abstol;
	
	int linearMultistepMethod; // CV_ADAMS for non-stiff problems, CV_BDF for stiff problems
	int nonlinearSolverMethod; // CV_FUNCTIONAL for non-stiff problems, CV_NEWTON for stiff problems
	bool findRoots; // Specify whether or not to use the function g for rootfinding
	
	realtype t0; // initial time
	realtype tInt; // time reached by integrator
	sdVector y0;
	sdVector y;
	std::vector<int> rootsFound;
	unsigned int nRoots;
private:

	static int f(realtype t, N_Vector yIn, N_Vector ydotIn, void *f_data);
	static int g(realtype t, N_Vector yIn, realtype *gout, void *g_data);
	static int Jac(long int N, DenseMat Jin, realtype t,
	               N_Vector yIn, N_Vector fy, void *jac_data,
	               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
	
	sdODE* theODE;
	void *sundialsMem;
	int flag;
	int flagr;
	int nEq;
};


class sdDAE {
public:
	virtual int f(realtype t, sdVector& y, sdVector& ydot, sdVector& res)=0;
	virtual int g(realtype t, sdVector& y, sdVector& ydot, realtype* gOut) {return -1;}

	virtual int Jac(realtype t, sdVector& y, sdVector& ydot, sdVector& res,
			        realtype c_j, sdMatrix& J) {return -1;}

	virtual int JvProd(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn, 
					   sdVector& vIn, sdVector& JvIn, realtype c_j) {return -1;}

	virtual int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn, 
		                            sdVector& resIn, realtype c_j) {return -1;}

	virtual int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn,
									sdVector& resIn, sdVector& rhs, sdVector& outVec,
									realtype c_j, realtype delta) {return -1;}

	virtual ~sdDAE(void) {}
};

// wrapper class for the Sundials IDA solver
// for Differential-Algebraic Equations
class sundialsIDA
{
public:
	sundialsIDA(unsigned int n);
	~sundialsIDA(void);

	void initialize(void);
	int integrateToTime(realtype t);
	void setDAE(sdDAE* newDAE);
	int getRootInfo(); 
	void printStats();
	
	// Check function return value...
	//   opt == 0 means SUNDIALS function allocates memory so check if
	//            returned NULL pointer
	//   opt == 1 means SUNDIALS function returns a flag so check if
	//            flag >= 0
	//   opt == 2 means function allocates memory so check if returned
	//            NULL pointer 
	static int check_flag(void *flagvalue, char *funcname, int opt);
	
	realtype reltol;
	sdVector abstol;
	
	bool findRoots; // Specify whether or not to use the function g for rootfinding
	bool calcIC;

	realtype t0; // initial time
	realtype tInt; // time reached by integrator
	sdVector y;
	sdVector ydot;
	sdVector y0; // initial condition
	sdVector ydot0;

	sdVector componentId; // elements are 1.0 for differential variables, 0.0 for algebraic
	std::vector<int> rootsFound;
	unsigned int nRoots;
private:

	static int f(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn, void *f_data); // f(y) = res
	static int g(realtype t, N_Vector yIn, N_Vector ydotIn, realtype *gout, void *g_data);

	static int Jac(long int N, realtype t, N_Vector yIn, N_Vector ydotIn, 
		           N_Vector res, realtype c_j, void* jac_data, DenseMat Jin, 
				   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

	static int JvProd(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn,
					  N_Vector vIn, N_Vector JvIn, realtype c_j, void* jac_data,
					  N_Vector tmp1, N_Vector tmp2);

	static int preconditionerSetup(realtype t, N_Vector yIn, N_Vector ydotIn, 
		                           N_Vector resIn, realtype c_j, void* p_data, 
								   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
	
	static int preconditionerSolve(realtype t, N_Vector yIn, N_Vector ydotIn,
								   N_Vector resIn, N_Vector rhs, N_Vector outVec,
								   realtype c_j, realtype delta, void* p_data,
								   N_Vector tmp);
	sdDAE* theDAE;
	void *sundialsMem;
	int flag;
	int flagr;
	int nEq;
};


