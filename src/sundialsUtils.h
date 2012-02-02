#pragma once

#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_band.h>

#include <sundials/sundials_types.h>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>
#include <ida/ida.h>
//#include <ida/ida_dense.h>
#include <ida/ida_spbcgs.h>

#include <iostream>
#include <vector>
#include <time.h>

typedef DlsMat DenseMat;
typedef DlsMat BandMat;

// wrapper class for Sundials "N_Vector"
class sdVector
{
public:
    sdVector(unsigned int n);
    sdVector(N_Vector other);
    sdVector(const sdVector& other);
    ~sdVector();

    realtype& operator[](unsigned int i);
    realtype& operator[](unsigned int i) const;
    N_Vector& forSundials() { return v; }
    unsigned int length() const { return n; }
    size_t size() const { return n; }

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
    sdMatrix();
    ~sdMatrix();
    realtype& operator()(unsigned int i, unsigned int j);
    realtype& operator()(unsigned int i, unsigned int j) const;
    realtype* forSundials() { return M->data; }

private:
    DenseMat M;
    bool alloc;
};

class sdBandMatrix
{
public:
    sdBandMatrix(long int N, long int bwUpper, long int bwLower);
    sdBandMatrix(BandMat other);
    sdBandMatrix();
    ~sdBandMatrix();
    realtype& operator()(long int i, long int j);
    realtype& operator()(long int i, long int j) const;
    BandMat& forSundials() { return M; }

private:
    BandMat M;
    bool alloc;
};

class sdODE {
public:
    virtual int f(const realtype t, const sdVector& y, sdVector& ydot) = 0;
    virtual int g(realtype t, sdVector& y, realtype* gOut) { return 0; }
    virtual int denseJacobian(const realtype t, const sdVector& y,
                              const sdVector& ydot, sdMatrix& J) { return -1; }
    virtual int bandedJacobian(const realtype t, const sdVector& y,
                               const sdVector& ydot, sdBandMatrix& J) { return -1; }
    virtual void initialize() {}
    virtual ~sdODE() {}
};

// wrapper class for the Sundials CVODE solver
class sundialsCVODE
{
public:
    sundialsCVODE(unsigned int n);
    ~sundialsCVODE();

    // *** INITIALIZATION ***
    // Call each of these functions and assign each of these
    // variables before calling initialize().

    void setBandwidth(int upper, int lower); // Required if using a banded Jacobian
    void setODE(sdODE* newODE); // Required

    realtype t0; // initial time
    sdVector y; // initial or current state vector

    realtype reltol;
    sdVector abstol;
    int linearMultistepMethod; // CV_ADAMS for non-stiff problems, CV_BDF for stiff problems
    int nonlinearSolverMethod; // CV_FUNCTIONAL for non-stiff problems, CV_NEWTON for stiff problems
    bool findRoots; // Specify whether or not to use the function g for rootfinding
    int maxNumSteps; // Maximum number of internal timesteps taken in one call to integrateToTime
    realtype minStep; // minimum step size [s]

    // Call initialize() before using integrateToTime to integrate.
    // Integration may be restarted by reassigning y and t0 and calling initialize again.
    void initialize();
    bool initialized() const;

    int integrateToTime(realtype t); // Take as many steps as needed to reach t (up to maxNumSteps)
    int integrateOneStep(realtype tf);  // Take one step towards tf without stepping past it

    // optional output functions
    int getRootInfo();
    void printStats();
    long int getNumSteps(); // cumulative number of internal steps
    int getLastOrder(); // order used during the last step
    realtype getLastStep(); // step size used for the last step

    // Check function return value...
    //   opt == 0 means SUNDIALS function allocates memory so check if
    //            returned NULL pointer
    //   opt == 1 means SUNDIALS function returns a flag so check if
    //            flag >= 0
    //   opt == 2 means function allocates memory so check if returned
    //            NULL pointer
    static int check_flag(void* flagvalue, const char* funcname, int opt);

    realtype tInt; // time reached by integrator
    std::vector<int> rootsFound;
    unsigned int nRoots;

    int errorCount;
    int errorStopCount;

private:
    static int f(realtype t, N_Vector yIn, N_Vector ydotIn, void* f_data);
    static int g(realtype t, N_Vector yIn, realtype *gout, void* g_data);
    static int denseJac(int N, realtype t, N_Vector yIn,
                        N_Vector fy, DenseMat Jin, void* jac_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int bandJac(int N, int mupper, int mLower, realtype t,
                       N_Vector y, N_Vector fy, DlsMat Jac, void* user_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    sdODE* theODE;
    void* sundialsMem;
    int nEq;

    int bandwidth_upper;
    int bandwidth_lower;

    bool _initialized;
};


class sdDAE {
public:
    virtual int f(realtype t, sdVector& y, sdVector& ydot, sdVector& res) = 0;
    virtual int g(realtype t, sdVector& y, sdVector& ydot, realtype* gOut) { return -1; }

    virtual int Jac(realtype t, sdVector& y, sdVector& ydot, sdVector& res,
                    realtype c_j, sdMatrix& J) { return -1; }

    virtual int JvProd(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                       sdVector& vIn, sdVector& JvIn, realtype c_j) { return -1; }

    virtual int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                                    sdVector& resIn, realtype c_j) { return -1; }

    virtual int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn,
                                    sdVector& resIn, sdVector& rhs, sdVector& outVec,
                                    realtype c_j, realtype delta) { return -1; }

    virtual ~sdDAE() {}
};

// wrapper class for the Sundials IDA solver
// for Differential-Algebraic Equations
class sundialsIDA
{
public:
    sundialsIDA(unsigned int n);
    ~sundialsIDA();

    void initialize();
    int integrateToTime(realtype t);
    int integrateOneStep();
    void setDAE(sdDAE* newDAE);
    int getRootInfo();
    void printStats(clock_t dt = 0);

    // Check function return value...
    //   opt == 0 means SUNDIALS function allocates memory so check if
    //            returned NULL pointer
    //   opt == 1 means SUNDIALS function returns a flag so check if
    //            flag >= 0
    //   opt == 2 means function allocates memory so check if returned
    //            NULL pointer
    static int check_flag(void* flagvalue, const char* funcname, int opt);

    double getStepSize();
    void setInitialStepSize(double dt);
    void setMaxStepSize(double dt);
    int getLastOrder();
    void disableErrorOutput();

    realtype reltol;
    sdVector abstol;

    bool findRoots; // Specify whether or not to use the function g for rootfinding
    bool calcIC;
    bool imposeConstraints;

    realtype t0; // initial time
    realtype tInt; // time reached by integrator
    sdVector y;
    sdVector ydot;
    sdVector y0; // initial condition (output only)
    sdVector ydot0;

    // elements are 1.0 for differential variables, 0.0 for algebraic
    // Needs to be set only if using the solver's IC finder (calcIC==true)
    sdVector componentId;
    sdVector constraints;

    std::vector<int> rootsFound;
    unsigned int nRoots;

private:
    static int f(realtype t, N_Vector yIn, N_Vector ydotIn, N_Vector resIn, void* f_data); // f(y) = res
    static int g(realtype t, N_Vector yIn, N_Vector ydotIn, realtype *gout, void* g_data);

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
    void* sundialsMem;
    int nEq;
};
