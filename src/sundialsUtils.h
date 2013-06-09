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

//! wrapper class for Sundials "N_Vector"
class sdVector
{
public:
    //! Construct a new vector of length `n`.
    sdVector(unsigned int n);

    //! Wrap an existing `N_Vector`.
    sdVector(N_Vector other);

    sdVector(const sdVector& other);
    ~sdVector();

    template <class T>
    realtype& operator[](T i)
    {
        return NV_Ith_S(v, static_cast<long int>(i));
    }

    template <class T>
    realtype& operator[](T i) const
    {
        return NV_Ith_S(v, static_cast<long int>(i));
    }

    //! Return the underlying `N_Vector` object needed by Sundials functions.
    N_Vector& forSundials() { return v; }
    unsigned int length() const { return n; }
    size_t size() const { return n; }

private:
    N_Vector v;
    bool alloc;
    unsigned int n;
};

std::ostream& operator<<(std::ostream& os, const sdVector& v);

//! Wrapper class for Sundials `DenseMat` objects
class sdMatrix
{
public:
    //! Create a new matrix with `n` rows and `m` columns.
    sdMatrix(unsigned int n, unsigned int m);

    //! Create a wrapper for an existing `DensMat` object.
    sdMatrix(DenseMat other);

    sdMatrix();
    ~sdMatrix();

    template <class T>
    realtype& operator()(T i, T j)
    {
        return DENSE_ELEM(M, static_cast<long int>(i), static_cast<long int>(j));
    }

    template <class T>
    realtype& operator()(T i, T j) const
    {
        return DENSE_ELEM(M, static_cast<long int>(i), static_cast<long int>(j));
    }

    //! Return a pointer to the underlying data, needed by Sundials functions.
    realtype* forSundials() { return M->data; }

private:
    DenseMat M;
    bool alloc;
};

//! Wrapper class for Sundials `BandMat` objects
class sdBandMatrix
{
public:
    //! Create and wrap a new `BandMat`.
    //! @param N dimension of the matrix (square)
    //! @param bwUpper upper bandwidth of the matrix
    //! @param bwLower lower bandwidth of the matrix
    sdBandMatrix(long int N, long int bwUpper, long int bwLower);

    //! Create a wrapper for an existing `BandMat` object.
    sdBandMatrix(BandMat other);

    sdBandMatrix();
    ~sdBandMatrix();

    realtype& operator()(long int i, long int j)
    {
        return BAND_ELEM(M,i,j);
    }

    realtype& operator()(long int i, long int j) const
    {
        return BAND_ELEM(M,i,j);
    }

    //! Return a pointer to the underlying `BandMat` needed by Sundials functions
    BandMat& forSundials() { return M; }

    void print(const std::string& name="M") const;

private:
    BandMat M;
    bool alloc;
};

//! Abstract base class for an %ODE to be integrated by SundialsCvode.
class sdODE {
public:
    //! Evaluate the %ODE function \f[ y' = f(y,t) \f].
    //! @param t the current time
    //! @param y the current vector of state variables
    //! @param[out] ydot the computed time derivatives
    virtual int f(const realtype t, const sdVector& y, sdVector& ydot) = 0;

    //! Optional function for CVODE to find roots of.
    //! @param t the current time
    //! @param y the current vector of state variables
    //! @param[out] gOut the value of \f[ g(y,t) \f]
    virtual int g(realtype t, sdVector& y, realtype* gOut) { return 0; }

    //! Optional function used to compute the dense Jacobian.
    //! Used in conjuction with `CVDense`.
    //! @param t the current time
    //! @param y the current vector of state variables
    //! @param ydot the current vector of time derivatives
    //! @param[out] J The computed Jacobian matrix
    virtual int denseJacobian(const realtype t, const sdVector& y,
                              const sdVector& ydot, sdMatrix& J) { return -1; }

    //! Optional function used to compute a banded Jacobian.
    //! Used in conjuction with `CVBand`.
    //! @param t the current time
    //! @param y the current vector of state variables
    //! @param ydot the current vector of time derivatives
    //! @param[out] J The computed Jacobian matrix
    virtual int bandedJacobian(const realtype t, const sdVector& y,
                               const sdVector& ydot, sdBandMatrix& J) { return -1; }

    //! Called at the start of integration to perform any required initialization.
    virtual void initialize() {}
    virtual ~sdODE() {}
};

//! Wrapper class for the Sundials CVODE solver
/*!
 *  Use of this class has the following restrictions:
 *  - Configuration variables for CVODE may only be set before calling
 *    initialize(). The following *must* be initialized:
 *    - setBandwidth() - if using a banded Jacobian
 *    - setODE()
 *    - #t0
 *    - #y
 *    - #reltol
 *    - #abstol
 *    - #linearMultistepMethod
 *    - #nonlinearSolverMethod
 *  - The following configuration variables are optional:
 *    - #findRoots
 *    - #maxNumSteps
 *    - #minStep
 *  - Methods which call the CVODE solver, such as integrateToTime() may only
 *    be called after calling initialize().
 */
class SundialsCvode
{
public:
    //! Create a solver for a problem with `n` variables.
    SundialsCvode(unsigned int n);
    ~SundialsCvode();

    //! Set the *upper* and *lower* bandwidth of the banded Jacobian.
    void setBandwidth(int upper, int lower);

    //! Set the %ODE to be solved
    void setODE(sdODE* newODE);

    realtype t0; //!< initial time
    sdVector y; //!< initial or current state vector

    realtype reltol; //!< scalar relative tolerance
    sdVector abstol; //!< absolute tolerance for each solution component

    //! CV_ADAMS for non-stiff problems, CV_BDF for stiff problems
    int linearMultistepMethod;

    //! CV_FUNCTIONAL for non-stiff problems, CV_NEWTON for stiff problems
    int nonlinearSolverMethod;

    //! Specify whether or not to use the function g for rootfinding
    bool findRoots;

    //! Maximum number of internal timesteps taken in one call to integrateToTime
    int maxNumSteps;

    //! minimum step size
    realtype minStep;

    //! Initialize the CVODE solver object using the previously-specified
    //! options. Integration may be restarted by reassigning #y and #t0 and
    //! calling initialize again, provided that the problem size has not
    //! changed.
    void initialize();

    //! `true` if initialize() has already been called.
    bool initialized() const;

    //! Take as many steps as needed to reach `t` (up to #maxNumSteps)
    int integrateToTime(realtype t);

    //! Take one step towards `tf` without stepping past it.
    int integrateOneStep(realtype tf);

    //! Update #rootsFound to indicate any roots that have been found.
    int getRootInfo();

    //! Print solver statistics.
    void printStats();

    //! Returns the cumulative number of internal steps taken.
    long int getNumSteps();

    //! Returns the order of the integration method used during the last step
    int getLastOrder();

    //! Get the step size used for the last step
    realtype getLastStep();

    //! Check function return value for Sundials functions.
    //! - `opt == 0` means SUNDIALS function allocates memory so check if
    //!   returned NULL pointer
    //! - `opt == 1` means SUNDIALS function returns a flag so check if
    //!   `flag >= 0`
    static int check_flag(void* flagvalue, const char* funcname, int opt);

    realtype tInt; //!< time reached by integrator
    std::vector<int> rootsFound; //!< Root information populated by getRootInfo()
    unsigned int nRoots; //!< number of root functions to be called

    int errorCount; //!< Number of integration failures
    int errorStopCount; //!< Maximum permissible number of integration failures

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

//! Abstract base class for an system of differential-algebraic equations to
//! be integrated by SundialsIda.
class sdDAE {
public:
    //! Evaluate the residual of the DAE function \f[ f(y,y',t) = 0 \f].
    //! @param t the current time
    //! @param y the current vector of state variables
    //! @param ydot the current time derivatives
    //! @param[out] res the computed residual
    virtual int f(realtype t, sdVector& y, sdVector& ydot, sdVector& res) = 0;

    //! Optional function for IDA to find roots of.
    //! @param t the current time
    //! @param y the current vector of state variables
    //! @param ydot the current time derivatives
    //! @param[out] gOut the value of \f[ g(y,t) \f]
    virtual int g(realtype t, sdVector& y, sdVector& ydot, realtype* gOut) { return -1; }

    //! Jacobian routine. Compute J(t,y,y') = df/dy + c_j * df/d(y').
    virtual int Jac(realtype t, sdVector& y, sdVector& ydot, sdVector& res,
                    realtype c_j, sdMatrix& J) { return -1; }

    //! Compute the product J*y
    virtual int JvProd(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                       sdVector& vIn, sdVector& JvIn, realtype c_j) { return -1; }

    //! Set up the preconditioner matrix.
    //! Typically, this involves forming an approximate Jacobian and factorizing it.
    virtual int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                                    sdVector& resIn, realtype c_j) { return -1; }

    //! Solve the linear system P*outVec = rhs where P is the preconditioner.
    virtual int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn,
                                    sdVector& resIn, sdVector& rhs, sdVector& outVec,
                                    realtype c_j, realtype delta) { return -1; }

    virtual ~sdDAE() {}
};

//! Wrapper class for the Sundials IDA solver for differential-algebraic
//! equations.
/*!
 *  Use of this class has the following restrictions:
 *  - Configuration variables for IDA may only be set before calling
 *    initialize(). The following *must* be initialized:
 *    - setDAE()
 *    - #t0
 *    - #y
 *    - #ydot, if `#calcIC == false`
 *    - #reltol
 *    - #abstol
 *    - #constraints, if `#imposeConstraints == true`
 *    - #componentId, if `#calcIC == true`
 *  - The following configuration variables are optional:
 *    - #findRoots
 *    - #calcIC
 *    - #imposeConstraints
 *  - The following optional configuration methods may only be called *after*
 *    initialize():
 *    - setInitialStepSize()
 *    - setMaxStepSize()
 *    - disableErrorOutput()
 *  - Methods which call the IDA solver, such as integrateToTime() may only
 *    be called after calling initialize().
 */
class SundialsIda
{
public:
    //! Create a solver for a problem with `n` variables.
    SundialsIda(unsigned int n);
    ~SundialsIda();

    //! Initialize the IDA solver object using the previously-specified
    //! options. This function should only be called once.
    void initialize();

    //! Take as many steps as needed to reach `t`.
    int integrateToTime(realtype t);

    //! Take one step.
    int integrateOneStep();

    //! Set the DAE to be solved
    void setDAE(sdDAE* newDAE);

    //! Update #rootsFound to indicate any roots that have been found.
    int getRootInfo();

    //! Print solver statistics.
    void printStats(clock_t dt = 0);

    //! Check function return value...
    //! - opt == 0 means SUNDIALS function allocates memory so check if
    //!            returned NULL pointer
    //! - opt == 1 means SUNDIALS function returns a flag so check if
    //!            flag >= 0
    static int check_flag(void* flagvalue, const char* funcname, int opt);

    //! Get the step size used for the last step
    double getStepSize();

    //! Set the initial step size to attempt.
    void setInitialStepSize(double dt);

    //! Set the maximum allowable step size.
    void setMaxStepSize(double dt);

    //! Returns the order of the integration method used during the last step
    int getLastOrder();

    //! Disable printing error messages to stdout
    void disableErrorOutput();

    realtype reltol; //!< scalar relative tolerance
    sdVector abstol; //!< absolute tolerance for each solution component

    //! Specify whether or not to use the function SdDAE::g for rootfinding
    bool findRoots;

    //! `true` if IDA should compute the initial value for #ydot
    bool calcIC;

    //! `true` if constraints set in #constraints should be imposed
    bool imposeConstraints;

    realtype t0; // initial time
    realtype tInt; // time reached by integrator
    sdVector y; //!< initial or current state vector
    sdVector ydot;  //!< initial or current derivatives of state vector
    sdVector y0; //!< initial condition (output only, set only if using #calcIC)
    sdVector ydot0; //! initial derivatives (output only, set only if using #calcIC)

    //! elements are 1.0 for differential variables, 0.0 for algebraic. Needs
    //! to be set only if using the solver's IC finder (#calcIC == true)
    sdVector componentId;

    //! Used to impose constraints on solution variables.
    //! The following constraints may be specified:
    //!
    //!      0.0 - no constraint imposed
    //!      1.0 - y[i] >= 0
    //!     -1.0 - y[i] <= 0
    //!      2.0 - y[i] > 0
    //!     -2.0 - y[i] < 0
    sdVector constraints;

    std::vector<int> rootsFound; //!< Root information set by getRootInfo()
    unsigned int nRoots; //!< Number of active root functions

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
