#pragma once

#include "mathUtils.h"

//! Abstract base class for an %ODE initial value problem that can be
//! integrated by an Integrator.
class ODE
{
public:
    virtual ~ODE() {}
    //! %ODE function defined as `ydot = f(t,y)`
    virtual void f(const double t, const dvec& y, dvec& ydot) = 0;
};

//! A system of linear ODEs represented by a tridiagonal matrix.
//! The ODE may be written as an inhomogeneous linear linear system:
//!     \f[ \dot{y} = Ay + k \f]
class TridiagonalODE
{
public:
    TridiagonalODE() {}
    virtual ~TridiagonalODE() {}

    //! Provide the matrix associated with the ODE to the integrator.
    //! @param[out] a elements of the first subdiagonal of *A*
    //! @param[out] b elements of the main diagonal of *A*
    //! @param[out] c elements of the first superdiagonal of *A*
    virtual void get_A(dvec& a, dvec& b, dvec& c) = 0;

    //! Provides the constant term *k* to the integrator
    virtual void get_k(dvec& k) = 0;

    //! Set the number of points in the domain.
    virtual void resize(size_t N) {}

    virtual void initialize() {}
};

//! Abstract base class for an %ODE integrator that can integrate instances of
//! class ODE.
class Integrator
{
public:
    Integrator();
    virtual ~Integrator() {}

    //! Set the initial condition for the problem. Must be called before
    //! starting integration.
    virtual void set_y0(const dvec& y0);

    //! Set up parameters and problem-dependent data structures for the
    //! solver. Must be called before starting integration.
    //! @param t0 Start time for the integration
    //! @param h initial step size
    virtual void initialize(const double t0, const double h);

    //! Get the last step size used
    double get_h() const;

    //! Get the current time reached by the integrator
    double get_t() const;

    //! Get the current state vector
    virtual const dvec& get_y() const;

    //! Get the time derivative of the current state vector
    virtual const dvec& get_ydot() = 0;

    //! Take a single step using the current step size.
    virtual void step() = 0;

    //! Take as many steps as necessary to reach *tEnd* without stepping past
    //! it.
    virtual void integrateToTime(double tEnd);

    dvec y; //!< solution vector
    dvec ydot; //!< derivative of state vector
    double t; //!< current time

protected:
    double h; //!< timestep
    size_t N; //!< Dimension of y
};

//! Integrates an %ODE defined as `ydot = f(t,y)` using the explicit Euler
//! method. This class exists mostly as a demonstration of how to implement
//! Integrator.
class ExplicitIntegrator : public Integrator
{
public:
    ExplicitIntegrator(ODE& ode);
    ~ExplicitIntegrator() {}

    void set_y0(const dvec& y0);
    const dvec& get_ydot();

    void step();

private:
    ODE& myODE;
    dvec ydot;
};

//! Integrator using 1st and 2nd order BDF for tridiagonal systems, with
//! matrix inversions done using the Thomas algorithm.
class TridiagonalIntegrator : public Integrator
{
public:
    TridiagonalIntegrator(TridiagonalODE& ode);
    virtual ~TridiagonalIntegrator() {}

    void set_y0(const dvec& y0);
    void resize(size_t N_in);
    void initialize(const double t0, const double h);

    const dvec& get_ydot();

    void step();

private:
    TridiagonalODE& myODE;

    dvec a; //!< subdiagonal components of right-hand-side
    dvec b; //!< diagonal components of right-hand-side
    dvec c; //!< superdiagonal components of right-hand-side
    dvec k; //!< constant contribution to `y'`
    dvec lu_b, lu_c, lu_d; // matrix elements after factorization
    dvec invDenom_; // internal
    unsigned int stepCount; //!< the number of steps taken since last initialization

    int N; //!< The system size
    dvec yprev; //!< previous solution
};
