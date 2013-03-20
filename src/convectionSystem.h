#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "grid.h"
#include "perfTimer.h"
#include "quasi2d.h"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class CanteraGas;
class ConvectionTermWrapper;

//! System representing the coupled convection equations for U, T, and Wmx.
/*!
 *  Solves equations for the tangential velocity gradient `U`, the temperature
 *  `T` and the mixture molecular weight `Wmx` which are coupled with the
 *  continuity equation.
 */
class ConvectionSystemUTW : public sdODE, public GridBased
{
public:
    ConvectionSystemUTW();

    //! The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);

    //! Use the contents of *y* to fill in the current state variables
    void unroll_y(const sdVector& y);

    //! Fill in *y* with the current values of the state variables
    void roll_y(sdVector& y) const;

    //! Fill in *ydot* with the current time derivatives
    void roll_ydot(sdVector& ydot) const;

    //! Set the size of the domain to *nPoints*
    void resize(const size_t nPoints);

    //! Set the term in each equation representing contributions from
    //! diffusion or production terms to zero.
    void resetSplitConstants();

    //! Determine the index of the grid point at which to impose the boundary
    //! condition for the continuity equation, #jContBC, based on the boundary
    //! condition type and the current mass flux #V.
    void updateContinuityBoundaryCondition(const dvec& qdot,
                                           ContinuityBoundaryCondition::BC newBC);

    dvec U; //!< Tangential velocity gradient [1/s]. State variable.
    dvec dUdt; //!< Time derivative of `U`, [1/s^2].
    dvec T; //!< Temperature [K]. State variable.
    dvec dTdt; //!< Time derivative of `T` [K/s].
    dvec Wmx; //!< Mixture molecular weight [kg/kmol]. State variable.
    dvec dWdt; //!< Time derivative of `Wmx`.

    double Tleft; //!< Temperature left boundary value
    double Wleft; //!< mixture molecular weight left boundary value
    double rVzero; //!< mass flux boundary value at j=0

    dvec drhodt; //!< Time derivative of the density [kg/m^3*s]

    dvec splitConstT; //!< Constant in `T` equation introduced by splitting method.
    dvec splitConstW; //!< Constant in `Wmx` equation introduced by splitting method.
    dvec splitConstU; //!< Constant in `U` equation introduced by splitting method.

    //! Cantera data
    CanteraGas* gas;

    // Auxiliary variables
    dvec V; //!< mass flux [kg/m^2*s]
    dvec rV; //!< (radial) mass flux (r*V) [kg/m^2*s or kg/m*rad*s]
    dvec rho; //!< mixture density [kg/m^3]

    dvec dUdx; //!< Gradient of `U` normal to the flame
    dvec dTdx; //!< Gradient of `T` normal to the flame
    dvec dWdx; //!< Gradient of `Wmx` normal to the flame

    //! The method used for determining the continuity boundary condition.
    ContinuityBoundaryCondition::BC continuityBC;
    size_t jContBC; //!< The point at which the continuity equation BC is applied.
    double xVzero; //!< Location of the stagnation point (if using fixedZero BC)

private:
    void V2rV(); //!< compute #rV from #V
    void rV2V(); //!< compute #V from #rV

    size_t nVars; //!< Number of state variables at each grid point (`== 3`)
};

typedef std::map<double, dvec> vecInterpolator;

//! System representing the convection equation for a single species.
//! The species is convected by a prescribed velocity field \f$v(x)\f$.
class ConvectionSystemY : public sdODE, public GridBased
{
public:
    ConvectionSystemY() : quasi2d(false) {}

    //! The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);

    //! Set the size of the domain to *nPoints*
    void resize(const size_t nPoints);

    //! Set the term representing contributions from diffusion or production
    //! terms to zero.
    void resetSplitConstants();

     //! Mass fraction to the left of the domain. Used only in conjuction with
     //! BoundaryCondition::ControlVolume.
    double Yleft;

    int k; //!< Species index. Used for debugging purposes only.
    dvec splitConst; //!< constant term introduced by splitting

    //! Axial (normal) velocity [m/s] as a function of time. This velocity
    //! data is computed by ConvectionSystemUTW. The keys are the times [s]
    //! for which the corresponding velocity field was computed.
    boost::shared_ptr<vecInterpolator> vInterp;

    //! Interpolator for the *z* velocity in the quasi-2d problem.
    boost::shared_ptr<BilinearInterpolator> vzInterp;

    //! Interpolator for the *r* velocity in the quasi-2d problem.
    boost::shared_ptr<BilinearInterpolator> vrInterp;

    //! `true` when solving the quasi-2d problem; `false` otherwise.
    bool quasi2d;

private:
    //! Compute the value of #v at time `t` by linearly interpolating between
    //! the adjacent velocity fields stored in #vInterp.
    void update_v(const double t);

    //! The velocity normal to the flame [m/s].
    dvec v;
};


class ConvectionSystemSplit : public GridBased
{
    // System which combines a ConvectionSystemUTW objects and several
    // ConvectionSystemY objects that together represent the complete
    // convection term for all components.
public:
    ConvectionSystemSplit();

    void setGrid(const OneDimGrid& grid);
    void setTolerances(const ConfigOptions& options);
    void setGas(CanteraGas& gas);
    void resize(const size_t nPoints, const size_t nSpec, dmatrix& state);
    void setState(double tInitial);
    void setLeftBC(const double Tleft, const dvec& Yleft);
    void set_rVzero(const double rVzero);
    void evaluate(); // evaluate time derivatives and mass flux at the current state

    // Time derivatives of species and temperature from the other split terms are needed
    // to correctly compute the density derivative appearing in the continuity equation
    void setDensityDerivative(const dvec& drhodt);

    // Constants introduced by the splitting method
    void setSplitConstants(const dmatrix& splitConst);

    void resetSplitConstants();

    void integrateToTime(const double tf);
    void unroll_y(); // convert the solver's solution vectors to the full U, Y, and T
    int getNumSteps();
    void setupQuasi2D(boost::shared_ptr<BilinearInterpolator>& vzInterp,
                      boost::shared_ptr<BilinearInterpolator>& vrInterp);

    VecMap U;
    VecMap T;
    MatrixMap Y;
    dvec Wmx;

    // Time derivatives and mass flux are updated by evaluate()
    dvec V;
    dvec dUdt;
    dvec dTdt;
    dvec dWdt;
    dmatrix dYdt;

    ConvectionSystemUTW utwSystem;

    boost::shared_ptr<vecInterpolator> vInterp;

    PerfTimer utwTimer, speciesTimer;

private:
    // set parameters of a new species solver
    void configureSolver(SundialsCvode& solver, const size_t k);

    // CVODE integration tolerances
    double reltol; // relative tolerance
    double abstolU; // velocity absolute tolerance
    double abstolT; // temperature absolute tolerance
    double abstolW; // molecular weight absolute tolerance
    double abstolY; // mass fraction absolute tolerance

    boost::shared_ptr<SundialsCvode> utwSolver;
    boost::ptr_vector<ConvectionSystemY> speciesSystems;
    boost::ptr_vector<SundialsCvode> speciesSolvers;

    dvec Yleft;
    dvec W;

    size_t nSpec;
    size_t nVars;

    CanteraGas* gas;

    bool quasi2d;
    friend class ConvectionTermWrapper;
};

class ConvectionTermWrapper {
public:
    ConvectionTermWrapper(ConvectionSystemSplit* parent_, double t_)
        : parent(parent_)
        , t(t_)
        {}
    void operator()(const tbb::blocked_range<size_t>& r) const;
private:
    ConvectionSystemSplit* parent;
    double t;
};
