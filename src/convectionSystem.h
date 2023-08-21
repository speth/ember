#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "grid.h"
#include "perfTimer.h"
#include "quasi2d.h"
#include "scalarFunction.h"

#include <boost/ptr_container/ptr_vector.hpp>

class CanteraGas;

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

    //! An object that provides the strain rate and its time derivative
    ScalarFunction* strainFunction;

    //! Set the function used to compute the strain rate as a function of time
    void setStrainFunction(ScalarFunction* f) { strainFunction = f; }

    //! Set the density of the unburned mixture.
    //! This value appears in the source term of the momentum equation.
    double rhou; //!< density of the unburned gas
    void setRhou(double _rhou) { rhou = _rhou; }

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

    size_t k; //!< Species index. Used for debugging purposes only.
    dvec splitConst; //!< constant term introduced by splitting

    //! Axial (normal) velocity [m/s] as a function of time. This velocity
    //! data is computed by ConvectionSystemUTW. The keys are the times [s]
    //! for which the corresponding velocity field was computed.
    std::shared_ptr<vecInterpolator> vInterp;

    //! Interpolator for the *z* velocity in the quasi-2d problem.
    std::shared_ptr<BilinearInterpolator> vzInterp;

    //! Interpolator for the *r* velocity in the quasi-2d problem.
    std::shared_ptr<BilinearInterpolator> vrInterp;

    //! `true` when solving the quasi-2d problem; `false` otherwise.
    bool quasi2d;

private:
    //! Compute the value of #v at time `t` by linearly interpolating between
    //! the adjacent velocity fields stored in #vInterp.
    void update_v(const double t);

    //! The velocity normal to the flame [m/s].
    dvec v;
};


//! Composite system representing the Convection term for all components.
/*!
 * System which combines a ConvectionSystemUTW object and several
 * ConvectionSystemY objects that together represent the complete convection
 * term for all components.
 */
class ConvectionSystemSplit : public GridBased
{
public:
    ConvectionSystemSplit();

    void setGrid(const OneDimGrid& grid);

    //! Set tolerances for the CVODE solvers.
    void setTolerances(const ConfigOptions& options);

    //! Set the Cantera object used for property evaluations.
    void setGas(CanteraGas& gas);

    //! Set the problem size and provide values for the current state variables.
    void resize(const size_t nPoints, const size_t nSpec, dmatrix& state);

    //! Set the state of the internal solvers from the state of the composite
    //! solver.
    void setState(double tInitial);

    //! Set boundary conditions for temperature and mass fractions left of the
    //! domain when using BoundaryCondition::ControlVolume or
    //! BoundaryCondition::WallFlux.
    void setLeftBC(const double Tleft, const dvec& Yleft);

    //! Set the mass flux boundary value at j=0.
    void set_rVzero(const double rVzero);

    //! Evaluate time derivatives for each state variable and the normal mass
    //! flux at the current state.
    void evaluate();

    //! Set the time derivative of the density. Time derivatives of species
    //! and temperature from the other split terms are needed to correctly
    //! compute the density derivative appearing in the continuity equation.
    void setDensityDerivative(const dvec& drhodt);

    //! @see ConvectionSystemUTW::updateContinuityBoundaryCondition
    void updateContinuityBoundaryCondition(const dvec& qdot,
                                           ContinuityBoundaryCondition::BC newBC);

    //! Set the constants introduced by the splitting method.
    void setSplitConstants(const dmatrix& splitConst);

    //! Set the constants introduced by the splitting method to zero.
    void resetSplitConstants();

    //! Integrate the the split convection equations to *tf*.
    void integrateToTime(const double tf);

    //! Integrate the species terms in the range [k1, k2)
    void integrateSpeciesTerms(size_t k1, size_t k2);


    //! convert the solver's solution vectors to the full *U*, *Y*, and *T*.
    void unroll_y();

    //! Compute the total number of timesteps taken by the solvers for the
    //! individual components.
    int getNumSteps();

    //! Set the velocity field data used in the Quasi2D case.
    void setupQuasi2D(std::shared_ptr<BilinearInterpolator>& vzInterp,
                      std::shared_ptr<BilinearInterpolator>& vrInterp);

    VecMap U; //!< normalized tangential velocity (u*a/u_inf) [1/s]
    VecMap T; //!< temperature [K]
    MatrixMap Y; //!< species mass fractions, Y(k,j) [-]
    dvec Wmx; //!< Mixture molecular weight [kg/kmol]

    // Time derivatives and mass flux are updated by evaluate()
    dvec V; //!< mass flux normal to the flame [kg/m^2*s]
    dvec dUdt; //!< Time derivative of #U [1/s^2]
    dvec dTdt; //!< Time derivative of #T [K/s]
    dvec dWdt; //!< Time derivative of #Wmx [kg/kmol*s]
    dmatrix dYdt; //!< Time derivative of #Y [1/s]

    //! System used to solve for #U, #T, and #Wmx
    ConvectionSystemUTW utwSystem;

    //! Interpolation data for #V, used for the species solvers
    std::shared_ptr<vecInterpolator> vInterp;

    PerfTimer utwTimer, speciesTimer;

private:
    // set parameters of a new species solver
    void configureSolver(SundialsCvode& solver, const size_t k);

    //! CVODE integration tolerances
    double reltol; //!< relative integrator tolerance
    double abstolU; //!< velocity absolute tolerance
    double abstolT; //!< temperature absolute tolerance
    double abstolW; //!< molecular weight absolute tolerance
    double abstolY; //!< mass fraction absolute tolerance

    //! Solver for #utwSystem
    std::shared_ptr<SundialsCvode> utwSolver;

    //! Systems used to solve the convection term for each species.
    boost::ptr_vector<ConvectionSystemY> speciesSystems;

    //! Solvers for system in #speciesSystems
    boost::ptr_vector<SundialsCvode> speciesSolvers;

    //! Mass fraction left of the domain. Used with
    //! BoundaryCondition::ControlVolume and BoundaryCondition::WallFlux.
    dvec Yleft;

    dvec W; //!< Molecular weight of each species [kg/kmol]

    size_t nSpec; //!< Number of species
    size_t nVars; //!< Number of state variables in the UTW system (`==3`)

    CanteraGas* gas; //!< Cantera object used for computing #Wmx

    //! `true` when solving the quasi-2d problem; `false` otherwise.
    bool quasi2d;

    double tStageStop; //!< end time of the current integration stage

    SundialsContext sunContext;
};
