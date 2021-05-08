#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "qssintegrator.h"
#include "quasi2d.h"
#include "callback.h"

class PerfTimer;
class CanteraGas;
class ScalarFunction;
class ConfigOptions;

//! Base class used to integrate the chemical source term at a single point.
class SourceSystem
{
public:
    SourceSystem();
    virtual ~SourceSystem() {}

    //! Set the initial condition for integrator
    //! @param tInitial  integrator start time
    //! @param uu        tangential velocity
    //! @param tt        temperature
    //! @param yy        vector of species mass fractions
    virtual void setState(double tInitial, double uu, double tt,
                          const dvec& yy) = 0;

    //! Take as many steps as needed to reach `tf`.
    //! `tf` is relative to `tInitial`.
    virtual int integrateToTime(double tf) = 0;

    //! Take one step toward `tf` without stepping past it.
    //! `tf` is relative to `tInitial`.
    virtual int integrateOneStep(double tf) = 0;

    //! Current integrator time, relative to `tInitial`.
    virtual double time() const = 0;

    //! Extract current internal integrator state into U, T, and Y
    virtual void unroll_y() = 0;

    virtual void setDebug(bool debug_) { debug = debug_; }

    //! Return a string indicating the number of internal timesteps taken
    virtual std::string getStats() = 0;

    //! Compute thermodynamic properties (density, enthalpies, heat capacities)
    //! from the current state variables.
    void updateThermo();

    //! Calculate the heat release rate associated with a possible external
    //! ignition source.
    double getQdotIgniter(double t);

    //! Set the CanteraGas object to use for thermodynamic and kinetic property
    //! calculations.
    void setGas(CanteraGas* _gas) { gas = _gas; }

    //! Resize internal arrays for a problem of the specified size (`nSpec+2`)
    virtual void initialize(size_t nSpec);

    //! Set integrator tolerances and other parameters
    virtual void setOptions(ConfigOptions& options_);

    void setTimers(PerfTimer* reactionRates, PerfTimer* thermo,
                   PerfTimer* jacobian);

    //! Set the index j and position x that are represented by this system.
    void setPosition(size_t j, double x);

    //! Set the function used to compute the strain rate as a function of time
    void setStrainFunction(ScalarFunction* f) { strainFunction = f; }

    //! Set the function used to compute the reaction rate multiplier
    void setRateMultiplierFunction(ScalarFunction* f) { rateMultiplierFunction = f; }

    //! Set the function used to compute the heat loss rate to the environment
    void setHeatLossFunction(IntegratorCallback* f) { heatLoss = f; }

    //! Set the density of the unburned mixture.
    //! This value appears in the source term of the momentum equation.
    void setRhou(double _rhou) { rhou = _rhou; }

    //! Set all the balanced splitting constants to zero.
    void resetSplitConstants() { splitConst.setZero(); }

    //! Assign the interpolators used for solving quasi-2D problems
    void setupQuasi2d(std::shared_ptr<BilinearInterpolator> vzInterp,
                      std::shared_ptr<BilinearInterpolator> TInterp);

    //! Write the current values of the state variables, formatted to be read by
    //! Python, to the specified stream. Call with `init=true` when first called
    //! to include initializers for the variables.
    virtual void writeState(std::ostream& out, bool init);

    virtual void writeJacobian(std::ostream& out) {};

    double U; //!< tangential velocity
    double T; //!< temperature
    dvec Y; //!< species mass fraction

    //! Extra constant term introduced by splitting
    dvec splitConst;

protected:
    bool debug;
    ConfigOptions* options;

    //! Cantera data
    CanteraGas* gas;

    //! Timer for time spent evaluating reaction rates
    PerfTimer* reactionRatesTimer;

    //! Timer for time spent evaluating thermodynamic properties
    PerfTimer* thermoTimer;

    //! Timer for time spent evaluating and factorizing the Jacobian
    PerfTimer* jacobianTimer;

    //! A class that provides the strain rate and its time derivative
    ScalarFunction* strainFunction;

    //! Provides a multiplier (optional) for the production terms
    ScalarFunction* rateMultiplierFunction;

    //! Heat loss rate
    IntegratorCallback* heatLoss;

    size_t nSpec; //!< number of species
    int j; //!< grid index for this system
    double x; //!< grid position for this system
    double rhou; //!< density of the unburned gas
    double qDot; //!< heat release rate per unit volume [W/m^3]
    double qLoss; //!< heat loss to the environment [W/m^3]

    // Physical properties
    double rho; //!< density [kg/m^3]
    double cp; //!< specific heat capacity (average) [J/kg*K]
    dvec cpSpec; //!< species specific heat capacity [J/mol*K]
    double Wmx; //!< mixture molecular weight [kg/mol]
    dvec W; //!< species molecular weights [kg/kmol]
    dvec hk; //!< species enthalpies [J/kmol]

    //! Flag set to 'true' when solving a quasi-2D problem with prescribed
    //! velocity and temperature fields.
    bool quasi2d;

    //! An interpolator for computing the axial (z) velocity when solving a
    //! quasi-2D problem
    std::shared_ptr<BilinearInterpolator> vzInterp;

    //! An interpolator for computing the temperature when solving a
    //! quasi-2D problem
    std::shared_ptr<BilinearInterpolator> TInterp;
};

//! Represents a system of equations used to integrate the (chemical) source
//! term at a single point using the CVODE integrator.
class SourceSystemCVODE : public SourceSystem, sdODE
{
public:
    SourceSystemCVODE() {}

    //! The ODE function: ydot = f(t,y)
    int f(const realtype t, const sdVector& y, sdVector& ydot);

    //! Calculate the Jacobian matrix: J = df/dy
    int denseJacobian(const realtype t, const sdVector& y,
                      const sdVector& ydot, sdMatrix& J);

    //! A simpler finite difference based Jacobian
    int fdJacobian(const realtype t, const sdVector& y,
                   const sdVector& ydot, sdMatrix& J);

    void setState(double tInitial, double uu, double tt, const dvec& yy);

    int integrateToTime(double tf);
    int integrateOneStep(double tf);
    double time() const;

    //! fill in the current state variables from `y`
    void unroll_y(const sdVector& y, double t);

    //! fill in the current state variables from the integrator state
    void unroll_y() {
        unroll_y(integrator->y, integrator->tInt);
    }

    //! fill in `y` with the values of the current state variables
    void roll_y(sdVector& y) const;

    //! fill in `ydot` with the values of the current time derivatives
    void roll_ydot(sdVector& ydot) const;

    std::string getStats();
    void initialize(size_t nSpec);
    void setOptions(ConfigOptions& options);

    virtual void writeState(std::ostream& out, bool init);

    //! Print the current Jacobian matrix ot the specified stream
    void writeJacobian(std::ostream& out);

    double dUdt; //!< time derivative of the tangential velocity
    double dTdt; //!< time derivative of the temperature
    dvec dYdt; //!< time derivative of the species mass fractions
    dvec wDot; //!< species net production rates [kmol/m^3*s]

private:
    std::unique_ptr<SundialsCvode> integrator;
};

//! This is the system representing the (chemical) source term at a point,
//! integrated with the QssIntegrator. The implementation of this class differs
//! from SourceSystemCVODE primarily in that the QSS integrator requires
//! separately calculating creation (Q) and destruction (D) rates of the state
//! variables.
class SourceSystemQSS : public SourceSystem, QssOde
{
public:
    SourceSystemQSS();

    //! The ODE function: ydot = f(t,y)
    void odefun(double t, const dvec& y, dvec& q, dvec& d,
                bool corrector=false);

    double time() const { return integrator.tn; }

    //! Assign the current state variables from `y`.
    //! The current value of the temperature is not updated during corrector
    //! iterations.
    void unroll_y(const dvec& y, bool corrector=false);

    //! Assign the state variables using the current integrator state.
    void unroll_y() { unroll_y(integrator.y); }

    //! fill in `y` with current state variables.
    void roll_y(dvec& y) const;

    //! fill in `q` and `d` with current creation and destruction rates for
    //! each component. The net time derivative is `ydot = q - d`.
    void roll_ydot(dvec& q, dvec& d) const;

    void initialize(size_t nSpec);
    void setOptions(ConfigOptions& options);

    void setState(double tStart, double uu, double tt, const dvec& yy);
    int integrateToTime(double tf) { return integrator.integrateToTime(tf); }
    int integrateOneStep(double tf) { return integrator.integrateOneStep(tf); }

    virtual std::string getStats();

    // creation and destruction rates of state variables
    double dUdtQ; //!< tangential velocity "creation" rate
    double dUdtD; //!< tangential velocity "destruction" rate
    double dTdtQ; //!< temperature "creation" rate
    double dTdtD; //!< temperature "destruction" rate
    dvec dYdtQ; //!< species mass fraction creation rate
    dvec dYdtD; //!< species mass fraction destruction rate

    double tCall; //!< the last time at which odefun was called
    dvec wDotQ, wDotD; //!< species production / destruction rates [kmol/m^3*s]

private:
    QssIntegrator integrator;
};
