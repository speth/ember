#pragma once

#undef NO_ERROR // fix for interaction between reaction_defs.h and winerror.h

#include "cantera/IdealGasMix.h"    // defines class IdealGasMix

#include "cantera/equilibrium.h"    // chemical equilibrium
#include "cantera/thermo.h"
#include "cantera/transport.h"      // transport properties
#include "cantera/kinetics.h"

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/base/ctml.h"

#include "mathUtils.h"

class ConfigOptions;

//! Mixture-averaged transport properties based on major species composition.
//! For a mixture containing N species, the standard formulas for computing the
//! mixture-averaged viscosity and diffusion coefficients are \f$ O(N^2) \f$.
//!
//! Species present in small concentrations contribute very little to the
//! mixture viscosity and the diffusion coefficients of other species, so they
//! can be excluded from the formulas without introducing significant errors.
//! If there are *M* species with mole fractions above the specified threshold,
//! (*M* < *N*) then excluding the contribution of the minor species reduces
//! the computational cost of evaluating the viscosity to \f$ O(M^2) \f$ and
//! the cost of evaluating the *N* diffusion coefficients to \f$ O(M N) \f$.
class ApproxMixTransport : public Cantera::MixTransport
{
public:
    ApproxMixTransport(Cantera::ThermoPhase& thermo,
                       Cantera::TransportFactory& factory);

    //! Set the mole fraction threshold below which species are not included
    //! in transport property calculations.
    void setThreshold(double threshold);

    double viscosity();
    void getMixDiffCoeffs(double* const d);
    void getMixDiffCoeffsMass(double* const d);
    void getMixDiffCoeffsMole(double* const d);

private:
    void updateViscosity_T();
    void updateDiff_T();
    void update_C();

    double _threshold;
    vector<size_t> _kMajor; //!< indices of the species where X[k] >= threshold
};


class InterpKinetics : public Cantera::GasKinetics
{
public:
    InterpKinetics(Cantera::ThermoPhase* phase);

    virtual void update_rates_T();
    void rebuildInterpData(size_t nTemps, double Tmin, double Tmax);

private:
    size_t m_ntemps; //!< number of temperatures over which to interpolate
    double m_tmin; //!< minimum temperature
    double m_tmax; //!< maximum temperature

    dmatrix m_rfn_const, m_rfn_slope;
    dmatrix m_rfn_low_const, m_rfn_low_slope;
    dmatrix m_rfn_high_const, m_rfn_high_slope;
    dmatrix m_falloff_work_const, m_falloff_work_slope;
    dmatrix m_rkcn_const, m_rkcn_slope;
};

//! A set of Cantera objects needed for calculating thermodynamic properties,
//! transport properties, and kinetic rates for a constant-pressure mixture.
class CanteraGas
{
public:
    CanteraGas();
    ~CanteraGas();

    double pressure; //!< thermodynamic pressure [Pa]
    size_t nSpec; //!< number of species

    void setOptions(const ConfigOptions& options);
    void initialize();
    bool initialized() const { return isInitialized; }

    //! Calculate the mole fractions of the reactant mixture from the
    //! compositions of the fuel and oxidizer mixtures and the equivalence ratio.
    dvec calculateReactantMixture(const std::string& fuel,
                                  const std::string& oxidizer,
                                  double equivalenceRatio);

    void setStateMass(const dvec& Y, const double T);
    void setStateMass(const double* Y, const double T);

    void setStateMole(const dvec& X, const double T);
    void setStateMole(const double* X, const double T);

    void getMoleFractions(dvec& X) const;
    void getMoleFractions(double* X) const;

    void getMassFractions(dvec& Y) const;
    void getMassFractions(double* Y) const;

    double getDensity() const;
    double getMixtureMolecularWeight() const;

    void getMolecularWeights(dvec& W) const;
    void getMolecularWeights(double* W) const;

    double getViscosity() const;
    double getThermalConductivity() const;

    //! Get diffusion coefficients for calculating mass diffusion velocities
    //! with respect to mole fraction gradients.
    //! `v[k] = -D[k] / X[k] * grad(X[k])`
    void getDiffusionCoefficientsMole(dvec& Dkm) const;
    void getDiffusionCoefficientsMole(double* Dkm) const;

    //! Get product of density and diffusion coefficients for calculating
    //! diffusive mass fluxes with respect to mole fraction gradients.
    //! `j[k] = - rhoD[k] / X[k] * grad(X[k])`
    void getWeightedDiffusionCoefficientsMole(dvec& rhoD) const;
    void getWeightedDiffusionCoefficientsMole(double* rhoD) const;

    //! Get product of density and diffusion coefficients for calculating
    //! diffusive mass fluxes with respect to mass fraction gradients.
    //! `j[k] = - rho * D[k] * grad(Y[k])`
    void getWeightedDiffusionCoefficientsMass(double* rhoD);
    void getWeightedDiffusionCoefficientsMass(dvec& rhoD);

    //! Get thermal diffusion coefficients.
    //! `j[k] = - Dt[k] / (T * Y[k]) * grad(T)`
    void getThermalDiffusionCoefficients(dvec& Dkt) const;
    void getThermalDiffusionCoefficients(double* Dkt) const;

    double getSpecificHeatCapacity() const;
    void getSpecificHeatCapacities(dvec& cpSpec) const;
    void getSpecificHeatCapacities(double* cpSpec) const;
    void getEnthalpies(dvec& hk) const;
    void getEnthalpies(double* hk) const;

    //! Get net molar reaction rates for each species
    void getReactionRates(dvec& wDot) const;
    void getReactionRates(double* wDot) const;

    void getCreationRates(dvec& wDot) const;
    void getCreationRates(double* wDot) const;

    void getDestructionRates(dvec& wDot) const;
    void getDestructionRates(double* wDot) const;

    Cantera::IdealGasPhase thermo;

private:
    std::string mechanismFile;
    std::string phaseID;
    std::string transportModel;
    std::string kineticsModel;
    double transportThreshold;
    bool isInitialized;

    Cantera::XML_Node* rootXmlNode;
    Cantera::XML_Node* phaseXmlNode;

    Cantera::GasKinetics* kinetics;
    Cantera::GasTransport* transport;

    dmatrix Dbin; //!< binary diffusion coefficients for species k
};
