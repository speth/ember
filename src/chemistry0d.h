#pragma once

#undef NO_ERROR // fix for interaction between reaction_defs.h and winerror.h

#include "config.h"
#include "mathUtils.h"

#include "cantera/transport/MixTransport.h"
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/UnityLewisTransport.h"
#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/base/Solution.h"

using namespace Cantera;

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
    ApproxMixTransport();

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


class InterpKinetics : public Cantera::BulkKinetics
{
public:
    InterpKinetics();
};

class MultiArrheniusInterp : public Cantera::MultiRateBase
{
public:
    std::string type() override {
        return m_base.type();
    }

    void add(size_t rxn_index, ReactionRate& rate) override {
        m_indices.emplace_back(rxn_index, m_nReactions++);
        m_base.add(rxn_index, rate);
    }

    bool replace(size_t rxn_index, ReactionRate& rate) override {
        throw DebugException("MultiArrheniusInterp::replace: operation not supported.");
    }

    void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        m_base.resize(nSpecies, nReactions, nPhases);
    }

    void getRateConstants(double* kf) override {
        for (const auto& [jGlobal, jLocal] : m_indices) {
            kf[jGlobal] = m_kf_const(jLocal, m_T_index) + m_kf_slope(jLocal, m_T_index) * m_dT;
        }
    }

    void processRateConstants_ddT(double* rop, const double* kf, double deltaT) override {
        m_base.processRateConstants_ddT(rop, kf, deltaT);
    }

    void processRateConstants_ddP(double* rop, const double* kf, double deltaP) override {
        m_base.processRateConstants_ddP(rop, kf, deltaP);
    }

    void processRateConstants_ddM(double* rop, const double* kf, double deltaM,
                                  bool overwrite=true) override
    {
        m_base.processRateConstants_ddM(rop, kf, deltaM, overwrite);
    }

    void update(double T) override {
        m_base.update(T);
    }

    void update(double T, double extra) override {
        m_base.update(T, extra);
    }

    void update(double T, const vector<double>& extra) override {
        m_base.update(T, extra);
    }

    bool update(const ThermoPhase& phase, const Kinetics& kin) override {
        bool changed = false;
        if (m_nReactions != m_kf_const.rows()) {
            rebuildInterpData(phase, kin);
            changed = true;
        }

        double Tnow = phase.temperature();
        if (Tnow != Tprev) {
            changed = true;
            Tprev = Tnow;
        }

        m_T_index = static_cast<size_t>(floor(
            (Tnow - m_Tmin) / (m_Tmax - m_Tmin) * (m_nTemps-1)));
        m_dT = Tnow - m_T_index * (m_Tmax - m_Tmin) / (m_nTemps - 1) - m_Tmin;
        assert(m_dT >= 0 && m_dT <= (m_Tmax - m_Tmin) / (m_nTemps - 1));
        assert(m_T_index < m_nTemps);

        m_base.update(phase, kin);
        return changed;
    }

    double evalSingle(ReactionRate& rate) override {
        return m_base.evalSingle(rate);
    }

private:
    void rebuildInterpData(const ThermoPhase& phase, const Kinetics& kin);

    Cantera::MultiRate<Cantera::ArrheniusRate, Cantera::ArrheniusData> m_base;
    vector<pair<size_t, size_t>> m_indices; //! list of (global index, local index)

    size_t m_nReactions = 0;
    size_t m_nTemps = 200; //!< number of temperatures over which to interpolate
    double m_Tmin = 273; //!< minimum temperature
    double m_Tmax = 3500; //!< maximum temperature

    dmatrix m_kf_const, m_kf_slope; // Interpolation data

    // Interpolation parameters for the current temperature
    size_t m_T_index;
    double m_dT;
    double Tprev = 0;
};


//! A set of Cantera objects needed for calculating thermodynamic properties,
//! transport properties, and kinetic rates for a constant-pressure mixture.
class CanteraGas
{
public:
    CanteraGas();

    double pressure; //!< thermodynamic pressure [Pa]
    size_t nSpec; //!< number of species

    void setOptions(const ConfigOptions& options);
    void initialize();
    bool initialized() const { return isInitialized; }

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

    //! Multiply the rate constant for all reactions by *m*.
    void setRateMultiplier(double m);

    //! Get net molar reaction rates for each species
    void getReactionRates(dvec& wDot) const;
    void getReactionRates(double* wDot) const;

    void getCreationRates(dvec& wDot) const;
    void getCreationRates(double* wDot) const;

    void getDestructionRates(dvec& wDot) const;
    void getDestructionRates(double* wDot) const;

    std::shared_ptr<Cantera::ThermoPhase> thermo;
private:
    std::shared_ptr<Cantera::Solution> soln;
    std::string mechanismFile;
    std::string phaseID;
    std::string transportModel;
    std::string kineticsModel;
    double transportThreshold;
    std::shared_ptr<Cantera::Kinetics> kinetics;
    std::shared_ptr<Cantera::Transport> transport;
    bool isInitialized;

    dmatrix Dbin; //!< binary diffusion coefficients for species k

    //! Value of the rate multiplier the last time it was set
    double lastRateMultiplier;
};
