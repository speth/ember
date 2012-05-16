#pragma once

#include "cantera/IdealGasMix.h"    // defines class IdealGasMix

#include "cantera/equilibrium.h"    // chemical equilibrium
#include "cantera/thermo.h"
#include "cantera/transport.h"      // transport properties
#include "cantera/kinetics.h"

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/base/ctml.h"

#include "mathUtils.h"
#include "readConfig.h"

class ApproxMixTransport : public Cantera::MixTransport
{
public:
    ApproxMixTransport(Cantera::ThermoPhase& thermo,
                       Cantera::TransportFactory& factory);
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
    vector<size_t> _kMajor; // indices of the species where X[k] >= threshold
};

class CanteraGas
{
    // This class groups together a set of Cantera objects needed for calculating
    // thermodynamic properties, transport properties, and kinetic rates

public:
    CanteraGas();
    ~CanteraGas();

    double pressure; // thermodynamic pressure
    size_t nSpec; // number of species

    void setOptions(const configOptions& options);
    void initialize();
    bool initialized() const { return isInitialized; }

    //! Calculate the mole fractions of the reactant mixture from the
    //! compositions of the fuel and oxidizer mixtures and the equivalence ratio.
    dvector calculateReactantMixture(const std::string& fuel,
                                     const std::string& oxidizer,
                                     double equivalenceRatio);

    void setStateMass(const dvector& Y, const double T);
    void setStateMass(const double* Y, const double T);

    void setStateMole(const dvector& X, const double T);
    void setStateMole(const double* X, const double T);

    void getMoleFractions(dvector& X) const;
    void getMoleFractions(double* X) const;
    void getMassFractions(dvector& Y) const;
    void getMassFractions(double* Y) const;

    double getDensity() const;
    double getMixtureMolecularWeight() const;
    void getMolecularWeights(dvector& W) const;
    void getMolecularWeights(double* W) const;

    double getViscosity() const;
    double getThermalConductivity() const;

    // Diffusion coefficients for calculating v[k] = -D[k] / X[k] * grad(X[k])
    void getDiffusionCoefficientsMole(dvector& Dkm) const;
    void getDiffusionCoefficientsMole(double* Dkm) const;
    void getWeightedDiffusionCoefficientsMole(dvector& rhoD) const;
    void getWeightedDiffusionCoefficientsMole(double* rhoD) const;

    // Diffusion coefficients for calculating j[k] = - rho * D[k] * grad(Y[k])
    void getWeightedDiffusionCoefficientsMass(double* rhoD);
    void getWeightedDiffusionCoefficientsMass(dvector& rhoD);

    void getThermalDiffusionCoefficients(dvector& Dkt) const;
    void getThermalDiffusionCoefficients(double* Dkt) const;

    double getSpecificHeatCapacity() const;
    void getSpecificHeatCapacities(dvector& cpSpec) const;
    void getSpecificHeatCapacities(double* cpSpec) const;
    void getEnthalpies(dvector& hk) const;
    void getEnthalpies(double* hk) const;

    void getReactionRates(dvector& wDot) const;
    void getReactionRates(double* wDot) const;

    void getCreationRates(dvector& wDot) const;
    void getCreationRates(double* wDot) const;
    void getDestructionRates(dvector& wDot) const;
    void getDestructionRates(double* wDot) const;

    Cantera::IdealGasPhase thermo;

private:
    std::string mechanismFile;
    std::string phaseID;
    std::string transportModel;
    double transportThreshold;
    bool isInitialized;

    Cantera::XML_Node* rootXmlNode;
    Cantera::XML_Node* phaseXmlNode;

    Cantera::GasKinetics* kinetics;
    Cantera::GasTransport* transport;

    dmatrix Dbin; // binary diffusion coefficients for species k
    dvector X; // mole fractions
    dvector Y; // mass fractions
    vector<size_t> kMajor; // indices of species with mole fractions above some threshold
};
