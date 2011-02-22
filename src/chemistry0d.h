#pragma once

#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/equilibrium.h>    // chemical equilibrium
#include <cantera/thermo.h>
#include <cantera/transport.h>      // transport properties
#include <cantera/kinetics.h>

#include <cantera/kernel/IdealGasPhase.h>
#include <cantera/kernel/GasKinetics.h>
#include <cantera/kernel/ctml.h>

#include "mathUtils.h"
#include "readConfig.h"

#ifdef WIN32
#define Cantera_CXX Cantera
#endif

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

    bool usingMultiTransport; // use multicomponent transport model? (vs. mixture-averaged)

    Cantera::XML_Node* rootXmlNode;
    Cantera::XML_Node* phaseXmlNode;

    Cantera::GasKinetics* kinetics;
    Cantera::Transport* transport;

    Array2D Dbin; // binary diffusion coefficients for species k
    dvector X; // mole fractions
    dvector Y; // mass fractions
};
