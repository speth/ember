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

    std::string mechanismFile;
    std::string phaseID;
    double pressure; // thermodynamic pressure
    size_t nSpec; // number of species

    bool usingMultiTransport; // use multicomponent transport model? (vs. mixture-averaged)

    void initialize(bool multiTransportFlag);

    void setStateMass(const dvector& Y, const double T);
    void setStateMass(const double* Y, const double T);
    void setStateMole(const dvector& X, const double T);

    void getMoleFractions(dvector& X) const;
    void getMassFractions(dvector& Y) const;
    double getDensity() const;
    double getMixtureMolecularWeight() const;
    void getMolecularWeights(dvector& W) const;

    double getViscosity() const;
    double getThermalConductivity() const;
    void getDiffusionCoefficients(dvector& Dkm) const;
    void getWeightedDiffusionCoefficients(dvector& rhoD) const;
    void getThermalDiffusionCoefficients(dvector& Dkt) const;

    double getSpecificHeatCapacity() const;
    void getSpecificHeatCapacities(dvector& cpSpec) const;
    void getEnthalpies(dvector& hk) const;

    void getReactionRates(dvector& wDot) const;

    Cantera::IdealGasPhase thermo;

private:
    Cantera::XML_Node* rootXmlNode;
    Cantera::XML_Node* phaseXmlNode;

    Cantera::GasKinetics* kinetics;
    Cantera::MultiTransport* multiTransport;
    Cantera::MixTransport* mixTransport;
};
