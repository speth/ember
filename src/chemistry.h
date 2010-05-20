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

class canteraGas
{
    // This class groups together a set of Cantera objects needed for calculating
    // thermodynamic properties, transport properties, and kinetic rates

public:
    canteraGas();
    ~canteraGas();

    std::string mechanismFile;
    std::string phaseID;
    double pressure; // thermodynamic pressure
    int nSpec; // number of species

    bool usingMultiTransport; // use multicomponent transport model? (vs. mixture-averaged)

    void initialize(bool multiTransportFlag);
    void resize(unsigned int n);

    void setStateMass(Array2D& Y, dvector& T);
    void setStateMole(Array2D& X, dvector& T);

    void getMoleFractions(Array2D& X);
    void getMassFractions(Array2D& Y);
    void getDensity(dvector& rho);
    void getMixtureMolecularWeight(dvector& Wmx);
    void getMolecularWeights(dvector& W);

    void getViscosity(dvector& mu);
    void getThermalConductivity(dvector& lambda);
    void getDiffusionCoefficients(Array2D& Dkm);
    void getWeightedDiffusionCoefficients(Array2D& rhoD);
    void getThermalDiffusionCoefficients(Array2D& Dkt);

    void getSpecificHeatCapacity(dvector& cp);
    void getSpecificHeatCapacities(Array2D& cpSpec);
    void getEnthalpies(Array2D& hk);

    void getReactionRates(Array2D& wDot);

    void getTransportProperties(dvector& mu, dvector& lambda, Array2D& rhoD, Array2D& Dkt);
    void getThermoProperties(dvector& rho, dvector& Wmx, dvector& cp, Array2D& cpSpec, Array2D& hk);

    Cantera::IdealGasPhase thermo;

private:
    Cantera::XML_Node* rootXmlNode;
    Cantera::XML_Node* phaseXmlNode;

    int nPoints;
    Array2D Y;
    dvector T;

    Cantera::GasKinetics* kinetics;
    Cantera::MultiTransport* multiTransport;
    Cantera::MixTransport* mixTransport;
};
