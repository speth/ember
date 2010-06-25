#include "chemistry0d.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

CanteraGas::CanteraGas()
    : rootXmlNode(NULL)
    , phaseXmlNode(NULL)
    , kinetics(NULL)
    , multiTransport(NULL)
    , mixTransport(NULL)
{
}

CanteraGas::~CanteraGas()
{
    Cantera::close_XML_File(mechanismFile);
    //phaseXmlNode->unlock();
    //rootXmlNode->unlock();
//    delete rootXmlNode; // this deletes all child nodes as well

    delete kinetics;
    delete multiTransport;
    delete mixTransport;
}


void CanteraGas::initialize(bool multiTransportFlag)
{
    usingMultiTransport = multiTransportFlag;

    // XML Information File
    if (!boost::filesystem::exists(mechanismFile)) {
        throw debugException("Error: Cantera input file \"" + mechanismFile + "\" not found.");
    }

    rootXmlNode = Cantera::get_XML_File(mechanismFile);
    phaseXmlNode = rootXmlNode->findNameID("phase", phaseID);

    // Initialize the default thermodynamic properties object
    Cantera::importPhase(*phaseXmlNode, &thermo);

    // Initialize the default chemical kinetics object
    kinetics = new Cantera::GasKinetics(&thermo);
    kinetics->init();
    Cantera::installReactionArrays(*phaseXmlNode, *kinetics, phaseID);
    kinetics->finalize();

    // Initialize the default transport properties object
    Cantera::TransportFactory* transFac = Cantera::TransportFactory::factory();
    if (usingMultiTransport) {
        multiTransport = dynamic_cast<Cantera::MultiTransport*>(transFac->newTransport("Multi",&thermo));
        transFac->initTransport(multiTransport, &thermo);
    } else {
        mixTransport = dynamic_cast<Cantera::MixTransport*>(transFac->newTransport("Mix",&thermo));
        transFac->initTransport(mixTransport, &thermo);
    }
    transFac->deleteFactory();

    nSpec = thermo.nSpecies();
}

void CanteraGas::setStateMass(const dvector& Y_in, const double T)
{
    dvector Y(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        Y[k] = max(Y[k], 0.0);
    }
    thermo.setState_TPY(T, pressure, &Y[0]);
}

void CanteraGas::setStateMole(const dvector& X_in, const double T)
{
    dvector X(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        X[k] = max(X[k], 0.0);
    }
    thermo.setState_TPX(T, pressure, &X[0]);
}

void CanteraGas::getMoleFractions(dvector& X) const
{
    thermo.getMoleFractions(&X[0]);
}

void CanteraGas::getMassFractions(dvector& Y) const
{
    thermo.getMassFractions(&Y[0]);
}

double CanteraGas::getDensity() const
{
    return thermo.density();
}

double CanteraGas::getMixtureMolecularWeight() const
{
    return thermo.meanMolecularWeight();
}

void CanteraGas::getMolecularWeights(dvector& W) const
{
    for (size_t k=0; k<nSpec; k++) {
        W[k] = thermo.molecularWeight(k);
    }
}

double CanteraGas::getViscosity() const
{
    if (usingMultiTransport) {
        return multiTransport->viscosity();
    } else {
        return mixTransport->viscosity();
    }
}

double CanteraGas::getThermalConductivity() const
{
    if (usingMultiTransport) {
        return multiTransport->thermalConductivity();
    } else {
        return mixTransport->thermalConductivity();
    }
}

void CanteraGas::getDiffusionCoefficients(dvector& Dkm) const
{
    if (usingMultiTransport) {
        multiTransport->getMixDiffCoeffs(&Dkm[0]);
    } else {
        mixTransport->getMixDiffCoeffs(&Dkm[0]);
    }
}

void CanteraGas::getWeightedDiffusionCoefficients(dvector& rhoD) const
{
    if (usingMultiTransport) {
        multiTransport->getMixDiffCoeffs(&rhoD[0]);
        double rho = thermo.density();
        for (size_t k=0; k<nSpec; k++) {
            rhoD[k] *= rho;
        }
    } else {
        mixTransport->getMixDiffCoeffs(&rhoD[0]);
        double rho = thermo.density();
        for (size_t k=0; k<nSpec; k++) {
            rhoD[k] *= rho;
        }
    }
}

void CanteraGas::getThermalDiffusionCoefficients(dvector& Dkt) const
{
    if (usingMultiTransport) {
        multiTransport->getThermalDiffCoeffs(&Dkt[0]);
    } else {
        mixTransport->getThermalDiffCoeffs(&Dkt[0]);
    }
}

double CanteraGas::getSpecificHeatCapacity() const
{
    return thermo.cp_mass();
}

void CanteraGas::getSpecificHeatCapacities(dvector& cpSpec) const
{
    thermo.getPartialMolarCp(&cpSpec[0]);
}

void CanteraGas::getEnthalpies(dvector& hk) const
{
    thermo.getPartialMolarEnthalpies(&hk[0]);
}

void CanteraGas::getReactionRates(dvector& wDot) const
{
    kinetics->getNetProductionRates(&wDot[0]);
}
