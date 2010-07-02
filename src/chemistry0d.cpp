#include "chemistry0d.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

CanteraGas::CanteraGas()
    : rootXmlNode(NULL)
    , phaseXmlNode(NULL)
    , kinetics(NULL)
    , transport(NULL)
{
    std:: cout << "CanteraGas::CanteraGas()" << std::endl;
}

CanteraGas::~CanteraGas()
{
    std:: cout << "CanteraGas::~CanteraGas()" << std::endl;
    Cantera::close_XML_File(mechanismFile);
    //phaseXmlNode->unlock();
    //rootXmlNode->unlock();
//    delete rootXmlNode; // this deletes all child nodes as well

    delete kinetics;
    delete transport;
}

void CanteraGas::setOptions(const configOptions& options)
{
    usingMultiTransport = options.usingMultiTransport;
    mechanismFile = options.gasMechanismFile;
    phaseID = options.gasPhaseID;
    pressure = options.pressure;
}

void CanteraGas::initialize()
{
    std:: cout << "CanteraGas::initialize()" << std::endl;
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
        transport = transFac->newTransport("Multi",&thermo);
        transFac->initTransport(transport, &thermo);
    } else {
        transport = transFac->newTransport("Mix",&thermo);
        transFac->initTransport(transport, &thermo);
    }
    transFac->deleteFactory();

    nSpec = thermo.nSpecies();
}

void CanteraGas::setStateMass(const dvector& Y_in, const double T)
{
    dvector Y(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        Y[k] = max(Y_in[k], 0.0);
    }
    thermo.setState_TPY(T, pressure, &Y[0]);
}

void CanteraGas::setStateMass(const double* Y_in, const double T)
{
    dvector Y(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        Y[k] = max(Y_in[k], 0.0);
    }
    thermo.setState_TPY(T, pressure, &Y[0]);
}

void CanteraGas::setStateMole(const dvector& X_in, const double T)
{
    dvector X(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        X[k] = max(X_in[k], 0.0);
    }
    thermo.setState_TPX(T, pressure, &X[0]);
}

void CanteraGas::setStateMole(const double* X_in, const double T)
{
    dvector X(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        X[k] = max(X_in[k], 0.0);
    }
    thermo.setState_TPX(T, pressure, &X[0]);
}

void CanteraGas::getMoleFractions(dvector& X) const
{
    thermo.getMoleFractions(&X[0]);
}

void CanteraGas::getMoleFractions(double* X) const
{
    thermo.getMoleFractions(X);
}

void CanteraGas::getMassFractions(double* Y) const
{
    thermo.getMassFractions(Y);
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

void CanteraGas::getMolecularWeights(double* W) const
{
    for (size_t k=0; k<nSpec; k++) {
        W[k] = thermo.molecularWeight(k);
    }
}

double CanteraGas::getViscosity() const
{
    if (usingMultiTransport) {
        return transport->viscosity();
    } else {
        return transport->viscosity();
    }
}

double CanteraGas::getThermalConductivity() const
{
    if (usingMultiTransport) {
        return transport->thermalConductivity();
    } else {
        return transport->thermalConductivity();
    }
}

void CanteraGas::getDiffusionCoefficients(dvector& Dkm) const
{
    if (usingMultiTransport) {
        transport->getMixDiffCoeffs(&Dkm[0]);
    } else {
        transport->getMixDiffCoeffs(&Dkm[0]);
    }
}

void CanteraGas::getDiffusionCoefficients(double* Dkm) const
{
    if (usingMultiTransport) {
        transport->getMixDiffCoeffs(Dkm);
    } else {
        transport->getMixDiffCoeffs(Dkm);
    }
}


void CanteraGas::getWeightedDiffusionCoefficients(dvector& rhoD) const
{
    if (usingMultiTransport) {
        transport->getMixDiffCoeffs(&rhoD[0]);
        double rho = thermo.density();
        for (size_t k=0; k<nSpec; k++) {
            rhoD[k] *= rho;
        }
    } else {
        transport->getMixDiffCoeffs(&rhoD[0]);
        double rho = thermo.density();
        for (size_t k=0; k<nSpec; k++) {
            rhoD[k] *= rho;
        }
    }
}

void CanteraGas::getWeightedDiffusionCoefficients(double* rhoD) const
{
    if (usingMultiTransport) {
        transport->getMixDiffCoeffs(rhoD);
        double rho = thermo.density();
        for (size_t k=0; k<nSpec; k++) {
            rhoD[k] *= rho;
        }
    } else {
        transport->getMixDiffCoeffs(rhoD);
        double rho = thermo.density();
        for (size_t k=0; k<nSpec; k++) {
            rhoD[k] *= rho;
        }
    }
}

void CanteraGas::getThermalDiffusionCoefficients(dvector& Dkt) const
{
    if (usingMultiTransport) {
        transport->getThermalDiffCoeffs(&Dkt[0]);
    } else {
        transport->getThermalDiffCoeffs(&Dkt[0]);
    }
}

void CanteraGas::getThermalDiffusionCoefficients(double* Dkt) const
{
    if (usingMultiTransport) {
        transport->getThermalDiffCoeffs(Dkt);
    } else {
        transport->getThermalDiffCoeffs(Dkt);
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

void CanteraGas::getSpecificHeatCapacities(double* cpSpec) const
{
    thermo.getPartialMolarCp(cpSpec);
}

void CanteraGas::getEnthalpies(dvector& hk) const
{
    thermo.getPartialMolarEnthalpies(&hk[0]);
}

void CanteraGas::getEnthalpies(double* hk) const
{
    thermo.getPartialMolarEnthalpies(&hk[0]);
}

void CanteraGas::getReactionRates(dvector& wDot) const
{
    kinetics->getNetProductionRates(&wDot[0]);
}

void CanteraGas::getReactionRates(double* wDot) const
{
    kinetics->getNetProductionRates(wDot);
}
