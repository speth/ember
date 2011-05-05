#include "chemistry0d.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

CanteraGas::CanteraGas()
    : rootXmlNode(NULL)
    , phaseXmlNode(NULL)
    , kinetics(NULL)
    , transport(NULL)
{
}

CanteraGas::~CanteraGas()
{
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
    // XML Information File
    if (!boost::filesystem::exists(mechanismFile)) {
        throw debugException((format(
            "Error: Cantera input file '%s' not found.") % mechanismFile).str());
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
    } else {
        transport = transFac->newTransport("Mix",&thermo);
    }
    transFac->initTransport(transport, &thermo);
    transFac->deleteFactory();

    nSpec = thermo.nSpecies();
    Dbin.resize(nSpec,nSpec,0);
    X.resize(nSpec,0);
    Y.resize(nSpec,0);
}

void CanteraGas::setStateMass(const dvector& Y_in, const double T)
{
    setStateMass(&Y_in[0], T);
}

void CanteraGas::setStateMass(const double* Y_in, const double T)
{
    for (size_t k=0; k<nSpec; k++) {
        Y[k] = max(Y_in[k], 0.0);
    }
    thermo.setState_TPY(T, pressure, &Y[0]);
}

void CanteraGas::setStateMole(const dvector& X_in, const double T)
{
    setStateMole(&X_in[0], T);
}

void CanteraGas::setStateMole(const double* X_in, const double T)
{
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
    getMolecularWeights(&W[0]);
}

void CanteraGas::getMolecularWeights(double* W) const
{
    for (size_t k=0; k<nSpec; k++) {
        W[k] = thermo.molecularWeight(k);
    }
}

double CanteraGas::getViscosity() const
{
    return transport->viscosity();
}

double CanteraGas::getThermalConductivity() const
{
    return transport->thermalConductivity();
}

void CanteraGas::getDiffusionCoefficientsMole(dvector& Dkm) const
{
    transport->getMixDiffCoeffs(&Dkm[0]);
}

void CanteraGas::getDiffusionCoefficientsMole(double* Dkm) const
{
    transport->getMixDiffCoeffs(Dkm);
}


void CanteraGas::getWeightedDiffusionCoefficientsMole(dvector& rhoD) const
{
    getWeightedDiffusionCoefficientsMole(&rhoD[0]);
}

void CanteraGas::getWeightedDiffusionCoefficientsMole(double* rhoD) const
{
    transport->getMixDiffCoeffs(rhoD);
    double rho = thermo.density();
    for (size_t k=0; k<nSpec; k++) {
        rhoD[k] *= rho;
    }
}

void CanteraGas::getWeightedDiffusionCoefficientsMass(double* rhoD)
{
    thermo.getMassFractions(&Y[0]);
    thermo.getMoleFractions(&X[0]);
    transport->getBinaryDiffCoeffs(nSpec, &Dbin(0,0));
    double rho = thermo.density();
    double eps = 1e-14/Y.size(); // Avoid 0/0 as Y[k] -> 1
    double Keps = 1e-14;
    // See Kee, p. 528, Eq. 12.178
    for (size_t k=0; k<nSpec; k++) {
        double sum1 = 0;
        double sum2 = 0;

        for (size_t i=0; i<nSpec; i++) {
            if (i==k) {
                continue;
            }
            sum1 += X[i]/Dbin(k,i);
            sum2 += (Y[i]+eps)/Dbin(k,i);
        }
        rhoD[k] = rho/(sum1 + X[k]/(1+Keps-Y[k])*sum2);
    }
}

void CanteraGas::getWeightedDiffusionCoefficientsMass(dvector& rhoD)
{
    getWeightedDiffusionCoefficientsMass(&rhoD[0]);
}

void CanteraGas::getThermalDiffusionCoefficients(dvector& Dkt) const
{
    transport->getThermalDiffCoeffs(&Dkt[0]);
}

void CanteraGas::getThermalDiffusionCoefficients(double* Dkt) const
{
    transport->getThermalDiffCoeffs(Dkt);
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

void CanteraGas::getCreationRates(dvector& wDot) const
{
    kinetics->getCreationRates(&wDot[0]);
}

void CanteraGas::getCreationRates(double* wDot) const
{
    kinetics->getCreationRates(wDot);
}

void CanteraGas::getDestructionRates(dvector& wDot) const
{
    kinetics->getDestructionRates(&wDot[0]);
}

void CanteraGas::getDestructionRates(double* wDot) const
{
    kinetics->getDestructionRates(wDot);
}
