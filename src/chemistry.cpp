#include "chemistry.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

canteraGas::canteraGas()
    : rootXmlNode(NULL)
    , phaseXmlNode(NULL)
    , kinetics(NULL)
    , multiTransport(NULL)
    , mixTransport(NULL)
{
}

canteraGas::~canteraGas()
{
    Cantera::close_XML_File(mechanismFile);
    //phaseXmlNode->unlock();
    //rootXmlNode->unlock();
//    delete rootXmlNode; // this deletes all child nodes as well

    delete kinetics;
    delete multiTransport;
    delete mixTransport;

    resize(0);
}

void canteraGas::resize(unsigned int n)
{
    Y.resize(nSpec,n);
    T.resize(n);
    nPoints = n;
}

void canteraGas::initialize(bool multiTransportFlag)
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

void canteraGas::setStateMass(Array2D& Y_in, dvector& T_in)
{
    for (int j=0; j<nPoints; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = max(Y_in(k,j),0.0);
        }
        T[j] = T_in[j];
    }
}

void canteraGas::setStateMole(Array2D& X, dvector& T_in)
{
    dvector Xj(nSpec);
    for (int j=0; j<nPoints; j++) {
        for (int k=0; k<nSpec; k++) {
            Xj[k] = max(X(k,j),0.0);
        }
        thermo.setState_TPX(T[j], pressure, &Xj[0]);
        thermo.getMassFractions(&Y(0,j));
        T[j] = T_in[j];
    }
}

void canteraGas::getMoleFractions(Array2D& X)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setMassFractions(&Y(0,j));
        thermo.getMoleFractions(&X(0,j));
    }
}

void canteraGas::getMassFractions(Array2D& Y_out)
{
    for (int j=0; j<nPoints; j++) {
        for (int k=0; k<nSpec; k++) {
            Y_out(k,j) = Y(k,j);
        }
    }
}

void canteraGas::getDensity(dvector& rho)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setState_TPY(T[j], pressure, &Y(0,j));
        rho[j] = thermo.density();
    }
}

void canteraGas::getMixtureMolecularWeight(dvector& Wmx)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setState_TPY(T[j], pressure, &Y(0,j));
        Wmx[j] = thermo.meanMolecularWeight();
    }
}

void canteraGas::getMolecularWeights(dvector& W)
{
    for (int k=0; k<nSpec; k++) {
        W[k] = thermo.molecularWeight(k);
    }
}

void canteraGas::getViscosity(dvector& mu)
{
    if (usingMultiTransport) {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            mu[j] = multiTransport->viscosity();
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            mu[j] = mixTransport->viscosity();
        }
    }
}

void canteraGas::getThermalConductivity(dvector& lambda)
{
    if (usingMultiTransport) {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            lambda[j] = multiTransport->thermalConductivity();
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            lambda[j] = mixTransport->thermalConductivity();
        }
    }
}

void canteraGas::getDiffusionCoefficients(Array2D& Dkm)
{
    if (usingMultiTransport) {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            multiTransport->getMixDiffCoeffs(&Dkm(0,j));
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            mixTransport->getMixDiffCoeffs(&Dkm(0,j));
        }
    }
}

void canteraGas::getWeightedDiffusionCoefficients(Array2D& rhoD)
{
    if (usingMultiTransport) {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            multiTransport->getMixDiffCoeffs(&rhoD(0,j));
            double rho = thermo.density();
            for (int k=0; k<nSpec; k++) {
                rhoD(k,j) *= rho;
            }
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            mixTransport->getMixDiffCoeffs(&rhoD(0,j));
            double rho = thermo.density();
            for (int k=0; k<nSpec; k++) {
                rhoD(k,j) *= rho;
            }
        }
    }
}

void canteraGas::getThermalDiffusionCoefficients(Array2D& Dkt)
{
    if (usingMultiTransport) {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            multiTransport->getThermalDiffCoeffs(&Dkt(0,j));
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            mixTransport->getThermalDiffCoeffs(&Dkt(0,j));
        }
    }
}

void canteraGas::getSpecificHeatCapacity(dvector& cp)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setState_TPY(T[j], pressure, &Y(0,j));
        cp[j] = thermo.cp_mass();
    }
}

void canteraGas::getSpecificHeatCapacities(Array2D& cpSpec)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setState_TPY(T[j], pressure, &Y(0,j));
        thermo.getPartialMolarCp(&cpSpec(0,j));
    }
}

void canteraGas::getEnthalpies(Array2D& hk)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setState_TPY(T[j], pressure, &Y(0,j));
        thermo.getPartialMolarEnthalpies(&hk(0,j));
    }
}

void canteraGas::getReactionRates(Array2D& wDot)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setState_TPY(T[j], pressure, &Y(0,j));
        kinetics->getNetProductionRates(&wDot(0,j));
    }
}

void canteraGas::getTransportProperties(dvector& mu, dvector& lambda, Array2D& rhoD, Array2D& Dkt)
{
    if (usingMultiTransport) {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            mu[j] = multiTransport->viscosity();
            lambda[j] = multiTransport->thermalConductivity();
            multiTransport->getMixDiffCoeffs(&rhoD(0,j));
            double rho = thermo.density();
            for (int k=0; k<nSpec; k++) {
                rhoD(k,j) *= rho;
            }
            multiTransport->getThermalDiffCoeffs(&Dkt(0,j));
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            thermo.setState_TPY(T[j], pressure, &Y(0,j));
            mu[j] = mixTransport->viscosity();
            lambda[j] = mixTransport->thermalConductivity();
            mixTransport->getMixDiffCoeffs(&rhoD(0,j));
            double rho = thermo.density();
            for (int k=0; k<nSpec; k++) {
                rhoD(k,j) *= rho;
            }
            mixTransport->getThermalDiffCoeffs(&Dkt(0,j));
        }
    }
}

void canteraGas::getThermoProperties(dvector& rho, dvector& Wmx, dvector& cp, Array2D& cpSpec, Array2D& hk)
{
    for (int j=0; j<nPoints; j++) {
        thermo.setState_TPY(T[j], pressure, &Y(0,j));
        rho[j] = thermo.density();
        Wmx[j] = thermo.meanMolecularWeight();
        cp[j] = thermo.cp_mass();
        thermo.getPartialMolarCp(&cpSpec(0,j));
        thermo.getPartialMolarEnthalpies(&hk(0,j));
    }
}
