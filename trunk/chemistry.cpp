#include "chemistry.h"
#include "boost/filesystem.hpp"

gasArray::gasArray()
	: rootXmlNode(NULL)
	, phaseXmlNode(NULL)
	, m_kineticsBase(NULL)
	, m_MultiTransportBase(NULL)
	, m_MixTransportBase(NULL)
{
}

gasArray::~gasArray()
{
	Cantera::close_XML_File(mechanismFile);
	//phaseXmlNode->unlock();
	//rootXmlNode->unlock();
//	delete rootXmlNode; // this deletes all child nodes as well

	delete m_kineticsBase;
	delete m_MultiTransportBase;
	delete m_MixTransportBase;

	resize(0);
}

void gasArray::resize(unsigned int n)
{
	if (n > m_thermo.size()) {
		// insert elements into the arrays
		while (m_thermo.size() != n) {
			m_thermo.insert(m_thermo.end(),new Cantera::IdealGasPhase);
			**m_thermo.rbegin() = m_thermoBase;

			m_kinetics.insert(m_kinetics.end(), new Cantera::GasKinetics(*m_thermo.rbegin()));
			(**m_kinetics.rbegin()).init();
			Cantera::installReactionArrays(*phaseXmlNode, **m_kinetics.rbegin(), phaseID);
			(**m_kinetics.rbegin()).finalize();

			if (usingMultiTransport) {
				m_MultiTransport.insert(m_MultiTransport.end(), new Cantera::MultiTransport);
				(**m_MultiTransport.rbegin()) = *m_MultiTransportBase;
				(**m_MultiTransport.rbegin()).m_thermo = *m_thermo.rbegin();
			} else {
				m_MixTransport.insert(m_MixTransport.end(), new Cantera::MixTransport);
				(**m_MixTransport.rbegin()) = *m_MixTransportBase;
				(**m_MixTransport.rbegin()).m_thermo = *m_thermo.rbegin();

			}
		}

	} else if (2*n < m_thermo.size()) {
		// remove elements from each of the arrays
		vector<Cantera::GasKinetics*>::iterator iter;
		for (iter=m_kinetics.begin()+n; iter!=m_kinetics.end(); iter++) {
			delete *iter;
		}
		m_kinetics.erase(m_kinetics.begin()+n,m_kinetics.end());

		vector<Cantera::IdealGasPhase*>::iterator iter2;
		for (iter2=m_thermo.begin()+n; iter2!=m_thermo.end(); iter2++) {
			delete *iter2;
		}
		m_thermo.erase(m_thermo.begin()+n,m_thermo.end());

		if (usingMultiTransport) {
			vector<Cantera::MultiTransport*>::iterator iter3;
			for (iter3=m_MultiTransport.begin()+n; iter3!=m_MultiTransport.end(); iter3++) {
				delete *iter3;
			}
			m_MultiTransport.erase(m_MultiTransport.begin()+n,m_MultiTransport.end());
		} else {
			vector<Cantera::MixTransport*>::iterator iter3;
			for (iter3=m_MixTransport.begin()+n; iter3!=m_MixTransport.end(); iter3++) {
				delete *iter3;
			}
			m_MixTransport.erase(m_MixTransport.begin()+n,m_MixTransport.end());

		}
	}
	nPoints = n;
	nSpec = m_thermoBase.nSpecies();
}

void gasArray::initialize(bool multiTransportFlag)
{
	usingMultiTransport = multiTransportFlag;

	// XML Information File
	if (!boost::filesystem::exists(mechanismFile)) {
		cout << "Error: Cantera input file \"" << mechanismFile << "\" not found." << endl;
		throw;
	}

	rootXmlNode = Cantera::get_XML_File(mechanismFile);
	phaseXmlNode = rootXmlNode->findNameID("phase",phaseID);

	// Initialize the default thermodynamic properties object
	Cantera::importPhase(*phaseXmlNode, &m_thermoBase);

	// Initialize the default chemical kinetics object
	m_kineticsBase = new Cantera::GasKinetics(&m_thermoBase);
	m_kineticsBase->init();
	Cantera::installReactionArrays(*phaseXmlNode, *m_kineticsBase, phaseID);
	m_kineticsBase->finalize();

	// Initialize the default transport properties object
	Cantera::TransportFactory* transFac = Cantera::TransportFactory::factory();
	if (usingMultiTransport) {
		m_MultiTransportBase = new Cantera::MultiTransport;
		transFac->initTransport(m_MultiTransportBase,&m_thermoBase);
	} else {
		m_MixTransportBase = new Cantera::MixTransport;
		transFac->initTransport(m_MixTransportBase,&m_thermoBase);
	}
	transFac->deleteFactory();
}

void gasArray::setStateMass(Cantera::Array2D& Y, dvector& T)
{
	dvector Yj(nSpec);
	for (int j=0; j<nPoints; j++) {
		for (int k=0; k<nSpec; k++) {
			Yj[k] = max(Y(k,j),0.0);
		}
		m_thermo[j]->setState_TPY(T[j], pressure, &Yj[0]);
	}
}

void gasArray::setStateMole(Cantera::Array2D& X, dvector& T)
{
	for (int j=0; j<nPoints; j++) {
		m_thermo[j]->setState_TPX(T[j], pressure, &X(0,j));
	}
}

void gasArray::getMoleFractions(Cantera::Array2D& X)
{
	for (int j=0; j<nPoints; j++) {
		m_thermo[j]->getMoleFractions(&X(0,j));
	}
}

void gasArray::getMassFractions(Cantera::Array2D& Y)
{
	for (int j=0; j<nPoints; j++) {
		m_thermo[j]->getMassFractions(&Y(0,j));
	}
}

void gasArray::getDensity(dvector& rho)
{
	for (int j=0; j<nPoints; j++) {
		rho[j] = m_thermo[j]->density();
	}
}

void gasArray::getMixtureMolecularWeight(dvector& Wmx)
{
	for (int j=0; j<nPoints; j++) {
		Wmx[j] = m_thermo[j]->meanMolecularWeight();
	}
}

void gasArray::getMolecularWeights(dvector& W)
{
	for (int k=0; k<nSpec; k++) {
		W[k] = m_thermoBase.molecularWeight(k);
	}
}

void gasArray::getViscosity(dvector& mu)
{
	if (usingMultiTransport) {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			mu[j] = m_MultiTransport[j]->viscosity();
		}
	} else {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			mu[j] = m_MixTransport[j]->viscosity();
		}
	}
}

void gasArray::getThermalConductivity(dvector& lambda)
{
	if (usingMultiTransport) {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			lambda[j] = m_MultiTransport[j]->thermalConductivity();
		}
	} else {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			lambda[j] = m_MixTransport[j]->thermalConductivity();
		}
	}

}

void gasArray::getDiffusionCoefficients(Cantera::Array2D& Dkm)
{
	if (usingMultiTransport) {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			m_MultiTransport[j]->getMixDiffCoeffs(&Dkm(0,j));
		}
	} else {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			m_MixTransport[j]->getMixDiffCoeffs(&Dkm(0,j));
		}
	}
}

void gasArray::getWeightedDiffusionCoefficients(Cantera::Array2D& rhoD)
{
	if (usingMultiTransport) {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			m_MultiTransport[j]->getMixDiffCoeffs(&rhoD(0,j));
			double rho = m_thermo[j]->density();
			for (int k=0; k<nSpec; k++) {
				rhoD(k,j) *= rho;
			}
		}
	} else {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			m_MixTransport[j]->getMixDiffCoeffs(&rhoD(0,j));
			double rho = m_thermo[j]->density();
			for (int k=0; k<nSpec; k++) {
				rhoD(k,j) *= rho;
			}
		}
	}
}

void gasArray::getThermalDiffusionCoefficients(Cantera::Array2D& Dkt)
{
	if (usingMultiTransport) {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			m_MultiTransport[j]->getThermalDiffCoeffs(&Dkt(0,j));
		}
	} else {
		#pragma omp parallel for
		for (int j=0; j<nPoints; j++) {
			m_MixTransport[j]->getThermalDiffCoeffs(&Dkt(0,j));
		}
	}
}

void gasArray::getSpecificHeatCapacity(dvector& cp)
{
	for (int j=0; j<nPoints; j++) {
		cp[j] = m_thermo[j]->cp_mass();
	}
}

void gasArray::getSpecificHeatCapacities(Cantera::Array2D& cpSpec)
{
	for (int j=0; j<nPoints; j++) {
		m_thermo[j]->getPartialMolarCp(&cpSpec(0,j));
	}
}

void gasArray::getEnthalpies(Cantera::Array2D& hk)
{
	for (int j=0; j<nPoints; j++) {
		m_thermo[j]->getPartialMolarEnthalpies(&hk(0,j));
	}
}

void gasArray::getReactionRates(Cantera::Array2D& wDot)
{
	#pragma omp parallel for
	for (int j=0; j<nPoints; j++) {
		m_kinetics[j]->getNetProductionRates(&wDot(0,j));
	}
}

Cantera::IdealGasPhase& gasArray::operator[](unsigned int i) const
{
	return *m_thermo[i];
}

Cantera::IdealGasPhase& gasArray::thermo(unsigned int i) const
{
	return *m_thermo[i];
}

Cantera::GasKinetics& gasArray::kinetics(unsigned int i) const
{
	return *m_kinetics[i];
}

Cantera::MultiTransport& gasArray::multiTrans(unsigned int i) const
{
	return *m_MultiTransport[i];
}

Cantera::MixTransport& gasArray::mixTrans(unsigned int i) const
{
	return *m_MixTransport[i];
}
