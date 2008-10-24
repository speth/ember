#include "chemistry.h"
#include "boost/filesystem.hpp"

gasArray::gasArray()
	: rootXmlNode(NULL)
	, phaseXmlNode(NULL)
	, m_kineticsBase(NULL)
	, m_transportBase(NULL)
{
}

gasArray::~gasArray()
{
	Cantera::close_XML_File(mechanismFile);
	//phaseXmlNode->unlock();
	//rootXmlNode->unlock();
//	delete rootXmlNode; // this deletes all child nodes as well

	delete m_kineticsBase;
	delete m_transportBase;

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

			//**m_kinetics.rbegin() = *m_kineticsBase;
			//(**m_kinetics.rbegin()).m_thermo[0] = *m_thermo.rbegin();

			m_transport.insert(m_transport.end(), new Cantera::MultiTransport);
			(**m_transport.rbegin()) = *m_transportBase;
			(**m_transport.rbegin()).m_thermo = *m_thermo.rbegin();
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

		vector<Cantera::MultiTransport*>::iterator iter3;
		for (iter3=m_transport.begin()+n; iter3!=m_transport.end(); iter3++) {
			delete *iter3;
		}
		m_transport.erase(m_transport.begin()+n,m_transport.end());
	}
	nPoints = n;
	nSpec = m_thermoBase.nSpecies();
}

void gasArray::initialize(void)
{
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
	m_transportBase = new Cantera::MultiTransport;
	Cantera::TransportFactory* transFac = Cantera::TransportFactory::factory();
	transFac->initTransport(m_transportBase,&m_thermoBase);
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
	#pragma omp parallel for
	for (int j=0; j<nPoints; j++) {
		mu[j] = m_transport[j]->viscosity();
	}
}

void gasArray::getThermalConductivity(dvector& lambda)
{
	#pragma omp parallel for
	for (int j=0; j<nPoints; j++) {
		lambda[j] = m_transport[j]->thermalConductivity();
	}
}

void gasArray::getDiffusionCoefficients(Cantera::Array2D& Dkm)
{
	#pragma omp parallel for
	for (int j=0; j<nPoints; j++) {
		m_transport[j]->getMixDiffCoeffs(&Dkm(0,j));
	}
}

void gasArray::getWeightedDiffusionCoefficients(Cantera::Array2D& rhoD)
{
	#pragma omp parallel for
	for (int j=0; j<nPoints; j++) {
		m_transport[j]->getMixDiffCoeffs(&rhoD(0,j));
		double rho = m_thermo[j]->density();
		for (int k=0; k<nSpec; k++) {
			rhoD(k,j) *= rho;
		}
	}
}

void gasArray::getThermalDiffusionCoefficients(Cantera::Array2D& Dkt)
{
	#pragma omp parallel for
	for (int j=0; j<nPoints; j++) {
		m_transport[j]->getThermalDiffCoeffs(&Dkt(0,j));
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

Cantera::MultiTransport& gasArray::trans(unsigned int i) const
{
	return *m_transport[i];
}

void gasArray::testFunction(void)
{
	Cantera::IdealGasPhase thermo1;
	Cantera::IdealGasPhase thermo2;

	Cantera::GasKinetics* kin1 = new Cantera::GasKinetics;
	Cantera::GasKinetics* kin2 = new Cantera::GasKinetics;

	Cantera::importPhase(*phaseXmlNode, &thermo1);
	Cantera::importPhase(*phaseXmlNode, &thermo2);

	kin1->addPhase(thermo1);
	kin2->addPhase(thermo1);

	kin1->init();
	kin2->init();

	Cantera::installReactionArrays(*phaseXmlNode,*kin1,phaseID);
	Cantera::installReactionArrays(*phaseXmlNode,*kin2,phaseID);

	kin1->finalize();
	kin2->finalize();

	kin2->m_thermo[0] = &thermo2;


	dvector wdot1(thermo1.nSpecies()), wdot2(thermo2.nSpecies());

	Cantera::MultiTransport* trans1 = new Cantera::MultiTransport();
	Cantera::MultiTransport* trans2 = new Cantera::MultiTransport();

	Cantera::TransportFactory* transFac = Cantera::TransportFactory::factory();
	transFac->initTransport(trans1,&thermo1);
	*trans2 = *trans1;
	trans2->m_thermo = &thermo2;


	thermo1.setState_TPX(300,101325,"O2:1.0, H2:0.5");
	thermo2.setState_TPX(300,101325,"AR:1.0");

	kin1->getNetProductionRates(&wdot1[0]);
	kin2->getNetProductionRates(&wdot2[0]);

	cout << wdot1 << endl;
	cout << wdot2 << endl;

	cout << trans1->viscosity() << endl;
	cout << trans2->viscosity() << endl;

	thermo1.setState_TPX(300,101325,"AR:1.0");
	thermo2.setState_TPX(300,101325,"O2:1.0, H2:0.5");

	cout << trans1->viscosity() << endl;
	cout << trans2->viscosity() << endl;

	delete kin1;
	delete kin2;

}
