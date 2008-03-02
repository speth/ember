#include "chemistry.h"
#include "mathUtils.h"

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
	phaseXmlNode->unlock();
	rootXmlNode->unlock();
	delete rootXmlNode; // this deletes all child nodes as well

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

			m_kinetics.insert(m_kinetics.end(), new Cantera::GasKinetics(&m_thermoBase));
			(**m_kinetics.rbegin()).m_thermo[0] = *m_thermo.rbegin();

			m_transport.insert(m_transport.end(), new Cantera::MultiTransport);
			(**m_transport.rbegin()) = *m_transportBase;
			(**m_transport.rbegin()).m_thermo = *m_thermo.rbegin();
		}

	} else if (n < m_thermo.size()) {
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
}

void gasArray::initialize(void) 
{
	// XML Information File
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

void gasArray::setState(Cantera::Array2D& Y, dvector& T)
{
	for (unsigned int j=0; j<Y.nColumns(); j++) {
		m_thermo[j]->setState_TPX(T[j], pressure, &Y(0,j));
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


	thermo1.setState_TPX(300,101325,"O2:1.0, CH4:0.5");
	thermo2.setState_TPX(300,101325,"AR:1.0");

	kin1->getNetProductionRates(&wdot1[0]);
	kin2->getNetProductionRates(&wdot2[0]);

	cout << wdot1 << endl;
	cout << wdot2 << endl;

	cout << trans1->viscosity() << endl;
	cout << trans2->viscosity() << endl;

	thermo1.setState_TPX(300,101325,"AR:1.0");
	thermo2.setState_TPX(300,101325,"O2:1.0, CH4:0.5");

	cout << trans1->viscosity() << endl;
	cout << trans2->viscosity() << endl;

	int blargh = 0;

		
	delete kin1;
	delete kin2;

}