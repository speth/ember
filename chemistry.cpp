#include "chemistry.h"

gasArray::gasArray()
	: rootXmlNode(NULL)
	, phaseXmlNode(NULL)
{
}

gasArray::~gasArray()
{
	Cantera::close_XML_File(mechanismFile);
	phaseXmlNode->unlock();
	rootXmlNode->unlock();
	delete rootXmlNode; // this deletes all child nodes as well
	
	resize(0);
}

void gasArray::resize(unsigned int n)
{
	if (n > d_thermo.size()) {
		// insert elements into the arrays
		while (d_thermo.size() != n) {
			d_thermo.insert(d_thermo.end(),new Cantera::IdealGasPhase());
			**(d_thermo.end()-1) = thermoBase;

			d_kinetics.insert(d_kinetics.end(), new Cantera::GasKinetics(*(d_thermo.end()-1)));
			(**(d_kinetics.end()-1)).init();
			Cantera::installReactionArrays(*phaseXmlNode,**(d_kinetics.end()-1),phaseID);
		}

	} else if (n < d_thermo.size()) {
		// remove elements of the arrays

		vector<Cantera::GasKinetics*>::iterator iter;
		for (iter=d_kinetics.begin()+n; iter!=d_kinetics.end(); iter++) {
			delete *iter;
		}
		d_kinetics.erase(d_kinetics.begin()+n,d_kinetics.end());

		vector<Cantera::IdealGasPhase*>::iterator iter2;
		for (iter2=d_thermo.begin()+n; iter2!=d_thermo.end(); iter2++) {
			delete *iter2;
		}
		d_thermo.erase(d_thermo.begin()+n,d_thermo.end());
	}
}

void gasArray::init(void) 
{
	rootXmlNode = Cantera::get_XML_File(mechanismFile);
	phaseXmlNode = rootXmlNode->findNameID("phase",phaseID);

	Cantera::importPhase(*phaseXmlNode, &thermoBase);
}


Cantera::IdealGasPhase& gasArray::operator[](unsigned int i) const
{
	return *d_thermo[i];
}

Cantera::IdealGasPhase& gasArray::thermo(unsigned int i) const
{
	return *d_thermo[i];
}

Cantera::GasKinetics& gasArray::kinetics(unsigned int i) const
{
	return *d_kinetics[i];
}
