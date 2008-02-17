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

class gasArray
{
public:
	gasArray();
	~gasArray();
	std::string mechanismFile;
	std::string phaseID;
	
	void resize(unsigned int n);
	void init(void);
	Cantera::IdealGasPhase& operator[](unsigned int i) const;
	Cantera::IdealGasPhase& thermo(unsigned int i) const;
	Cantera::GasKinetics& kinetics(unsigned int i) const;

private:
	Cantera::XML_Node* rootXmlNode;
	Cantera::XML_Node* root;
	Cantera::XML_Node* phaseXmlNode;
	Cantera::IdealGasPhase thermoBase;

	vector<Cantera::IdealGasPhase*> d_thermo;
	vector<Cantera::GasKinetics*> d_kinetics;
};