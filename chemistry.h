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

// NOTE: Use of this class requires a slightly modified installation of Cantera
// gasArray modifies protected members of Transport objects, requiring gasArray
// to be declared as a friend of Cantera's transport-related classes.
// In the files "transport/TransportBase.h", "transport/MixTransport.h" and
// "transport/MultiTransport.h", add the following declaration near the top of
// the file (before the start of namespace Cantera):
//
// class gasArray;
//
// And, as the first line of class Transport, class MixTransport and 
// class MultiTransport, add the friend declaration (before "public:")
// 
// friend class ::gasArray;
//

class gasArray
{
public:
	gasArray();
	~gasArray();
	std::string mechanismFile;
	std::string phaseID;
	double pressure; // thermodynamic pressure
	
	void initialize(void);
	void resize(unsigned int n);
	
	void setState(Cantera::Array2D& Y, dvector& T);

	void getMoleFractions(Cantera::Array2D& X);
	void getDensity(dvector& rho);
	void getMixtureMolecularWeight(dvector& Wmx);
	void getMolecularWeights(dvector& W);

	void getViscosity(dvector& mu);
	void getThermalConductivity(dvector& lambda);
	void getDiffusionCoefficients(Cantera::Array2D& Dkm);
	void getThermalDiffusionCoefficients(Cantera::Array2D& Dkt);
	
	void getSpecificHeatCapacity(dvector& cp);
	void getSpecificHeatCapacities(Cantera::Array2D& cpSpec);
	void getEnthalpies(Cantera::Array2D& hk);
	
	void getReactionRates(Cantera::Array2D& wDot);

	Cantera::IdealGasPhase& operator[](unsigned int i) const;
	Cantera::IdealGasPhase& thermo(unsigned int i) const;
	Cantera::GasKinetics& kinetics(unsigned int i) const;
	Cantera::MultiTransport& trans(unsigned int i) const;

	void testFunction(void);

private:
	Cantera::XML_Node* rootXmlNode;
	Cantera::XML_Node* phaseXmlNode;

	int nPoints;
	int nSpec;

	vector<Cantera::IdealGasPhase*> m_thermo;
	vector<Cantera::GasKinetics*> m_kinetics;
	vector<Cantera::MultiTransport*> m_transport;

	// Default objects
	Cantera::IdealGasPhase m_thermoBase;
	Cantera::GasKinetics* m_kineticsBase;
	Cantera::MultiTransport* m_transportBase;
};
