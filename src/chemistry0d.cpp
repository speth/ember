#include "chemistry0d.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#ifdef CANTERA_EXTENDED_TRANSPORT
ApproxMixTransport::ApproxMixTransport(Cantera::ThermoPhase& thermo,
                                       Cantera::TransportFactory& factory)
    : _threshold(0.0)
{
    dvector state;
    thermo.saveState(state);
    factory.initTransport(this, &thermo, 0, 0);
    thermo.restoreState(state);
}

void ApproxMixTransport::setThreshold(double threshold)
{
    _threshold = threshold;
}

double ApproxMixTransport::viscosity()
{
    update_T();
    update_C();

    if (m_viscmix_ok) {
        return m_viscmix;
    }

    doublereal vismix = 0.0;

    // update m_visc and m_phi if necessary
    if (!m_viscwt_ok) {
        updateViscosity_T();
    }

    for (size_t ik=0; ik<_kMajor.size(); ik++) {
        size_t k = _kMajor[ik];
        double sum = 0;
        for (size_t ij=0; ij<_kMajor.size(); ij++) {
            size_t j = _kMajor[ij];
            sum += m_molefracs[j]*m_phi(k,j);
        }
        vismix += m_molefracs[k]*m_visc[k] / sum;
    }

    m_viscmix = vismix;
    return vismix;
}

void ApproxMixTransport::getMixDiffCoeffs(double* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    int k, j;
    double mmw = m_thermo->meanMolecularWeight();
    double sumxw = 0.0, sum2;
    double p = pressure_ig();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (k = 0; k < m_nsp; k++) {
            sumxw += m_molefracs[k] * m_mw[k];
        }
        for (k = 0; k < m_nsp; k++) {
            sum2 = 0.0;
            if (m_molefracs[k] >= _threshold) {
                for (j = 0; j < m_nsp; j++) {
                    if (j != k) {
                        sum2 += m_molefracs[j] / m_bdiff(j,k);
                    }
                }
            } else {
                foreach (j, _kMajor) {
                    if (j != k) {
                        sum2 += m_molefracs[j] / m_bdiff(j,k);
                    }
                }
            }

            if (sum2 <= 0.0) {
                d[k] = m_bdiff(k,k) / p;
            } else {
                d[k] = (sumxw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
            }
        }
    }
}

void ApproxMixTransport::getMixDiffCoeffsMass(double* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    double mmw = m_thermo->meanMolecularWeight();
    double p = pressure_ig();

    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (int k=0; k<m_nsp; k++) {
            double sum1 = 0;
            double sum2 = 0;
            if (m_molefracs[k] >= _threshold) {
                for (int i=0; i<m_nsp; i++) {
                    if (i==k) {
                        continue;
                    }
                    sum1 += m_molefracs[i] / m_bdiff(k,i);
                    sum2 += m_molefracs[i] * m_mw[i] / m_bdiff(k,i);
                }
            } else {
                foreach (int i, _kMajor) {
                    if (i==k) {
                        continue;
                    }
                    sum1 += m_molefracs[i] / m_bdiff(k,i);
                    sum2 += m_molefracs[i] * m_mw[i] / m_bdiff(k,i);
                }
            }
            sum1 *= p;
            sum2 *= p * m_molefracs[k] / (mmw - m_mw[k]*m_molefracs[k]);
            d[k] = 1.0 / (sum1 +  sum2);
        }
    }
}

void ApproxMixTransport::getMixDiffCoeffsMole(double* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    doublereal sum2;
    doublereal p = pressure_ig();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (int k = 0; k < m_nsp; k++) {
            sum2 = 0.0;
            if (m_molefracs[k] > _threshold) {
                for (int j = 0; j < m_nsp; j++) {
                    if (j != k) {
                        sum2 += m_molefracs[j] / m_bdiff(j,k);
                    }
                }
            } else {
                foreach (int j, _kMajor) {
                    if (j != k) {
                        sum2 += m_molefracs[j] / m_bdiff(j,k);
                    }
                }
            }

            if (sum2 <= 0.0) {
                d[k] = m_bdiff(k,k) / p;
            } else {
                d[k] = (1 - m_molefracs[k]) / (p * sum2);
            }
        }
    }
}

void ApproxMixTransport::updateViscosity_T()
{
    double vratiokj, wratiojk, factor1;

    if (!m_spvisc_ok) {
        updateSpeciesViscosities();
    }

    // see Eq. (9-5.15) of Reid, Prausnitz, and Poling
    for (size_t ij=0; ij<_kMajor.size(); ij++) {
        size_t j = _kMajor[ij];
        for (size_t ik=ij; ik<_kMajor.size(); ik++) {
            size_t k = _kMajor[ik];
            vratiokj = m_visc[k]/m_visc[j];
            wratiojk = m_mw[j]/m_mw[k];

            // Note that m_wratjk(k,j) holds the square root of
            // m_wratjk(j,k)!
            factor1 = 1.0 + (m_sqvisc[k]/m_sqvisc[j]) * m_wratjk(k,j);
            m_phi(k,j) = factor1*factor1 / (Cantera::SqrtEight * m_wratkj1(j,k));
            m_phi(j,k) = m_phi(k,j)/(vratiokj * wratiojk);
        }
    }
    m_viscwt_ok = true;
}

void ApproxMixTransport::updateDiff_T()
{
    // Evaluate binary diffusion coefficients at unit pressure for
    // the species pairs where at least one has a mole fraction above
    // the specified threshold.
    if (m_mode == Cantera::CK_Mode) {
        for (size_t ij=0; ij<_kMajor.size(); ij++) {
            int j = _kMajor[ij];
            for (int k=0; k<m_nsp; k++) {
                int ic = (k > j) ? m_nsp*j - (j-1)*j/2 + k - j
                                 : m_nsp*k - (k-1)*k/2 + j - k;
                m_bdiff(j,k) = exp(Cantera::dot4(m_polytempvec, m_diffcoeffs[ic]));
                m_bdiff(k,j) = m_bdiff(j,k);
            }
        }
    } else {
        for (size_t ij=0; ij<_kMajor.size(); ij++) {
            int j = _kMajor[ij];
            for (int k=0; k<m_nsp; k++) {
                int ic = (k > j) ? m_nsp*j - (j-1)*j/2 + k - j
                                 : m_nsp*k - (k-1)*k/2 + j - k;
                m_bdiff(j,k) = m_temp * m_sqrt_t*Cantera::dot5(m_polytempvec, m_diffcoeffs[ic]);
                m_bdiff(k,j) = m_bdiff(j,k);
            }
        }
    }

    m_bindiff_ok = true;
    m_diffmix_ok = false;
}

void ApproxMixTransport::update_C()
{
    MixTransport::update_C();

    _kMajor.clear();
    for (int k=0; k<m_nsp; k++) {
        if (m_molefracs[k] >= _threshold) {
            _kMajor.push_back(k);
        }
    }
}
#endif

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
    transportModel = options.transportModel;
    transportThreshold = options.transportThreshold;
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
    if (rootXmlNode == NULL) {
        throw debugException((format(
            "Error parsing Cantera XML file '%s'.") % mechanismFile).str());
    }

    phaseXmlNode = rootXmlNode->findNameID("phase", phaseID);
    if (phaseXmlNode == NULL) {
        throw debugException((format(
            "Error finding phase '%s' in '%s'.") % phaseID % mechanismFile).str());
    }

    // Initialize the default thermodynamic properties object
    Cantera::importPhase(*phaseXmlNode, &thermo);

    // Initialize the default chemical kinetics object
    kinetics = new Cantera::GasKinetics(&thermo);
    kinetics->init();
    Cantera::installReactionArrays(*phaseXmlNode, *kinetics, phaseID);
    kinetics->finalize();

    // Initialize the default transport properties object
    Cantera::TransportFactory* transFac = Cantera::TransportFactory::factory();
    if (transportModel == "Multi") {
        transport = transFac->newTransport("Multi",&thermo);
    } else if (transportModel == "Mix") {
        transport = transFac->newTransport("Mix",&thermo);
#ifdef CANTERA_EXTENDED_TRANSPORT
    } else if (transportModel == "Approx") {
        ApproxMixTransport* atran = new ApproxMixTransport(thermo, *transFac);
        atran->setThreshold(transportThreshold);
        transport = atran;
#endif
    } else {
        throw debugException("Error: Invalid transport model specified.");
    }
    transFac->initTransport(transport, &thermo);
    transFac->deleteFactory();

    nSpec = thermo.nSpecies();
    Dbin.resize(nSpec,nSpec);
    Dbin.setZero();
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
        Y[k] = std::max(Y_in[k], 0.0);
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
        X[k] = std::max(X_in[k], 0.0);
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
    double rho = thermo.density();

#ifdef CANTERA_EXTENDED_TRANSPORT
    transport->getMixDiffCoeffsMass(rhoD);
    for (size_t k=0; k<nSpec; k++) {
        rhoD[k] *= rho;
    }
#else
    thermo.getMassFractions(&Y[0]);
    thermo.getMoleFractions(&X[0]);
    transport->getBinaryDiffCoeffs(nSpec, &Dbin(0,0));
    double Keps = 1e-14;
    double eps = Keps/Y.size(); // Avoid 0/0 as Y[k] -> 1

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
#endif
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
