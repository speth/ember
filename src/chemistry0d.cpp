#include "chemistry0d.h"
#include "boost/filesystem.hpp"
#include "debugUtils.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

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

    if (m_visc_ok) {
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

    m_visc_ok = true;
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

    double mmw = m_thermo->meanMolecularWeight();
    double sumxw = 0.0, sum2;
    double p = m_thermo->pressure();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            sumxw += m_molefracs[k] * m_mw[k];
        }
        for (size_t k = 0; k < m_nsp; k++) {
            sum2 = 0.0;
            if (m_molefracs[k] >= _threshold) {
                for (size_t j = 0; j < m_nsp; j++) {
                    if (j != k) {
                        sum2 += m_molefracs[j] / m_bdiff(j,k);
                    }
                }
            } else {
                foreach (size_t j, _kMajor) {
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
    double p = m_thermo->pressure();

    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (size_t k=0; k<m_nsp; k++) {
            double sum1 = 0;
            double sum2 = 0;
            if (m_molefracs[k] >= _threshold) {
                for (size_t i=0; i<m_nsp; i++) {
                    if (i==k) {
                        continue;
                    }
                    sum1 += m_molefracs[i] / m_bdiff(k,i);
                    sum2 += m_molefracs[i] * m_mw[i] / m_bdiff(k,i);
                }
            } else {
                foreach (size_t i, _kMajor) {
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
    doublereal p = m_thermo->pressure();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            sum2 = 0.0;
            if (m_molefracs[k] > _threshold) {
                for (size_t j = 0; j < m_nsp; j++) {
                    if (j != k) {
                        sum2 += m_molefracs[j] / m_bdiff(j,k);
                    }
                }
            } else {
                foreach (size_t j, _kMajor) {
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
            size_t j = _kMajor[ij];
            for (size_t k=0; k<m_nsp; k++) {
                size_t ic = (k > j) ? m_nsp*j - (j-1)*j/2 + k - j
                                    : m_nsp*k - (k-1)*k/2 + j - k;
                m_bdiff(j,k) = exp(Cantera::dot4(m_polytempvec, m_diffcoeffs[ic]));
                m_bdiff(k,j) = m_bdiff(j,k);
            }
        }
    } else {
        for (size_t ij=0; ij<_kMajor.size(); ij++) {
            size_t j = _kMajor[ij];
            for (size_t k=0; k<m_nsp; k++) {
                size_t ic = (k > j) ? m_nsp*j - (j-1)*j/2 + k - j
                                    : m_nsp*k - (k-1)*k/2 + j - k;
                m_bdiff(j,k) = m_temp * m_sqrt_t*Cantera::dot5(m_polytempvec, m_diffcoeffs[ic]);
                m_bdiff(k,j) = m_bdiff(j,k);
            }
        }
    }

    m_bindiff_ok = true;
}

void ApproxMixTransport::update_C()
{
    MixTransport::update_C();

    _kMajor.clear();
    for (size_t k=0; k<m_nsp; k++) {
        if (m_molefracs[k] >= _threshold) {
            _kMajor.push_back(k);
        }
    }
}

CanteraGas::CanteraGas()
    : isInitialized(false)
    , rootXmlNode(NULL)
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
        transport = dynamic_cast<Cantera::GasTransport*>(
            transFac->newTransport("Multi",&thermo));
    } else if (transportModel == "Mix") {
        transport = dynamic_cast<Cantera::GasTransport*>(
            transFac->newTransport("Mix",&thermo));
    } else if (transportModel == "Approx") {
        ApproxMixTransport* atran = new ApproxMixTransport(thermo, *transFac);
        atran->setThreshold(transportThreshold);
        transport = atran;
    } else {
        throw debugException("Error: Invalid transport model specified.");
    }
    transFac->initTransport(transport, &thermo);
    transFac->deleteFactory();

    nSpec = thermo.nSpecies();
    Dbin.setZero(nSpec,nSpec);
    isInitialized = true;
}

dvec CanteraGas::calculateReactantMixture(const std::string& fuel,
                                             const std::string& oxidizer,
                                             double equivalenceRatio)
{
    int mC = thermo.elementIndex("C");
    int mO = thermo.elementIndex("O");
    int mH = thermo.elementIndex("H");

    double Cf(0), Hf(0), Of(0); // moles of C/H/O in fuel
    double Co(0), Ho(0), Oo(0); // moles of C/H/O in oxidizer

    dvec Xf(nSpec), Xo(nSpec), Xr(nSpec);

    thermo.setState_TPX(300.0, pressure, fuel);
    getMoleFractions(Xf);
    thermo.setState_TPX(300.0, pressure, oxidizer);
    getMoleFractions(Xo);

    dvec a(thermo.nElements());
    for (size_t k=0; k<nSpec; k++) {
        thermo.getAtoms(k, &a[0]);
        if (mC != -1) {
            Cf += a[mC]*Xf[k];
            Co += a[mC]*Xo[k];
        }
        if (mH != -1) {
            Hf += a[mH]*Xf[k];
            Ho += a[mH]*Xo[k];
        }
        if (mO != -1) {
            Of += a[mO]*Xf[k];
            Oo += a[mO]*Xo[k];
        }
    }
    double stoichAirFuelRatio = -(Of-2*Cf-Hf/2)/(Oo-2*Co-Ho/2);
    Xr = Xf * equivalenceRatio + stoichAirFuelRatio * Xo;
    Xr /= Xr.sum();

    return Xr;
}

void CanteraGas::setStateMass(const dvec& Y, const double T)
{
    thermo.setState_TPY(T, pressure, Y.data());
}

void CanteraGas::setStateMass(const double* Y, const double T)
{
    thermo.setState_TPY(T, pressure, Y);
}

void CanteraGas::setStateMole(const dvec& X, const double T)
{
    thermo.setState_TPX(T, pressure, X.data());
}

void CanteraGas::setStateMole(const double* X, const double T)
{
    thermo.setState_TPX(T, pressure, X);
}

void CanteraGas::getMoleFractions(dvec& X) const
{
    thermo.getMoleFractions(X.data());
}

void CanteraGas::getMoleFractions(double* X) const
{
    thermo.getMoleFractions(X);
}

void CanteraGas::getMassFractions(double* Y) const
{
    thermo.getMassFractions(Y);
}

void CanteraGas::getMassFractions(dvec& Y) const
{
    thermo.getMassFractions(Y.data());
}

double CanteraGas::getDensity() const
{
    return thermo.density();
}

double CanteraGas::getMixtureMolecularWeight() const
{
    return thermo.meanMolecularWeight();
}

void CanteraGas::getMolecularWeights(dvec& W) const
{
    getMolecularWeights(W.data());
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

void CanteraGas::getDiffusionCoefficientsMole(dvec& Dkm) const
{
    transport->getMixDiffCoeffs(Dkm.data());
}

void CanteraGas::getDiffusionCoefficientsMole(double* Dkm) const
{
    transport->getMixDiffCoeffs(Dkm);
}

void CanteraGas::getWeightedDiffusionCoefficientsMole(dvec& rhoD) const
{
    getWeightedDiffusionCoefficientsMole(rhoD.data());
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

    transport->getMixDiffCoeffsMass(rhoD);
    for (size_t k=0; k<nSpec; k++) {
        rhoD[k] *= rho;
    }
}

void CanteraGas::getWeightedDiffusionCoefficientsMass(dvec& rhoD)
{
    getWeightedDiffusionCoefficientsMass(rhoD.data());
}

void CanteraGas::getThermalDiffusionCoefficients(dvec& Dkt) const
{
    transport->getThermalDiffCoeffs(Dkt.data());
}

void CanteraGas::getThermalDiffusionCoefficients(double* Dkt) const
{
    transport->getThermalDiffCoeffs(Dkt);
}

double CanteraGas::getSpecificHeatCapacity() const
{
    return thermo.cp_mass();
}

void CanteraGas::getSpecificHeatCapacities(dvec& cpSpec) const
{
    thermo.getPartialMolarCp(cpSpec.data());
}

void CanteraGas::getSpecificHeatCapacities(double* cpSpec) const
{
    thermo.getPartialMolarCp(cpSpec);
}


void CanteraGas::getEnthalpies(dvec& hk) const
{
    thermo.getPartialMolarEnthalpies(hk.data());
}

void CanteraGas::getEnthalpies(double* hk) const
{
    thermo.getPartialMolarEnthalpies(&hk[0]);
}

void CanteraGas::getReactionRates(dvec& wDot) const
{
    kinetics->getNetProductionRates(wDot.data());
}

void CanteraGas::getReactionRates(double* wDot) const
{
    kinetics->getNetProductionRates(wDot);
}

void CanteraGas::getCreationRates(dvec& wDot) const
{
    kinetics->getCreationRates(wDot.data());
}

void CanteraGas::getCreationRates(double* wDot) const
{
    kinetics->getCreationRates(wDot);
}

void CanteraGas::getDestructionRates(dvec& wDot) const
{
    kinetics->getDestructionRates(wDot.data());
}

void CanteraGas::getDestructionRates(double* wDot) const
{
    kinetics->getDestructionRates(wDot);
}
