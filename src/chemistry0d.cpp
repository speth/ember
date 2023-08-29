#include "chemistry0d.h"
#include "debugUtils.h"
#include "readConfig.h"

#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/transport/MixTransport.h"
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/UnityLewisTransport.h"
#include "cantera/transport/TransportFactory.h"

ApproxMixTransport::ApproxMixTransport()
    : _threshold(0.0)
{
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
    assert(m_viscmix > 0 && m_viscmix < 1e100);
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
                for (size_t j : _kMajor) {
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
                for (size_t i : _kMajor) {
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
                for (size_t j : _kMajor) {
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
            m_phi(k,j) = factor1*factor1 / (std::sqrt(8.0) * m_wratkj1(j,k));
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
    m_viscwt_ok = false;
}


InterpKinetics::InterpKinetics()
{
    //! Preempt the creation of the standard MultiRate evaluator for Arrhenius reactions
    m_bulk_types["Arrhenius"] = m_bulk_rates.size();
    m_bulk_rates.push_back(std::make_unique<MultiArrheniusInterp>());
}

void MultiArrheniusInterp::rebuildInterpData(const ThermoPhase& thermo, const Kinetics& kin)
{
    // We need to modify the ThermoPhase state, but we will put it back the way
    // we found it.
    ThermoPhase& mthermo = const_cast<ThermoPhase&>(thermo);
    double Tnow = thermo.temperature();
    if (Tnow >= m_Tmax) {
        m_nTemps += 10;
        m_Tmax += 100.0;
    } else if (Tnow <= m_Tmin) {
        m_nTemps += 2;
        m_Tmin = std::max(Tnow - 20.0, 10.0);
    }

    m_kf_const = dmatrix::Zero(m_nReactions, m_nTemps);
    m_kf_slope = dmatrix::Zero(m_nReactions, m_nTemps);

    double dens_save = thermo.density();
    vector<double> kf_global(kin.nReactions());
    for (size_t n = 0; n < m_nTemps; n++) {
        mthermo.setState_TD(m_Tmin + ((double) n)/(m_nTemps - 1) * (m_Tmax - m_Tmin),
                            dens_save);

        m_base.update(thermo, kin);
        m_base.getRateConstants(kf_global.data());
        for (const auto& [jGlobal, jLocal] : m_indices) {
            m_kf_const(jLocal, n) = kf_global[jGlobal];
        }
    }

    double dT = (m_Tmax - m_Tmin) / (m_nTemps - 1);
    for (size_t n = 0; n < m_nTemps - 1; n++) {
        m_kf_slope.col(n) = (m_kf_const.col(n+1) - m_kf_const.col(n)) / dT;
    }
    mthermo.setState_TD(Tnow, dens_save);
}


CanteraGas::CanteraGas()
    : isInitialized(false)
    , lastRateMultiplier(1.0)
{
}

void CanteraGas::setOptions(const ConfigOptions& options)
{
    transportModel = options.transportModel;
    kineticsModel = options.kineticsModel;
    transportThreshold = options.transportThreshold;
    mechanismFile = options.gasMechanismFile;
    phaseID = options.gasPhaseID;
    pressure = options.pressure;
}

void CanteraGas::initialize()
{
    // suppress thermo warnings to simplify output. Without this, warnings will
    // be printed for every thread created
    Cantera::suppress_thermo_warnings();

    auto kin_fac = Cantera::KineticsFactory::factory();
    if (kineticsModel == "interp") {
        // Overwrite the normal factory function
        kin_fac->reg("gas", []() { return new InterpKinetics(); });
    } else {
        kin_fac->reg("gas", []() { return new Cantera::BulkKinetics(); });
    }

    auto tran_fac = Cantera::TransportFactory::factory();
    tran_fac->reg("Approx", []() { return new ApproxMixTransport(); });

    soln = Cantera::newSolution(mechanismFile, phaseID, transportModel);
    thermo = soln->thermo();
    kinetics = soln->kinetics();
    transport = soln->transport();

    if (transportModel == "Approx") {
        auto& atran = dynamic_cast<ApproxMixTransport&>(*transport);
        atran.setThreshold(transportThreshold);
    }

    nSpec = thermo->nSpecies();
    Dbin.setZero(nSpec,nSpec);
    isInitialized = true;
}

void CanteraGas::setStateMass(const dvec& Y, const double T)
{
    thermo->setState_TPY(T, pressure, Y.data());
}

void CanteraGas::setStateMass(const double* Y, const double T)
{
    thermo->setState_TPY(T, pressure, Y);
}

void CanteraGas::setStateMole(const dvec& X, const double T)
{
    thermo->setState_TPX(T, pressure, X.data());
}

void CanteraGas::setStateMole(const double* X, const double T)
{
    thermo->setState_TPX(T, pressure, X);
}

void CanteraGas::getMoleFractions(dvec& X) const
{
    thermo->getMoleFractions(X.data());
}

void CanteraGas::getMoleFractions(double* X) const
{
    thermo->getMoleFractions(X);
}

void CanteraGas::getMassFractions(double* Y) const
{
    thermo->getMassFractions(Y);
}

void CanteraGas::getMassFractions(dvec& Y) const
{
    thermo->getMassFractions(Y.data());
}

double CanteraGas::getDensity() const
{
    return thermo->density();
}

double CanteraGas::getMixtureMolecularWeight() const
{
    return thermo->meanMolecularWeight();
}

void CanteraGas::getMolecularWeights(dvec& W) const
{
    getMolecularWeights(W.data());
}

void CanteraGas::getMolecularWeights(double* W) const
{
    for (size_t k=0; k<nSpec; k++) {
        W[k] = thermo->molecularWeight(k);
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
    double rho = thermo->density();
    for (size_t k=0; k<nSpec; k++) {
        rhoD[k] *= rho;
    }
}

void CanteraGas::getWeightedDiffusionCoefficientsMass(double* rhoD)
{
    double rho = thermo->density();

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
    return thermo->cp_mass();
}

void CanteraGas::getSpecificHeatCapacities(dvec& cpSpec) const
{
    thermo->getPartialMolarCp(cpSpec.data());
}

void CanteraGas::getSpecificHeatCapacities(double* cpSpec) const
{
    thermo->getPartialMolarCp(cpSpec);
}


void CanteraGas::getEnthalpies(dvec& hk) const
{
    thermo->getPartialMolarEnthalpies(hk.data());
}

void CanteraGas::getEnthalpies(double* hk) const
{
    thermo->getPartialMolarEnthalpies(&hk[0]);
}

void CanteraGas::setRateMultiplier(double m)
{
    if (m == lastRateMultiplier) {
        return;
    }
    for (size_t k = 0; k < kinetics->nReactions(); k++) {
        kinetics->setMultiplier(k, m);
    }
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
