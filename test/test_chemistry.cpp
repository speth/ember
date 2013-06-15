#include "../src/chemistry0d.h"
#include "gtest/gtest.h"

#ifdef EMBER_EXTENDED_MULTITRANSPORT

class EigenTransportTest : public ::testing::Test
{
public:
    EigenTransportTest() {
        multi_gas.reset(Cantera::newPhase("gri30.xml", ""));
        eigen_gas.reset(Cantera::newPhase("gri30.xml", ""));
        multi.reset(dynamic_cast<Cantera::MultiTransport*>(Cantera::newTransportMgr("Multi", multi_gas.get())));
        eigen.reset(new MultiTransportEigen(*eigen_gas, *Cantera::TransportFactory::factory()));
        K = multi_gas->nSpecies();
    }

    void check_thermalDiffCoeffs() {
        dvec Dkt_multi = dvec::Constant(K, 0.0);
        dvec Dkt_eigen = dvec::Constant(K, 0.0);
        multi->getThermalDiffCoeffs(&Dkt_multi[0]);
        eigen->getThermalDiffCoeffs(&Dkt_eigen[0]);
        for (size_t j = 0; j < K; j++) {
            ASSERT_NEAR(Dkt_multi[j], Dkt_eigen[j], 1e-12);
        }
    }

    std::auto_ptr<Cantera::ThermoPhase> multi_gas, eigen_gas;
    std::auto_ptr<Cantera::MultiTransport> multi;
    std::auto_ptr<MultiTransportEigen> eigen;
    size_t K;
};


TEST_F(EigenTransportTest, ThermalConductivityPure)
{
    for (size_t k = 0; k < K; k++) {
        vector<double> X(K, 0.0);
        X[k] = 1.0;
        for (int i = 0; i < 2; i++) {
            multi_gas->setState_TPX(300 + 600 * i, 101325 + 50000 * i, &X[0]);
            eigen_gas->setState_TPX(300 + 600 * i, 101325 + 50000 * i, &X[0]);
            ASSERT_NEAR(multi->thermalConductivity(),
                        eigen->thermalConductivity(), 1e-11);
        }
    }
}

TEST_F(EigenTransportTest, ArbitraryMixture)
{
    for (size_t k = 0; k < K; k++) {
        // Set up an arbitrary mixture composition
        vector<double> X(K, 0.0);
        for (size_t j = 0; j < 10; j++) {
            X[(k+j) % K] = (j + 1) / 54.0;
        }
        multi_gas->setState_TPX(300, 101325, &X[0]);
        Cantera::equilibrate(*multi_gas, "TP");
        multi_gas->getMoleFractions(&X[0]);
        eigen_gas->setMoleFractions(&X[0]);

        for (int i = 0; i < 2; i++) {
            multi_gas->setState_TP(500 + 600 * i, 101325 + 50000 * i);
            eigen_gas->setState_TP(500 + 600 * i, 101325 + 50000 * i);
            ASSERT_NEAR(multi->thermalConductivity(),
                        eigen->thermalConductivity(), 1e-11);
        check_thermalDiffCoeffs();
        }
    }
}

TEST_F(EigenTransportTest, ThermalDiffusionCoeffsBinary)
{
    size_t kN2 = multi_gas->speciesIndex("N2");
    for (size_t k = 0; k < K; k++) {
        dvec X = dvec::Constant(K, 0.0);
        X[k] = 0.3;
        X[kN2] = 0.7;
        for (int i = 0; i < 2; i++) {
            multi_gas->setState_TPX(300 + 600 * i, 101325, &X[0]);
            eigen_gas->setState_TPX(300 + 600 * i, 101325, &X[0]);
            check_thermalDiffCoeffs();
        }
    }
}

#endif
