#include "sourceSystem.h"
#include "readConfig.h"

#include <boost/format.hpp>

SourceSystem::SourceSystem()
    : U(NaN)
    , T(NaN)
    , updateDiagonalJac(false)
    , gas(NULL)
    , quasi2d(false)
{
}

void SourceSystem::resize(size_t new_nSpec)
{
    nSpec = new_nSpec;
    Y.resize(nSpec, NaN);
    dYdt.resize(nSpec);
    cpSpec.resize(nSpec);
    wDot.resize(nSpec);
    hk.resize(nSpec);
    splitConst.resize(nSpec+2);
    diagonalJac.resize(nSpec+2, 0);
}

void SourceSystem::resetSplitConstants()
{
    splitConst.assign(splitConst.size(), 0);
}

void SourceSystem::setupQuasi2d(boost::shared_ptr<BilinearInterpolator> vzInterp,
                                boost::shared_ptr<BilinearInterpolator> TInterp)
{
    quasi2d = true;
    vzInterp_ = vzInterp;
    TInterp_ = TInterp;
}

int SourceSystem::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y, t);

    // *** Update auxiliary data ***
    reactionRatesTimer->start();
    gas->setStateMass(Y, T);
    gas->getReactionRates(wDot);
    reactionRatesTimer->stop();

    thermoTimer->start();
    gas->getEnthalpies(hk);
    rho = gas->getDensity();
    cp = gas->getSpecificHeatCapacity();
    thermoTimer->stop();

    qDot = 0.0;
    for (size_t k=0; k<nSpec; k++) {
        qDot -= wDot[k]*hk[k];
    }

    if (t >= options->ignition_tStart && t < options->ignition_tStart + options->ignition_duration) {
        qDot += options->ignition_energy /
                (options->ignition_stddev * sqrt(2 * M_PI) * options->ignition_duration) *
                exp(-pow(x - options->ignition_center, 2) / (2 * pow(options->ignition_stddev, 2)));
    }

    double a = strainFunction.a(t);
    double dadt = strainFunction.dadt(t);

    // *** Calculate the time derivatives
    if (!quasi2d) {
        dUdt = - U*U + rhou/rho*(dadt + a*a) + splitConst[kMomentum];
        dTdt = qDot/(rho*cp) + splitConst[kEnergy];
    } else {
        dUdt = splitConst[kMomentum];
        dTdt = splitConst[kEnergy];
    }
    double scale = (quasi2d) ? 1.0/vzInterp_->get(x, t) : 1.0;
    for (size_t k=0; k<nSpec; k++) {
        dYdt[k] = scale*wDot[k]*W[k]/rho + splitConst[kSpecies+k];
    }

    roll_ydot(ydot);
    return 0;
}

int SourceSystem::denseJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdMatrix& J)
{
    // TODO: Verify that f has just been called so that we don't need to
    // unroll y and compute all the transport properties.

    fdJacobian(t, y, ydot, J);
    return 0;

    double a = strainFunction.a(t);
    double dadt = strainFunction.dadt(t);
    double A = a*a + dadt;

    // Additional properties not needed for the normal function evaluations:
    thermoTimer->start();
    gas->getSpecificHeatCapacities(cpSpec);
    Wmx = gas->getMixtureMolecularWeight();
    thermoTimer->stop();

    // The constant "800" here has been empirically determined to give
    // good performance for typical test cases. This value can have
    // a substantial impact on the convergence rate of the solver.
    double eps = sqrt(DBL_EPSILON)*800;

    // *** Derivatives with respect to temperature
    double TplusdT = T*(1+eps);

    double dwhdT = 0;
    dvector dwdT(nSpec);
    dvector wDot2(nSpec);

    reactionRatesTimer->start();
    gas->setStateMass(Y, TplusdT);
    gas->getReactionRates(wDot2);
    reactionRatesTimer->stop();

    for (size_t k=0; k<nSpec; k++) {
        dwdT[k] = (wDot2[k]-wDot[k])/(TplusdT-T);
        dwhdT += hk[k]*dwdT[k] + cpSpec[k]*wDot[k];
    }

    double drhodT = -rho/T;

    // *** Derivatives with respect to species concentration
    dmatrix dwdY(nSpec, nSpec);
    dvector hdwdY(nSpec, 0);
    dvector YplusdY(nSpec);
    dvector drhodY(nSpec);

    double scale = (quasi2d) ? 1.0/vzInterp_->get(x, t) : 1.0;

    for (size_t k=0; k<nSpec; k++) {
        YplusdY[k] = (abs(Y[k]) > eps/2) ? Y[k]*(1+eps) : eps;
        reactionRatesTimer->start();
        gas->setStateMass(YplusdY, T);
        gas->getReactionRates(wDot2);
        reactionRatesTimer->stop();

        for (size_t i=0; i<nSpec; i++) {
            dwdY(i,k) = (wDot2[i]-wDot[i])/(YplusdY[k]-Y[k]);
            hdwdY[k] += hk[i]*dwdY(i,k);
        }

        drhodY[k] = rho*(W[k]-Wmx)/(W[k]*(1-Y[k]*(1-eps)));
    }

    for (size_t k=0; k<nSpec; k++) {
        for (size_t i=0; i<nSpec; i++) {
            // dSpecies/dY
            J(kSpecies+k, kSpecies+i) = scale * (dwdY(k,i)*W[k]/rho - wDot[k]*W[k]*drhodY[i]/(rho*rho));
        }
        if (!quasi2d) {
            // dSpecies/dT
            J(kSpecies+k, kEnergy) = dwdT[k]*W[k]/rho - wDot[k]*W[k]*drhodT/(rho*rho);

            // dEnergy/dY
            J(kEnergy, kSpecies+k) = -hdwdY[k]/(rho*cp) - qDot*drhodY[k]/(rho*rho*cp);

            // dMomentum/dY
            J(kMomentum, kSpecies+k) = -A*drhodY[k]/(rho*rho);
        }
    }

    if (!quasi2d) {
        // dEnergy/dT
        J(kEnergy, kEnergy) = -dwhdT/(rho*cp) - qDot*drhodT/(rho*rho*cp);

        // dMomentum/dU
        J(kMomentum, kMomentum) = -2*U;

        // dMomentum/dT
        J(kMomentum, kEnergy) = -A*drhodT/(rho*rho);
    }

    return 0;
}

int SourceSystem::fdJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdMatrix& J)
{
    jacobianTimer->start();
    sdVector yplusdy(y.length());
    sdVector ydot2(y.length());
    size_t nVars = nSpec+2;
    double eps = sqrt(DBL_EPSILON);
    double atol = DBL_EPSILON;

    for (size_t i=0; i<nVars; i++) {
        for (size_t k=0; k<nVars; k++) {
            yplusdy[k] = y[k];
        }
        double dy = (abs(y[i]) > atol) ? abs(y[i])*(eps) : abs(y[i])*eps + atol;
        yplusdy[i] += dy;
        f(t, yplusdy, ydot2);
        for (size_t k=0; k<nVars; k++) {
            J(k,i) = (ydot2[k]-ydot[k])/dy;
        }
    }

    if (updateDiagonalJac) {
        for (size_t k=0; k<nSpec+2; k++) {
            diagonalJac[k] = J(k,k);
        }
        updateDiagonalJac = false;
    }
    jacobianTimer->stop();

    return 0;
}

void SourceSystem::unroll_y(const sdVector& y, double t)
{
    if (!quasi2d) {
        T = y[kEnergy];
        U = y[kMomentum];
    } else {
        T = TInterp_->get(x, t);
        U = 0;
    }
    for (size_t k=0; k<nSpec; k++) {
        Y[k] = y[kSpecies+k];
    }
}

void SourceSystem::roll_y(sdVector& y) const
{
    y[kEnergy] = T;
    y[kMomentum] = U;
    for (size_t k=0; k<nSpec; k++) {
        y[kSpecies+k] = Y[k];
    }
}

void SourceSystem::roll_ydot(sdVector& ydot) const
{
    ydot[kEnergy] = dTdt;
    ydot[kMomentum] = dUdt;
    for (size_t k=0; k<nSpec; k++) {
        ydot[kSpecies+k] = dYdt[k];
    }
}

void SourceSystem::writeJacobian(sundialsCVODE& solver, std::ostream& out)
{
    size_t N = solver.y.length();
    double t = solver.tInt;
    sdMatrix J(N,N);
    sdVector ydot(N);
    f(t, solver.y, ydot);
    denseJacobian(t, solver.y, ydot, J);

    out << "J = []" << std::endl;
    for (size_t i=0; i<N; i++) {
        out << "J.append([";
        for (size_t k=0; k<N; k++) {
            out << boost::format("%.5e, ") % J(i,k);
        }
        out << "])" << std::endl;
    }
}

void SourceSystem::writeState(sundialsCVODE& solver, std::ostream& out, bool init)
{
    if (init) {
        out << "T = []" << std::endl;
        out << "Y = []" << std::endl;
        out << "t = []" << std::endl;
        out << "wDot = []" << std::endl;
    }

    out << "T.append(" << T << ")" << std::endl;
    out << "Y.append(" << Y << ")" << std::endl;
    out << "t.append(" << solver.tInt << ")" << std::endl;
    out << "wDot.append(" << wDot << ")" << std::endl;
}

// ********************


SourceSystemQSS::SourceSystemQSS()
    : U(NaN)
    , T(NaN)
    , gas(NULL)
{
    dUdtQ = 0;
    dUdtD = 0;
    dTdtQ = 0;
    dTdtD = 0;
}

void SourceSystemQSS::initialize(size_t new_nSpec)
{
    QSSIntegrator::initialize(new_nSpec + 2);
    nSpec = new_nSpec;
    Y.resize(nSpec, NaN);
    dYdtQ.resize(nSpec, 0);
    dYdtD.resize(nSpec, 0);
    cpSpec.resize(nSpec);
    splitConstY.resize(nSpec);
    wDotD.resize(nSpec);
    wDotQ.resize(nSpec);
    hk.resize(nSpec);

    enforce_ymin[kMomentum] = false;
}

void SourceSystemQSS::setOptions(configOptions& options_)
{
    options = &options_;
    epsmin = options->qss_epsmin;
    epsmax = options->qss_epsmax;
    dtmin = options->qss_dtmin;
    dtmax = options->qss_dtmax;
    itermax = options->qss_iterationCount;
    abstol = options->qss_abstol;
    stabilityCheck = options->qss_stabilityCheck;
}

void SourceSystemQSS::resetSplitConstants()
{
    splitConstY.assign(nSpec, 0);
    splitConstT = 0;
    splitConstU = 0;
}

void SourceSystemQSS::odefun(double t, const dvector& y, dvector& q, dvector& d, bool corrector)
{
    unroll_y(y, corrector);

    // *** Update auxiliary data ***
    reactionRatesTimer->start();
    gas->setStateMass(Y, T);
    gas->getCreationRates(wDotQ);
    gas->getDestructionRates(wDotD);
    reactionRatesTimer->stop();

    if (!corrector) {
        thermoTimer->start();
        gas->getEnthalpies(hk);
        rho = gas->getDensity();
        cp = gas->getSpecificHeatCapacity();
        thermoTimer->stop();
    }

    qDot = 0.0;
    for (size_t k=0; k<nSpec; k++) {
        qDot -= (wDotQ[k]-wDotD[k])*hk[k];
    }

    if (t >= options->ignition_tStart && t < options->ignition_tStart + options->ignition_duration) {
        qDot += options->ignition_energy /
                (options->ignition_stddev * sqrt(2 * M_PI) * options->ignition_duration) *
                exp(-pow(x - options->ignition_center, 2) / (2 * pow(options->ignition_stddev, 2)));
    }

    double a = strainFunction.a(t);
    double dadt = strainFunction.dadt(t);

    // *** Calculate the time derivatives
    dUdtQ = rhou/rho*(dadt + a*a) - U*U + splitConstU;
    dUdtD = 0;
    dTdtQ = qDot/(rho*cp) + splitConstT;
    for (size_t k=0; k<nSpec; k++) {
        dYdtQ[k] = wDotQ[k] * W[k] / rho + splitConstY[k];
        dYdtD[k] = wDotD[k] * W[k] / rho;
    }

    assert(rhou > 0);
    assert(rho > 0);
    assert(dadt > -1e100 && dadt < 1e100);
    assert(a > -1e100 && a < 1e100);
    assert(U > -1e100 && U < 1e100);
    assert(splitConstU > -1e100 && splitConstU < 1e100);

    assert(dUdtQ > -1e100 && dUdtQ < 1e100);
    assert(dUdtD > -1e100 && dUdtD < 1e100);

    roll_ydot(q, d);
}

void SourceSystemQSS::unroll_y(const dvector& y, bool corrector)
{
    if (!corrector) {
        T = y[kEnergy];
    }
    U = y[kMomentum];
    for (size_t k=0; k<nSpec; k++) {
        Y[k] = y[kSpecies+k];
    }
}

void SourceSystemQSS::roll_y(dvector& y) const
{
    y[kEnergy] = T;
    y[kMomentum] = U;
    for (size_t k=0; k<nSpec; k++) {
        y[kSpecies+k] = Y[k];
    }
}

void SourceSystemQSS::roll_ydot(dvector& q, dvector& d) const
{
    q[kEnergy] = dTdtQ;
    d[kEnergy] = dTdtD;
    q[kMomentum] = dUdtQ;
    d[kMomentum] = dUdtD;
    for (size_t k=0; k<nSpec; k++) {
        q[kSpecies+k] = dYdtQ[k];
        d[kSpecies+k] = dYdtD[k];
    }
}

void SourceSystemQSS::writeState(std::ostream& out, bool init)
{
    if (init) {
        out << "T = []" << std::endl;
        out << "Y = []" << std::endl;
        out << "t = []" << std::endl;
    }

    out << "T.append(" << T << ")" << std::endl;
    out << "Y.append(" << Y << ")" << std::endl;
    out << "t.append(" << tn << ")" << std::endl;
}


