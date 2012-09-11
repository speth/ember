#include "sourceSystem.h"
#include "readConfig.h"

#include <boost/format.hpp>

SourceSystem::SourceSystem()
    : U(NaN)
    , T(NaN)
    , debug(false)
    , options(NULL)
    , gas(NULL)
    , quasi2d(false)
{
}

void SourceSystem::updateThermo()
{
    thermoTimer->start();
    gas->getEnthalpies(hk);
    rho = gas->getDensity();
    cp = gas->getSpecificHeatCapacity();
    thermoTimer->stop();
}

double SourceSystem::getQdotIgniter(double t)
{
    ConfigOptions& opt = *options;
    if (t >= opt.ignition_tStart &&
        t < opt.ignition_tStart + opt.ignition_duration) {
        return opt.ignition_energy /
                (opt.ignition_stddev * sqrt(2 * M_PI) * opt.ignition_duration) *
                exp(-pow(x - opt.ignition_center, 2) /
                    (2 * pow(opt.ignition_stddev, 2)));
    } else {
        return 0.0;
    }
}

void SourceSystem::initialize(size_t new_nSpec)
{
    nSpec = new_nSpec;

    Y.setConstant(nSpec, NaN);
    cpSpec.resize(nSpec);
    splitConst.resize(nSpec + 2);
    hk.resize(nSpec);

    W.resize(gas->nSpec); // move this to initialize
    gas->getMolecularWeights(W);
}

void SourceSystem::setTimers
(PerfTimer* reactionRates, PerfTimer* thermo, PerfTimer* jacobian)
{
    reactionRatesTimer = reactionRates;
    thermoTimer = thermo;
    jacobianTimer = jacobian;
}

void SourceSystem::setPosition(size_t _j, double _x)
{
    j = _j;
    x = _x;
}

void SourceSystem::setupQuasi2d(boost::shared_ptr<BilinearInterpolator> vz_int,
                                boost::shared_ptr<BilinearInterpolator> T_int)
{
    quasi2d = true;
    vzInterp = vz_int;
    TInterp = T_int;
}

void SourceSystem::writeState(std::ostream& out, bool init)
{
    if (init) {
        out << "T = []" << std::endl;
        out << "Y = []" << std::endl;
        out << "t = []" << std::endl;
    }

    out << "T.append(" << T << ")" << std::endl;
    out << "Y.append(" << Y << ")" << std::endl;
    out << "t.append(" << time() << ")" << std::endl;
}

// ----------------------------------------------------------------------------

void SourceSystemCVODE::initialize(size_t new_nSpec)
{
    SourceSystem::initialize(new_nSpec);
    dYdt.resize(nSpec);
    wDot.resize(nSpec);

    integrator.reset(new SundialsCvode(nSpec+2));
    integrator->setODE(this);
    integrator->linearMultistepMethod = CV_BDF;
    integrator->nonlinearSolverMethod = CV_NEWTON;
    integrator->maxNumSteps = 1000000;
}

void SourceSystemCVODE::setOptions(ConfigOptions& opts)
{
    SourceSystem::setOptions(opts);
    integrator->abstol[kMomentum] = options->integratorMomentumAbsTol;
    integrator->abstol[kEnergy] = options->integratorEnergyAbsTol;
    for (size_t k=0; k<nSpec; k++) {
        integrator->abstol[kSpecies+k] = options->integratorSpeciesAbsTol;
    }
    integrator->reltol = options->integratorRelTol;
    integrator->minStep = options->integratorMinTimestep;
}

int SourceSystemCVODE::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y, t);

    reactionRatesTimer->start();
    gas->setStateMass(Y, T);
    gas->getReactionRates(wDot);
    reactionRatesTimer->stop();

    updateThermo();

    double a = strainFunction.a(t);
    double dadt = strainFunction.dadt(t);

    // *** Calculate the time derivatives
    double scale;
    if (!quasi2d) {
        scale = 1.0;
        dUdt = - U*U + rhou/rho*(dadt + a*a) + splitConst[kMomentum];
        qDot = - (wDot * hk).sum() + getQdotIgniter(t);
        dTdt = qDot/(rho*cp) + splitConst[kEnergy];
    } else {
        scale = 1.0/vzInterp->get(x, t);
        dUdt = splitConst[kMomentum];
        dTdt = splitConst[kEnergy];
    }
    dYdt = scale * wDot * W / rho + splitConst.tail(nSpec);

    roll_ydot(ydot);
    return 0;
}

int SourceSystemCVODE::denseJacobian(const realtype t, const sdVector& y,
                                     const sdVector& ydot, sdMatrix& J)
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
    dvec dwdT(nSpec);
    dvec wDot2(nSpec);

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
    dvec hdwdY = dvec::Zero(nSpec);
    dvec YplusdY(nSpec);
    dvec drhodY(nSpec);

    double scale = (quasi2d) ? 1.0/vzInterp->get(x, t) : 1.0;

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
            J(kSpecies+k, kSpecies+i) = scale *
                (dwdY(k,i)*W[k]/rho - wDot[k]*W[k]*drhodY[i]/(rho*rho));
        }
        if (!quasi2d) {
            // dSpecies/dT
            J(kSpecies+k, kEnergy) = dwdT[k]*W[k]/rho -
                wDot[k]*W[k]*drhodT/(rho*rho);

            // dEnergy/dY
            J(kEnergy, kSpecies+k) = -hdwdY[k]/(rho*cp) -
                qDot*drhodY[k]/(rho*rho*cp);

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

int SourceSystemCVODE::fdJacobian(const realtype t, const sdVector& y,
                                  const sdVector& ydot, sdMatrix& J)
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

    jacobianTimer->stop();

    return 0;
}

void SourceSystemCVODE::setState
(double tInitial, double uu, double tt, const dvec& yy)
{
    integrator->t0 = tInitial;
    integrator->y[kMomentum] = uu;
    integrator->y[kEnergy] = tt;
    Eigen::Map<dvec>(&integrator->y[kSpecies], nSpec) = yy;
    integrator->initialize();
}

int SourceSystemCVODE::integrateToTime(double tf)
{
    return integrator->integrateToTime(integrator->t0 + tf);
}

int SourceSystemCVODE::integrateOneStep(double tf)
{
    return integrator->integrateOneStep(integrator->t0 + tf);
}

double SourceSystemCVODE::time() const
{
    return integrator->tInt - integrator->t0;
}

void SourceSystemCVODE::unroll_y(const sdVector& y, double t)
{
    if (!quasi2d) {
        T = y[kEnergy];
        U = y[kMomentum];
    } else {
        T = TInterp->get(x, t);
        U = 0;
    }
    Y = Eigen::Map<dvec>(&y[kSpecies], nSpec);
}

void SourceSystemCVODE::roll_y(sdVector& y) const
{
    y[kEnergy] = T;
    y[kMomentum] = U;
    Eigen::Map<dvec>(&y[kSpecies], nSpec) = Y;
}

void SourceSystemCVODE::roll_ydot(sdVector& ydot) const
{
    ydot[kEnergy] = dTdt;
    ydot[kMomentum] = dUdt;
    Eigen::Map<dvec>(&ydot[kSpecies], nSpec) = dYdt;
}

std::string SourceSystemCVODE::getStats()
{
    return (format("%i") % integrator->getNumSteps()).str();
}

void SourceSystemCVODE::writeJacobian(std::ostream& out)
{
    size_t N = integrator->y.length();
    double t = integrator->tInt;
    sdMatrix J(N,N);
    sdVector ydot(N);
    f(t, integrator->y, ydot);
    denseJacobian(t, integrator->y, ydot, J);

    out << "J = []" << std::endl;
    for (size_t i=0; i<N; i++) {
        out << "J.append([";
        for (size_t k=0; k<N; k++) {
            out << boost::format("%.5e, ") % J(i,k);
        }
        out << "])" << std::endl;
    }
}

// ----------------------------------------------------------------------------

SourceSystemQSS::SourceSystemQSS()
{
    integrator.setOde(this);
    dUdtQ = 0;
    dUdtD = 0;
    dTdtQ = 0;
    dTdtD = 0;
}

void SourceSystemQSS::initialize(size_t new_nSpec)
{
    SourceSystem::initialize(new_nSpec);
    integrator.initialize(new_nSpec + 2);

    dYdtQ.setConstant(nSpec, 0);
    dYdtD.setConstant(nSpec, 0);
    wDotD.resize(nSpec);
    wDotQ.resize(nSpec);

    integrator.enforce_ymin[kMomentum] = false;
}

void SourceSystemQSS::setOptions(ConfigOptions& opts)
{
    SourceSystem::setOptions(opts);
    integrator.epsmin = options->qss_epsmin;
    integrator.epsmax = options->qss_epsmax;
    integrator.dtmin = options->qss_dtmin;
    integrator.dtmax = options->qss_dtmax;
    integrator.itermax = options->qss_iterationCount;
    integrator.abstol = options->qss_abstol;
    integrator.stabilityCheck = options->qss_stabilityCheck;
    integrator.ymin.setConstant(nSpec + 2, options->qss_minval);
    integrator.ymin[kMomentum] = -1e4;
}

void SourceSystemQSS::setState
(double tStart, double uu, double tt, const dvec& yy)
{
    dvec yIn(nSpec + 2);
    yIn << uu, tt, yy;
    integrator.setState(yIn, tStart);
}

void SourceSystemQSS::odefun(double t, const dvec& y, dvec& q, dvec& d,
                             bool corrector)
{
    tCall = t;
    unroll_y(y, corrector);

    // *** Update auxiliary data ***
    reactionRatesTimer->start();
    gas->setStateMass(Y, T);
    gas->getCreationRates(wDotQ);
    gas->getDestructionRates(wDotD);
    reactionRatesTimer->stop();

    if (!corrector) {
        updateThermo();
    }

    qDot = - ((wDotQ - wDotD) * hk).sum() + getQdotIgniter(t);

    double a = strainFunction.a(t);
    double dadt = strainFunction.dadt(t);

    // *** Calculate the time derivatives
    double scale;
    if (!quasi2d) {
        scale = 1.0;
        dUdtQ = rhou/rho*(dadt + a*a) - U*U + splitConst[kMomentum];
        dUdtD = 0;
        dTdtQ = qDot/(rho*cp) + splitConst[kEnergy];
    } else {
        scale = 1.0/vzInterp->get(x, t);
        dUdtQ = 0;
        dUdtD = 0;
        dTdtQ = 0;
    }

    dYdtQ = scale * wDotQ * W / rho + splitConst.tail(nSpec);
    dYdtD = scale * wDotD * W / rho;

    assert(rhou > 0);
    assert(rho > 0);
    assert(dadt > -1e100 && dadt < 1e100);
    assert(a > -1e100 && a < 1e100);
    assert(U > -1e100 && U < 1e100);
    assert(splitConst[kMomentum] > -1e100 && splitConst[kMomentum] < 1e100);

    assert(dUdtQ > -1e100 && dUdtQ < 1e100);
    assert(dUdtD > -1e100 && dUdtD < 1e100);

    roll_ydot(q, d);
}

void SourceSystemQSS::unroll_y(const dvec& y, bool corrector)
{
    if (!quasi2d) {
        if (!corrector) {
            T = y[kEnergy];
        }
        U = y[kMomentum];
    } else {
        if (!corrector) {
            T = TInterp->get(x, tCall);
        }
        U = 0;
    }
    Y = y.tail(nSpec);
}

void SourceSystemQSS::roll_y(dvec& y) const
{
    y << U, T, Y;
}

void SourceSystemQSS::roll_ydot(dvec& q, dvec& d) const
{
    q << dUdtQ, dTdtQ, dYdtQ;
    d << dUdtD, dTdtD, dYdtD;
}

std::string SourceSystemQSS::getStats()
{
    return (format("%i/%i") % integrator.gcount % integrator.rcount).str();
}
