#include "sourceSystem.h"
#include "readConfig.h"

SourceSystem::SourceSystem()
    : gas(NULL)
{
}

void SourceSystem::resize(size_t new_nSpec)
{
    nSpec = new_nSpec;
    Y.resize(nSpec);
    dYdt.resize(nSpec);
    cpSpec.resize(nSpec);
    wDot.resize(nSpec);
    hk.resize(nSpec);
    C.resize(nSpec+2);
}

int SourceSystem::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y);

    // *** Update auxiliary data ***
    gas->setStateMass(Y, T);
    gas->getReactionRates(wDot);
    gas->getEnthalpies(hk);
    rho = gas->getDensity();
    cp = gas->getSpecificHeatCapacity();

    qDot = 0.0;
    for (size_t k=0; k<nSpec; k++) {
        qDot -= wDot[k]*hk[k];
    }

    double a = strainFunction.a(t);
    double dadt = strainFunction.dadt(t);

    // *** Calculate the time derivatives
    dUdt = - U*U + rhou/rho*(dadt + a*a) + C[kMomentum];
    dTdt = qDot/(rho*cp) + C[kEnergy];
    for (size_t k=0; k<nSpec; k++) {
        dYdt[k] = wDot[k]*W[k]/rho + C[kSpecies+k];
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
    gas->getSpecificHeatCapacities(cpSpec);
    Wmx = gas->getMixtureMolecularWeight();

    // The constant "800" here has been empirically determined to give
    // good performance for typical test cases. This value can have
    // a substantial impact on the convergence rate of the solver.
    double eps = sqrt(DBL_EPSILON)*800;

    // *** Derivatives with respect to temperature
    double TplusdT = T*(1+eps);

    double dwhdT = 0;
    dvector dwdT(nSpec);
    dvector wDot2(nSpec);

    gas->setStateMass(Y, TplusdT);
    gas->getReactionRates(wDot2);

    for (size_t k=0; k<nSpec; k++) {
        dwdT[k] = (wDot2[k]-wDot[k])/(TplusdT-T);
        dwhdT += hk[k]*dwdT[k] + cpSpec[k]*wDot[k];
    }

    double drhodT = -rho/T;

    // *** Derivatives with respect to species concentration
    Array2D dwdY(nSpec, nSpec);
    dvector hdwdY(nSpec, 0);
    dvector YplusdY(nSpec);
    dvector drhodY(nSpec);

    for (size_t k=0; k<nSpec; k++) {
        YplusdY[k] = (abs(Y[k]) > eps/2) ? Y[k]*(1+eps) : eps;
        gas->setStateMass(YplusdY, T);
        gas->getReactionRates(wDot2);

        for (size_t i=0; i<nSpec; i++) {
            dwdY(i,k) = (wDot2[i]-wDot[i])/(YplusdY[k]-Y[k]);
            hdwdY[k] += hk[i]*dwdY(i,k);
        }

        drhodY[k] = rho*(W[k]-Wmx)/(W[k]*(1-Y[k]*(1-eps)));
    }

    for (size_t k=0; k<nSpec; k++) {
        for (size_t i=0; i<nSpec; i++) {
            // dSpecies/dY
            J(kSpecies+k, kSpecies+i) = dwdY(k,i)*W[k]/rho - wDot[k]*W[k]*drhodY[i]/(rho*rho);
        }
        // dSpecies/dT
        J(kSpecies+k, kEnergy) = dwdT[k]*W[k]/rho - wDot[k]*W[k]*drhodT/(rho*rho);

        // dEnergy/dY
        J(kEnergy, kSpecies+k) = -hdwdY[k]/(rho*cp) - qDot*drhodY[k]/(rho*rho*cp);

        // dMomentum/dY
        J(kMomentum, kSpecies+k) = -A*drhodY[k]/(rho*rho);
    }

    // dEnergy/dT
    J(kEnergy, kEnergy) = -dwhdT/(rho*cp) - qDot*drhodT/(rho*rho*cp);

    // dMomentum/dU
    J(kMomentum, kMomentum) = -2*U;

    // dMomentum/dT
    J(kMomentum, kEnergy) = -A*drhodT/(rho*rho);

    return 0;
}

int SourceSystem::fdJacobian(const realtype t, const sdVector& y, const sdVector& ydot, sdMatrix& J)
{
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
    return 0;
}

void SourceSystem::unroll_y(const sdVector& y)
{
    T = y[kEnergy];
    U = y[kMomentum];
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
