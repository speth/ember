#include "convectionSystem.h"

int ConvectionSystem::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y);

    // *** Update auxiliary data ***
    for (size_t j=0; j<nPoints; j++) {
        gas->setStateMass(Y[j],T[j]);
        rho[j] = gas->getDensity();
        Wmx[j] = gas->getMixtureMolecularWeight();
    }
    gas->getMolecularWeights(W);

    // *** Calculate V ***
    // rV[0] is a constant
    for (size_t j=0; j<nPoints-1; j++) {
        // Compute the upwinded convective derivative
        if (rV[j] < 0 || j == 0) {
            dTdx[j] = (T[j+1] - T[j]) / hh[j];
            dUdx[j] = (U[j+1] - U[j]) / hh[j];
            for (size_t k=0; k<nSpec; k++) {
                dYdx[j][k] = (Y[j+1][k] - Y[j][k]) / hh[j];
            }
        } else {
            dTdx[j] = (T[j] - T[j-1]) / hh[j-1];
            dUdx[j] = (U[j] - U[j-1]) / hh[j-1];
            for (size_t k=0; k<nSpec; k++) {
                dYdx[j][k] = (Y[j][k] - Y[j-1][k]) / hh[j-1];
            }
        }
        rV[j+1] = rV[j] - hh[j] * ((rV[j] * dTdx[j] + rphalf[j] * Tconst[j]) * rho[j] / T[j]);
        rV[j+1] -= hh[j] * rho[j] * U[j] * rphalf[j];
        for (size_t k=0; k<nSpec; k++) {
            rV[j+1] -= hh[j] * (rV[j] * dYdx[j][k] + rphalf[j] * Yconst[j][k])
                       * rho[j] * Wmx[j] / W[k];
        }
    }

    rV2V();

    // *** Calculate dY/dt, dU/dt, dT/dt

    // Left boundary conditions
    dUdt[0] = -Uconst[0]; // zero-gradient condition for U no matter what

    if (grid.leftBC == BoundaryCondition::FixedValue) {
        dTdt[0] = 0;
        for (size_t k=0; k<nSpec; k++) {
            dYdt[0][k] = 0;
        }
    } else if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero = std::max(rV[0], 0.0);

        dTdt[0] = -rVzero * (T[0] - Tleft) / (rho[0] * centerVol) - Tconst[0];
        for (size_t k=0; k<nSpec; k++) {
            dYdt[0][k] = -rVzero * (Y[0][k] - Yleft[k]) / (rho[0] * centerVol) - Yconst[0][k];
        }

    } else { // grid.leftBC == BoundaryCondition::ZeroGradient
        dTdt[0] = -Tconst[0];
        for (size_t k=0; k<nSpec; k++) {
            dYdt[0][k] = -Yconst[0][k];
        }
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        dUdt[j] = -(V[j] * dUdx[j]) / rho[j] - Uconst[j];
        dTdt[j] = -(V[j] * dTdx[j]) / rho[j] - Tconst[j];
        for (size_t k=0; k<nSpec; k++) {
            dYdt[j][k] = -(V[j] * dYdx[j][k]) / rho[j] - Yconst[j][k];
        }
    }

    // Right boundary values
    dUdt[jj] = -Uconst[jj]; // zero-gradient condition for U no matter what
    if (grid.rightBC == BoundaryCondition::FixedValue) {
        dTdt[jj] = 0;
        for (size_t k=0; k<nSpec; k++) {
            dYdt[jj][k] = 0;
        }
    } else { // grid.leftBC == BoundaryCondition::ZeroGradient
        dTdt[jj] = -Tconst[jj];
        for (size_t k=0; k<nSpec; k++) {
            dYdt[jj][k] = -Yconst[jj][k];
        }
    }

    roll_ydot(ydot);
    return 0;
}

void ConvectionSystem::unroll_y(const sdVector& y)
{
    for (size_t j=0; j<nPoints; j++) {
        T[j] = y[j*nVars+kEnergy];
        U[j] = y[j*nVars+kMomentum];
        for (size_t k=0; k<nSpec; k++) {
            Y[j][k] = y[j*nVars+kSpecies+k];
        }
    }
}

void ConvectionSystem::roll_y(sdVector& y) const
{
    for (size_t j=0; j<nPoints; j++) {
        y[j*nVars+kEnergy] = T[j];
        y[j*nVars+kMomentum] = U[j];
        for (size_t k=0; k<nSpec; k++) {
            y[j*nVars+kSpecies+k] = Y[j][k];
        }
    }
}

void ConvectionSystem::roll_ydot(sdVector& ydot) const
{
    for (size_t j=0; j<nPoints; j++) {
        ydot[j*nVars+kEnergy] = dTdt[j];
        ydot[j*nVars+kMomentum] = dUdt[j];
        for (size_t k=0; k<nSpec; k++) {
            ydot[j*nVars+kSpecies+k] = dYdt[j][k];
        }
    }

}

void ConvectionSystem::resize(const size_t new_nSpec, const size_t new_nPoints)
{
    if (new_nSpec == nSpec && new_nPoints == nPoints) {
        return;
    }

    nSpec = new_nSpec;
    nPoints = new_nPoints;
    U.resize(nPoints);
    dUdt.resize(nPoints);
    dUdx.resize(nPoints);
    Uconst.resize(nPoints);

    T.resize(nPoints);
    dTdt.resize(nPoints);
    dTdx.resize(nPoints);
    Tconst.resize(nPoints);

    rho.resize(nPoints);
    Wmx.resize(nPoints);
    rV.resize(nPoints);
    W.resize(nSpec);
    V.resize(nPoints);

    Y.resize(nPoints);
    dYdt.resize(nPoints);
    dYdx.resize(nPoints);
    Yconst.resize(nPoints);

    Yleft.resize(nSpec);
    Yright.resize(nSpec);
    for (size_t k=0; k<nPoints; k++) {
        Y[k].resize(nSpec);
        dYdt[k].resize(nSpec);
        dYdx[k].resize(nSpec);
        Yconst[k].resize(nSpec);
    }
}

void ConvectionSystem::V2rV(void)
{
    if (alpha == 0) {
        for (size_t j=0; j<nPoints; j++) {
            rV[j] = V[j];
        }
    } else {
        for (size_t j=0; j<nPoints; j++) {
            rV[j] = x[j]*V[j];
        }
    }

    // In the case where the centerline is part of the domain,
    // V[0] actually stores rV[0] (since V = Inf there)
    if (x[0] == 0) {
        rV[0] = V[0];
    }
}

void ConvectionSystem::rV2V(void)
{
    if (alpha == 0) {
        for (size_t j=0; j<nPoints; j++) {
            V[j] = rV[j];
        }
    } else {
        for (size_t j=0; j<nPoints; j++) {
            V[j] = rV[j]/x[j];
        }
    }

    // In the case where the centerline is part of the domain,
    // V[0] actually stores rV[0] (since V = Inf there)
    if (x[0] == 0) {
        V[0] = rV[0];
    }
}
