#include "convectionSystem.h"

ConvectionSystem::ConvectionSystem()
    : gas(NULL)
    , nSpec(0)
    , nPoints(0)
    , nVars(0)
{
}

int ConvectionSystem::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y);

    // *** Update auxiliary data ***
    for (size_t j=0; j<nPoints; j++) {
        gas->setStateMass(&Y(0,j), T[j]);
        rho[j] = gas->getDensity();
        Wmx[j] = gas->getMixtureMolecularWeight();
    }
    gas->getMolecularWeights(W);

    // *** Calculate V ***
    rV[0] = rVzero;
    for (size_t j=0; j<nPoints-1; j++) {
        // Compute the upwinded convective derivative
        if (rV[j] < 0 || j == 0) {
            dTdx[j] = (T[j+1] - T[j]) / hh[j];
            dUdx[j] = (U[j+1] - U[j]) / hh[j];
            for (size_t k=0; k<nSpec; k++) {
                dYdx(k,j) = (Y(k,j+1) - Y(k,j)) / hh[j];
            }
        } else {
            dTdx[j] = (T[j] - T[j-1]) / hh[j-1];
            dUdx[j] = (U[j] - U[j-1]) / hh[j-1];
            for (size_t k=0; k<nSpec; k++) {
                dYdx(k,j) = (Y(k,j) - Y(k,j-1)) / hh[j-1];
            }
        }
        rV[j+1] = rV[j] - hh[j] * ((rV[j] * dTdx[j] + rphalf[j] * Tconst[j]) * rho[j] / T[j]);
        rV[j+1] -= hh[j] * rho[j] * U[j] * rphalf[j];
        for (size_t k=0; k<nSpec; k++) {
            rV[j+1] -= hh[j] * (rV[j] * dYdx(k,j) + rphalf[j] * Yconst(k,j))
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
            dYdt(k,0) = 0;
        }
    } else if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero_mod = std::max(rV[0], 0.0);

        dTdt[0] = -rVzero_mod * (T[0] - Tleft) / (rho[0] * centerVol) - Tconst[0];
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,0) = -rVzero_mod * (Y(k,0) - Yleft[k]) / (rho[0] * centerVol) - Yconst(k,0);
        }

    } else { // grid.leftBC == BoundaryCondition::ZeroGradient
        dTdt[0] = -Tconst[0];
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,0) = -Yconst(k,0);
        }
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        dUdt[j] = -(V[j] * dUdx[j]) / rho[j] - Uconst[j];
        dTdt[j] = -(V[j] * dTdx[j]) / rho[j] - Tconst[j];
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,j) = -(V[j] * dYdx(k,j)) / rho[j] - Yconst(k,j);
        }
    }

    // Right boundary values
    dUdt[jj] = -Uconst[jj]; // zero-gradient condition for U no matter what
    if (grid.rightBC == BoundaryCondition::FixedValue) {
        dTdt[jj] = 0;
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,jj) = 0;
        }
    } else { // grid.leftBC == BoundaryCondition::ZeroGradient
        dTdt[jj] = -Tconst[jj];
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,jj) = -Yconst(k,jj);
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
            Y(k,j) = y[j*nVars+kSpecies+k];
        }
    }
}

void ConvectionSystem::roll_y(sdVector& y) const
{
    for (size_t j=0; j<nPoints; j++) {
        y[j*nVars+kEnergy] = T[j];
        y[j*nVars+kMomentum] = U[j];
        for (size_t k=0; k<nSpec; k++) {
            y[j*nVars+kSpecies+k] = Y(k,j);
        }
    }
}

void ConvectionSystem::roll_ydot(sdVector& ydot) const
{
    for (size_t j=0; j<nPoints; j++) {
        ydot[j*nVars+kEnergy] = dTdt[j];
        ydot[j*nVars+kMomentum] = dUdt[j];
        for (size_t k=0; k<nSpec; k++) {
            ydot[j*nVars+kSpecies+k] = dYdt(k,j);
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
    nVars = nSpec + 2;
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

    Y.resize(nSpec, nPoints);
    dYdt.resize(nSpec, nPoints);
    dYdx.resize(nSpec, nPoints);
    Yconst.resize(nSpec, nPoints);
}

void ConvectionSystem::V2rV(void)
{
    if (alpha == 0) {
        for (size_t j=1; j<nPoints; j++) {
            rV[j] = V[j];
        }
    } else {
        for (size_t j=1; j<nPoints; j++) {
            rV[j] = x[j]*V[j];
        }
    }
}

void ConvectionSystem::rV2V(void)
{
    if (alpha == 0) {
        for (size_t j=1; j<nPoints; j++) {
            V[j] = rV[j];
        }
    } else {
        for (size_t j=1; j<nPoints; j++) {
            V[j] = rV[j]/x[j];
        }
    }
}
