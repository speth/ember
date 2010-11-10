#include "convectionSystem.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

ConvectionSystem::ConvectionSystem()
    : gas(NULL)
    , nSpec(0)
    , nVars(0)
{
}

int ConvectionSystem::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y);
    // *** Update auxiliary data ***
    thermoTimer->start();
    for (size_t j=0; j<nPoints; j++) {
        gas->setStateMass(&Y(0,j), T[j]);
        rho[j] = gas->getDensity();
        Wmx[j] = gas->getMixtureMolecularWeight();
    }
    gas->getMolecularWeights(W);
    thermoTimer->stop();

    // Update split terms
    for (size_t j=0; j<nPoints; j++) {
        Uconst[j] = splitConstU[j] + splitLinearU[j] * U[j];
        Tconst[j] = splitConstT[j] + splitLinearT[j] * T[j];
        for (size_t k=0; k<nSpec; k++) {
            Yconst(k,j) = splitConstY(k,j) + splitLinearY(k,j) * Y(k,j);
        }
    }

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

        rV[j+1] = rV[j] - hh[j] * (rV[j] * dTdx[j] - rphalf[j] * rho[j] * Tconst[j]) / T[j];
        rV[j+1] -= hh[j] * rho[j] * U[j] * rphalf[j];
        for (size_t k=0; k<nSpec; k++) {
            rV[j+1] -= hh[j] * (rV[j] * dYdx(k,j) - rphalf[j] * rho[j] * Yconst(k,j))
                       * Wmx[j] / W[k];
        }
    }

    rV2V();

    // *** Calculate dY/dt, dU/dt, dT/dt

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case
    dUdt[0] = Uconst[0]; // zero-gradient condition for U no matter what

    if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero_mod = std::max(rV[0], 0.0);

        dTdt[0] = -rVzero_mod * (T[0] - Tleft) / (rho[0] * centerVol) + Tconst[0];
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,0) = -rVzero_mod * (Y(k,0) - Yleft[k]) / (rho[0] * centerVol) + Yconst(k,0);
        }

    } else { // FixedValue or ZeroGradient
        dTdt[0] = Tconst[0];
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,0) = Yconst(k,0);
        }
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        dUdt[j] = -V[j] * dUdx[j] / rho[j] + Uconst[j];
        dTdt[j] = -V[j] * dTdx[j] / rho[j] + Tconst[j];
        for (size_t k=0; k<nSpec; k++) {
            dYdt(k,j) = -V[j] * dYdx(k,j) / rho[j] + Yconst(k,j);
        }
    }

    // Right boundary values
    // Convection term has nothing to contribute in any case,
    // So only the value from the other terms remains
    dUdt[jj] = Uconst[jj];
    dTdt[jj] = Tconst[jj];
    for (size_t k=0; k<nSpec; k++) {
        dYdt(k,jj) = Yconst(k,jj);
    }

    roll_ydot(ydot);
    return 0;
}

void ConvectionSystem::get_diagonal(const realtype t, dvector& linearU,
                                    dvector& linearT, Array2D& linearY)
{
    // Assume that f has just been called and that all auxiliary
    // arrays are in a consistent state.

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case

    // dUdot/dU
    linearU[0] = 0;

    if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero_mod = std::max(rV[0], 0.0);

        // dTdot/dT
        linearT[0] = -rVzero_mod / (rho[0] * centerVol);

        for (size_t k=0; k<nSpec; k++) {
            // dYdot/dY
            linearY(k,0) = -rVzero_mod / (rho[0] * centerVol);
        }

    } else { // FixedValue or ZeroGradient
        linearT[0] = 0;
        for (size_t k=0; k<nSpec; k++) {
            linearY(k,0) = 0;
        }
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        // depends on upwinding to calculated dT/dx etc.
        double value = (rV[j] < 0 || j == 0) ? V[j] / (rho[j] * hh[j])
                                             : -V[j] / (rho[j] * hh[j-1]);
        linearU[j] = value;
        linearT[j] = value;
        for (size_t k=0; k<nSpec; k++) {
            linearY(k,j) = value;
        }
    }
}

int ConvectionSystem::bandedJacobian(const realtype t, const sdVector& y,
                                     const sdVector& ydot, sdBandMatrix& J)
{
    // Assume that f has just been called and that all auxiliary
    // arrays are in a consistent state.

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case

    // dUdot/dU
    J(0*nVars+kMomentum, 0*nVars+kMomentum) = splitLinearU[0];

    if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero_mod = std::max(rV[0], 0.0);

        // dTdot/dT
        J(0*nVars+kEnergy, 0*nVars+kEnergy) =
                -rVzero_mod / (rho[0] * centerVol) + splitLinearT[0];

        for (size_t k=0; k<nSpec; k++) {
            // dYdot/dY
            J(0*nVars+kSpecies+k, 0*nVars+kSpecies+k) =
                    -rVzero_mod / (rho[0] * centerVol) + splitLinearY(k,0);
        }

    } else { // FixedValue or ZeroGradient
        J(0*nVars+kEnergy, 0*nVars+kEnergy) = splitLinearT[0];
        for (size_t k=0; k<nSpec; k++) {
            J(0*nVars+kSpecies+k, 0*nVars+kSpecies+k) = splitLinearY(k,0);
        }
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        // depends on upwinding to calculated dT/dx etc.
        double value = (rV[j] < 0 || j == 0) ? V[j] / (rho[j] * hh[j])
                                             : -V[j] / (rho[j] * hh[j-1]);
        J(j*nVars+kMomentum, j*nVars+kMomentum) = value + splitLinearU[j];
        J(j*nVars+kEnergy, j*nVars+kEnergy) = value + splitLinearT[j];
        for (size_t k=0; k<nSpec; k++) {
            J(j*nVars+kSpecies+k, j*nVars+kSpecies+k) = value + splitLinearY(k,j);
        }
    }

    // Right boundary conditions
    J(jj*nVars+kMomentum, jj*nVars+kMomentum) = splitLinearU[jj];
    J(jj*nVars+kEnergy, jj*nVars+kEnergy) = splitLinearT[jj];
    for (size_t k=0; k<nSpec; k++) {
        J(jj*nVars+kSpecies+k, jj*nVars+kSpecies+k) = splitLinearY(k,jj);
    }

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
    nSpec = new_nSpec;
    grid.setSize(new_nPoints);
    nVars = nSpec + 2;
    U.resize(nPoints);
    dUdt.resize(nPoints);
    dUdx.resize(nPoints);
    Uconst.resize(nPoints);
    splitConstU.resize(nPoints);
    splitLinearU.resize(nPoints);

    T.resize(nPoints);
    dTdt.resize(nPoints);
    dTdx.resize(nPoints);
    Tconst.resize(nPoints);
    splitConstT.resize(nPoints);
    splitLinearT.resize(nPoints);

    rho.resize(nPoints);
    Wmx.resize(nPoints);
    rV.resize(nPoints);
    W.resize(nSpec);
    V.resize(nPoints);

    Y.resize(nSpec, nPoints);
    dYdt.resize(nSpec, nPoints);
    dYdx.resize(nSpec, nPoints);
    Yconst.resize(nSpec, nPoints);
    splitConstY.resize(nSpec, nPoints);
    splitLinearY.resize(nSpec, nPoints);
}

void ConvectionSystem::initialize()
{
    splitConstU.assign(nPoints, 0);
    splitLinearU.assign(nPoints, 0);

    splitConstT.assign(nPoints, 0);
    splitLinearT.assign(nPoints, 0);
    splitConstY.data().assign(nPoints*nSpec, 0);
    splitLinearY.data().assign(nPoints*nSpec, 0);
}

void ConvectionSystem::V2rV(void)
{
    V[0] = rV[0];
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
    V[0] = rV[0];
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



ConvectionSystemUTW::ConvectionSystemUTW()
    : gas(NULL)
    , nVars(3)
{
}

int ConvectionSystemUTW::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    unroll_y(y);
    // *** Update auxiliary data ***
    for (size_t j=0; j<nPoints; j++) {
        rho[j] = gas->pressure*Wmx[j]/(Cantera::GasConstant*T[j]);
    }

    // Update split terms
    for (size_t j=0; j<nPoints; j++) {
        Uconst[j] = splitConstU[j] + splitLinearU[j] * U[j];
        Tconst[j] = splitConstT[j] + splitLinearT[j] * T[j];
    }

    // *** Calculate V ***
    rV[0] = rVzero;
    for (size_t j=0; j<nPoints-1; j++) {
        // Compute the upwinded convective derivative
        if (rV[j] < 0 || j == 0) {
            dTdx[j] = (T[j+1] - T[j]) / hh[j];
            dUdx[j] = (U[j+1] - U[j]) / hh[j];
            dWdx[j] = (Wmx[j+1] - Wmx[j]) / hh[j];
        } else {
            dTdx[j] = (T[j] - T[j-1]) / hh[j-1];
            dUdx[j] = (U[j] - U[j-1]) / hh[j-1];
            dWdx[j] = (Wmx[j] - Wmx[j-1]) / hh[j-1];
        }

        rV[j+1] = rV[j] - hh[j] * (rV[j] * dTdx[j] - rphalf[j] * rho[j] * Tconst[j]) / T[j];
        rV[j+1] += hh[j] * (rV[j] * dWdx[j] - rphalf[j] * rho[j] * splitConstW[j]) / Wmx[j];
        rV[j+1] -= hh[j] * rho[j] * U[j] * rphalf[j];
    }

    rV2V();

    // *** Calculate dW/dt, dU/dt, dT/dt

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case
    dUdt[0] = Uconst[0]; // zero-gradient condition for U no matter what

    if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero_mod = std::max(rV[0], 0.0);

        dTdt[0] = -rVzero_mod * (T[0] - Tleft) / (rho[0] * centerVol) + Tconst[0];
        dWdt[0] = -rVzero_mod * (Wmx[0] - Wleft) / (rho[0] * centerVol) + splitConstW[0];

    } else { // FixedValue or ZeroGradient
        dTdt[0] = Tconst[0];
        dWdt[0] = splitConstW[0];
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        dUdt[j] = -V[j] * dUdx[j] / rho[j] + Uconst[j];
        dTdt[j] = -V[j] * dTdx[j] / rho[j] + Tconst[j];
        dWdt[j] = -V[j] * dWdx[j] / rho[j] + splitConstW[0];
    }

    // Right boundary values
    // Convection term has nothing to contribute in any case,
    // So only the value from the other terms remains
    dUdt[jj] = Uconst[jj];
    dTdt[jj] = Tconst[jj];
    dWdt[jj] = splitConstW[jj];

    roll_ydot(ydot);
    return 0;
}

void ConvectionSystemUTW::get_diagonal
(const realtype t, dvector& linearU, dvector& linearT)
{
    // Assume that f has just been called and that all auxiliary
    // arrays are in a consistent state.

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case

    // dUdot/dU
    linearU[0] = 0;

    if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero_mod = std::max(rV[0], 0.0);

        // dTdot/dT
        linearT[0] = -rVzero_mod / (rho[0] * centerVol);

    } else { // FixedValue or ZeroGradient
        linearT[0] = 0;
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        // depends on upwinding to calculated dT/dx etc.
        double value = (rV[j] < 0 || j == 0) ? V[j] / (rho[j] * hh[j])
                                             : -V[j] / (rho[j] * hh[j-1]);
        linearU[j] = value;
        linearT[j] = value;
    }
}

int ConvectionSystemUTW::bandedJacobian(const realtype t, const sdVector& y,
                                     const sdVector& ydot, sdBandMatrix& J)
{
    // Assume that f has just been called and that all auxiliary
    // arrays are in a consistent state.

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case

    // dUdot/dU
    J(0*nVars+kMomentum, 0*nVars+kMomentum) = splitLinearU[0];

    if (grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rVzero_mod = std::max(rV[0], 0.0);

        // dTdot/dT
        J(0*nVars+kEnergy, 0*nVars+kEnergy) =
                -rVzero_mod / (rho[0] * centerVol) + splitLinearT[0];

        // dWdot/dW
        J(0*nVars+kWmx, 0*nVars+kWmx) = -rVzero_mod / (rho[0] * centerVol);

    } else { // FixedValue or ZeroGradient
        J(0*nVars+kEnergy, 0*nVars+kEnergy) = splitLinearT[0];
        // J(0*nVars+kWmx, 0*nVars+kWmx) = 0;
    }

    // Intermediate points
    for (size_t j=1; j<jj; j++) {
        // depends on upwinding to calculated dT/dx etc.
        double value = (rV[j] < 0 || j == 0) ? V[j] / (rho[j] * hh[j])
                                             : -V[j] / (rho[j] * hh[j-1]);
        J(j*nVars+kMomentum, j*nVars+kMomentum) = value + splitLinearU[j];
        J(j*nVars+kEnergy, j*nVars+kEnergy) = value + splitLinearT[j];
        J(j*nVars+kWmx, j*nVars+kWmx) = value;
    }

    // Right boundary conditions
    J(jj*nVars+kMomentum, jj*nVars+kMomentum) = splitLinearU[jj];
    J(jj*nVars+kEnergy, jj*nVars+kEnergy) = splitLinearT[jj];
    // J(jj*nVars+kWmx, jj*nVars+kWmx) = 0;

    return 0;
}

void ConvectionSystemUTW::unroll_y(const sdVector& y)
{
    for (size_t j=0; j<nPoints; j++) {
        T[j] = y[j*nVars+kEnergy];
        U[j] = y[j*nVars+kMomentum];
        Wmx[j] = y[j*nVars+kWmx];
    }
}

void ConvectionSystemUTW::roll_y(sdVector& y) const
{
    for (size_t j=0; j<nPoints; j++) {
        y[j*nVars+kEnergy] = T[j];
        y[j*nVars+kMomentum] = U[j];
        y[j*nVars+kWmx] = Wmx[j];
    }
}

void ConvectionSystemUTW::roll_ydot(sdVector& ydot) const
{
    for (size_t j=0; j<nPoints; j++) {
        ydot[j*nVars+kEnergy] = dTdt[j];
        ydot[j*nVars+kMomentum] = dUdt[j];
        ydot[j*nVars+kWmx] = dWdt[j];
    }
}

void ConvectionSystemUTW::resize(const size_t new_nPoints)
{
    grid.setSize(new_nPoints);
    rho.resize(nPoints);
    rV.resize(nPoints);
    V.resize(nPoints);

    U.resize(nPoints);
    dUdt.resize(nPoints);
    dUdx.resize(nPoints);
    Uconst.resize(nPoints);
    splitConstU.resize(nPoints);
    splitLinearU.resize(nPoints);

    T.resize(nPoints);
    dTdt.resize(nPoints);
    dTdx.resize(nPoints);
    Tconst.resize(nPoints);
    splitConstT.resize(nPoints);
    splitLinearT.resize(nPoints);

    Wmx.resize(nPoints);
    dWdt.resize(nPoints);
    dWdx.resize(nPoints);
    splitConstW.resize(nPoints);
}

void ConvectionSystemUTW::initialize()
{
    splitConstU.assign(nPoints, 0);
    splitLinearU.assign(nPoints, 0);
    splitConstT.assign(nPoints, 0);
    splitLinearT.assign(nPoints, 0);
    splitConstW.assign(nPoints, 0);
}

void ConvectionSystemUTW::V2rV(void)
{
    V[0] = rV[0];
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

void ConvectionSystemUTW::rV2V(void)
{
    V[0] = rV[0];
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


ConvectionSystemY::ConvectionSystemY()
{
}

int ConvectionSystemY::f(const realtype t, const sdVector& y, sdVector& ydot)
{
    // *** Calculate v (= V/rho) ***
    update_v(t);

    // *** Calculate dY/dt

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case
    if (startIndex == 0 && grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        // Note: v[0] actually contains r*v[0] in this case
        double rvzero_mod = std::max(v[0], 0.0);
        ydot[0] = -rvzero_mod * (y[0] - Yleft) / centerVol
                + splitConstY[0] + splitLinearY[0] * y[0];
    } else { // FixedValue, ZeroGradient, or truncated domain
        ydot[0] = splitConstY[0] + splitLinearY[0] * y[0];
    }

    // Intermediate points
    double dYdx;
    size_t i = 1;
    for (size_t j=startIndex+1; j<stopIndex; j++) {
        if (v[i] < 0) {
            dYdx = (y[i+1] - y[i]) / hh[j];
        } else {
            dYdx = (y[i] - y[i-1]) / hh[j-1];
        }
        ydot[i] = -v[i] * dYdx + splitConstY[i] + splitLinearY[i] * y[i];
        i++;
    }

    // Right boundary values
    // Convection term has nothing to contribute in any case,
    // So only the value from the other terms remains
    ydot[i] = splitConstY[i] + splitLinearY[i] * y[i];

    return 0;
}

void ConvectionSystemY::get_diagonal(const realtype t, dvector& linearY)
{
    // Assume that f has just been called and that all auxiliary
    // arrays are in a consistent state.

    // Left boundary conditions.
    // Convection term only contributes in the ControlVolume case

    if (startIndex == 0 && grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        // Note: v[0] actually contains r*v[0] in this case
        double rvzero_mod = std::max(v[0], 0.0);

        linearY[0] = -rvzero_mod / centerVol;

    } else { // FixedValue, ZeroGradient, or truncated domain
        linearY[0] = 0;
    }

    // Intermediate points
    size_t i = 1;
    for (size_t j=startIndex+1; j<stopIndex; j++) {
        // depends on upwinding to calculated dT/dx etc.
        linearY[i] = (v[i] < 0) ? v[i] / hh[j]
                                : -v[i] / hh[j-1];
        i++;
    }
    linearY[i] = 0;
}

int ConvectionSystemY::bandedJacobian(const realtype t, const sdVector& y,
                                     const sdVector& ydot, sdBandMatrix& J)
{
    // Assume that f has just been called and that all auxiliary
    // arrays are in a consistent state.

    // Left boundary condition.
    // Convection term only contributes in the ControlVolume case

    if (startIndex == 0 && grid.leftBC == BoundaryCondition::ControlVolume) {
        double centerVol = pow(x[1],alpha+1) / (alpha+1);
        double rvzero_mod = std::max(v[0], 0.0);

        J(0,0) = -rvzero_mod / centerVol + splitLinearY[0];

    } else { // FixedValue, ZeroGradient or truncated domain
        J(0,0) = splitLinearY[0];
    }

    // Intermediate points
    size_t i = 1;
    for (size_t j=startIndex+1; j<stopIndex; j++) {
        // depends on upwinding to calculated dT/dx etc.
        double value = (v[i] < 0) ? v[i] / hh[j]
                                  : -v[i] / hh[j-1];
        J(i,i) = value + splitLinearY[i];
        i++;
    }

    // Right boundary condition
    J(i,i) = splitLinearY[i];

    return 0;
}

void ConvectionSystemY::resize(const size_t new_nPoints)
{
    grid.setSize(new_nPoints);

    v.resize(nPoints);
    splitConstY.resize(nPoints);
    splitLinearY.resize(nPoints);
}

void ConvectionSystemY::initialize()
{
    splitConstY.assign(nPoints, 0);
    splitLinearY.assign(nPoints, 0);
}


void ConvectionSystemY::update_v(const double t)
{
    assert(vInterp->size() > 0);

    if (vInterp->size() == 1) {
        // vInterp only has one data point
        const dvector& vLeft = vInterp->begin()->second;
        size_t i = 0;
        for (size_t j=startIndex; j<=stopIndex; j++) {
            v[i] = vLeft[j];
            i++;
        }
        return;
    }

    // Find the value of v by interpolating between values contained in vInterp
    vecInterpolator::iterator iLeft = vInterp->lower_bound(t);
    if (iLeft == vInterp->end()) {
        // In this case, we're actually extrapolating past the rightmost point
        iLeft--;
    }
    vecInterpolator::iterator iRight = iLeft;
    iRight++;
    if (iRight == vInterp->end()) {
        iLeft--;
        iRight--;
    }
    const dvector& vLeft = iLeft->second;
    const dvector& vRight = iRight->second;

    // Linear interpolation
    double s = (t-iLeft->first)/(iRight->first - iLeft->first);

    size_t i = 0;
    for (size_t j=startIndex; j<=stopIndex; j++) {
        v[i] = vLeft[j]*(1-s) + vRight[j]*s;
    }
}

ConvectionSystemSplit::ConvectionSystemSplit()
    : nSpec(0)
    , nVars(3)
    , nPointsUTW(0)
{
}

void ConvectionSystemSplit::get_diagonal
(const realtype t, dvector& dU, dvector& dT, Array2D& dY)
{
    // Called after evaluate()
    utwSystem.get_diagonal(t, dU, dT);

    for (size_t k=0; k<nSpec; k++) {
        dvector dYk(nPointsSpec[k], 0);
        speciesSystems[k].get_diagonal(t, dYk);
        size_t i = 0;
        for (size_t j=(*startIndices)[k]; j<(*stopIndices)[k]; j++) {
            dY(k,j) = dYk[i];
            i++;
        }
    }
}

void ConvectionSystemSplit::setGrid(const oneDimGrid& grid)
{
    GridBased::setGrid(grid);
    utwSystem.setGrid(grid);
    foreach (ConvectionSystemY& system, speciesSystems) {
        system.setGrid(grid);
    }
}

void ConvectionSystemSplit::setTolerances(const configOptions& options)
{
    reltol = options.idaRelTol;
    abstolU = options.idaMomentumAbsTol;
    abstolT = options.idaEnergyAbsTol;
    abstolW = options.idaSpeciesAbsTol * 20;
    abstolY = options.idaSpeciesAbsTol;
}

void ConvectionSystemSplit::setGas(CanteraGas& gas_)
{
    gas = &gas_;
    utwSystem.gas = &gas_;
}

void ConvectionSystemSplit::setThermoTimer(perfTimer& timer)
{
    thermoTimer = &timer;
}

void ConvectionSystemSplit::resize
(const size_t nPointsUTWNew, const vector<size_t>& nPointsSpecNew, const size_t nSpecNew)
{
    // Create or destroy the necessary speciesSystems if nSpec has changed
    if (nSpec != nSpecNew) {
        speciesSystems.resize(nSpecNew);
        nSpec = nSpecNew;
    }

    if (speciesSolvers.size() != nSpec) {
        // Create speciesSolvers from scratch if the number of species has changed
        speciesSolvers.clear();
        nPointsSpec = nPointsSpecNew;
        for (size_t k=0; k<nSpec; k++) {
            speciesSolvers.push_back(new sundialsCVODE(nPointsSpec[k]));
            configureSolver(speciesSolvers[k], k);
        }
    } else {
        // Replace the solvers where the number of points has changed
        for (size_t k=0; k<nSpec; k++) {
            if (nPointsSpec[k] != nPointsSpecNew[k]) {
                speciesSolvers.replace(k, new sundialsCVODE(nPointsSpecNew[k]));
                nPointsSpec[k] = nPointsSpecNew[k];
                configureSolver(speciesSolvers[k], k);
            }
        }
        nPointsSpec = nPointsSpecNew;
    }

    // Recreate the UTW solver if necessary
    if (nPointsUTW != nPointsUTWNew) {
        nPointsUTW = nPointsUTWNew;
        utwSolver.reset(new sundialsCVODE(3*nPointsUTW));
        utwSolver->setODE(&utwSystem);
        utwSolver->setBandwidth(0,0);
        utwSolver->reltol = reltol;
        utwSolver->linearMultistepMethod = CV_BDF;
        utwSolver->nonlinearSolverMethod = CV_NEWTON;
        for (size_t j=0; j<nPointsUTW; j++) {
            utwSolver->abstol[3*j+kMomentum] = abstolU;
            utwSolver->abstol[3*j+kEnergy] = abstolT;
            utwSolver->abstol[3*j+kWmx] = abstolW;
        }
        utwSystem.resize(nPointsUTW);
    }

    nPoints = nPointsUTWNew;
}

void ConvectionSystemSplit::setSpeciesDomains
(vector<size_t>& startIndices_, vector<size_t>& stopIndices_)
{
    startIndices = &startIndices_;
    stopIndices = &stopIndices_;
}

void ConvectionSystemSplit::setState(const dvector& U_, const dvector& T_, Array2D& Y_)
{
    U = U_;
    T = T_;
    Y = Y_;

    for (size_t j=0; j<nPointsUTW; j++) {
        utwSolver->y[3*j+kMomentum] = U[j];
        utwSolver->y[3*j+kEnergy] = T[j];
        gas->setStateMass(&Y(0,j), T[j]);
        utwSolver->y[3*j+kWmx] = gas->getMixtureMolecularWeight();
    }

    for (size_t k=0; k<nSpec; k++) {
        size_t i = 0;
        for (size_t j=(*startIndices)[k]; j<=(*stopIndices)[k]; j++) {
            speciesSolvers[k].y[i] = Y(k,j);
            i++;
        }
    }
}

void ConvectionSystemSplit::setLeftBC(const double Tleft, const dvector& Yleft_)
{
    utwSystem.Tleft = Tleft;
    Yleft = Yleft_;
}

void ConvectionSystemSplit::set_rVzero(const double rVzero)
{
    utwSystem.rVzero = rVzero;
}

void ConvectionSystemSplit::initialize(const double t0)
{
    // Initialize systems
    utwSystem.initialize();

    foreach (ConvectionSystemY& system, speciesSystems) {
        system.initialize();
    }

    // Initialize solvers
    utwSolver->t0 = t0;
    utwSolver->maxNumSteps = 1000000;
    utwSolver->minStep = 1e-16;
    utwSolver->initialize();

    foreach (sundialsCVODE& solver, speciesSolvers) {
        solver.t0 = t0;
        solver.maxNumSteps = 1000000;
        solver.minStep = 1e-16;
        solver.initialize();
    }
}

void ConvectionSystemSplit::evaluate()
{
    sdVector ydotUTW(nVars*nPoints);
    utwSystem.f(utwSolver->tInt, utwSolver->y, ydotUTW);

    boost::shared_ptr<vecInterpolator> vInterp(new vecInterpolator());
    vInterp->insert(std::make_pair(utwSolver->tInt, utwSystem.V));

    V = utwSystem.V;
    dUdt = utwSystem.dUdt;
    dTdt = utwSystem.dTdt;

    dYdt.data().clear();
    dYdt.resize(nSpec, nPoints, 0);
    sdVector ydotk(nPoints);
    for (size_t k=0; k<nSpec; k++) {
        speciesSystems[k].vInterp = vInterp;
        speciesSystems[k].f(speciesSolvers[k].tInt, speciesSolvers[k].y, ydotk);
        size_t i = 0;
        for (size_t j=(*startIndices)[k]; j<=(*stopIndices)[k]; j++) {
            dYdt(k,j) = ydotk[i];
            i++;
        }
    }
}

void ConvectionSystemSplit::setSplitConstU(const dvector& constU)
{
    utwSystem.splitConstU = constU;
}

void ConvectionSystemSplit::setSplitConstT(const dvector& constT)
{
    utwSystem.splitConstT = constT;
}

void ConvectionSystemSplit::setSplitConstY(const Array2D& constY)
{
    for (size_t k=0; k<nSpec; k++) {
        size_t i = 0;
        for (size_t j=(*startIndices)[k]; j<=(*stopIndices)[k]; j++) {
            speciesSystems[k].splitConstY[i] = constY(k,j);
            i++;
        }
    }
}

void ConvectionSystemSplit::setSplitConst
(const dvector& constU, const dvector& constT, const Array2D& constY)
{
    setSplitConstU(constU);
    setSplitConstT(constT);
    setSplitConstY(constY);
}

void ConvectionSystemSplit::setSplitLinearU(const dvector& LinearU)
{
    utwSystem.splitLinearU = LinearU;
}

void ConvectionSystemSplit::setSplitLinearT(const dvector& LinearT)
{
    utwSystem.splitLinearT = LinearT;
}

void ConvectionSystemSplit::setSplitLinearY(const Array2D& LinearY)
{
    for (size_t k=0; k<nSpec; k++) {
        size_t i = 0;
        for (size_t j=(*startIndices)[k]; j<=(*stopIndices)[k]; j++) {
            speciesSystems[k].splitLinearY[i] = LinearY(k,j);
            i++;
        }
    }
}

void ConvectionSystemSplit::setSplitLinear
(const dvector& linearU, const dvector& linearT, const Array2D& linearY)
{
    setSplitLinearU(linearU);
    setSplitLinearT(linearT);
    setSplitLinearY(linearY);
}

void ConvectionSystemSplit::integrateToTime(const double tf)
{
    // Integrate the UTW system while storing the value of v after each timestep
    boost::shared_ptr<vecInterpolator> vInterp(new vecInterpolator());
    vInterp->insert(std::make_pair(utwSolver->tInt, utwSystem.V/utwSystem.rho));
    int cvode_flag = CV_SUCCESS;
    while (cvode_flag != CV_TSTOP_RETURN) {
        cvode_flag = utwSolver->integrateOneStep(tf);
        vInterp->insert(std::make_pair(utwSolver->tInt, utwSystem.V/utwSystem.rho));
    }

    // Integrate the species systems
    for (size_t k=0; k<nSpec; k++) {
        speciesSystems[k].vInterp = vInterp;
        speciesSolvers[k].integrateToTime(tf);
    }
}

int ConvectionSystemSplit::getNumSteps()
{
    int nSteps = utwSolver->getNumSteps();
    foreach (sundialsCVODE& solver, speciesSolvers) {
        nSteps += solver.getNumSteps();
    }
    return nSteps;
}

void ConvectionSystemSplit::unroll_y()
{
    for (size_t j=0; j<nPoints; j++) {
        T[j] = utwSolver->y[j*nVars+kEnergy];
        U[j] = utwSolver->y[j*nVars+kMomentum];
    }

    for (size_t k=0; k<nSpec; k++) {
        size_t i = 0;
        for (size_t j=(*startIndices)[k]; j<=(*stopIndices)[k]; j++) {
            Y(k,j) = speciesSolvers[k].y[i];
            i++;
        }
    }
}

void ConvectionSystemSplit::configureSolver(sundialsCVODE& solver, const size_t k)
{
    solver.setODE(&speciesSystems[k]);
    solver.setBandwidth(0,0);
    solver.reltol = reltol;
    for (size_t j=0; j<nPointsSpec[k]; j++) {
        solver.abstol[j] = abstolY;
    }
    solver.linearMultistepMethod = CV_BDF;
    solver.nonlinearSolverMethod = CV_NEWTON;

    speciesSystems[k].resize(nPointsSpec[k]);
    speciesSystems[k].startIndex = (*startIndices)[k];
    speciesSystems[k].stopIndex = (*stopIndices)[k];
    speciesSystems[k].Yleft = Yleft[k];
}
