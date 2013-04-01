#include "diffusionSystem.h"
#include <assert.h>

DiffusionSystem::DiffusionSystem()
    : yInf(0)
    , wallConst(0)
{
}

void DiffusionSystem::get_A(dvec& a, dvec& b, dvec& c)
{
    assert(mathUtils::notnan(D));
    assert(mathUtils::notnan(B));

    for (size_t j=1; j<=N-2; j++) {
        c1[j] = 0.5*B[j]/(dlj[j]*r[j]);
        c2[j] = rphalf[j]*(D[j]+D[j+1])/hh[j];
    }

    assert(mathUtils::notnan(c1));
    assert(mathUtils::notnan(c2));

    // Left boundary value
    size_t jStart;
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        jStart = 1;
    } else if (grid.leftBC == BoundaryCondition::ControlVolume) {
        jStart =  1;
        double c0 = B[0] * (grid.alpha + 1) * (D[0]+D[1]) / (2 * hh[0] * hh[0]);
        b[0] = -c0;
        c[0] = c0;
    } else if (grid.leftBC == BoundaryCondition::WallFlux) {
        jStart = 1;
        double c0 = B[0] * (grid.alpha + 1) / hh[0];
        double d = 0.5 * (D[0]+D[1]);
        b[0] = - c0 * (d / hh[0] + wallConst);
        c[0] = d * c0 / hh[0];
    } else  { // (leftBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[0] = y[1] in the finite difference formula.
        jStart = 2;
        b[1] = -c1[1]*c2[1];
        c[1] = c1[1]*c2[1];
    }

    // Right boundary value
    size_t jStop;
    if (grid.rightBC == BoundaryCondition::FixedValue) {
        jStop = N-1;
    } else { // (rightBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[N-1] = y[N-2] in the finite difference formula.
        jStop = N-2;
        a[N-2] = c1[N-2]*c2[N-3];
        b[N-2] = -c1[N-2]*c2[N-3];
    }

    // Intermediate points
    for (size_t j=jStart; j<jStop; j++) {
        a[j] = c1[j]*c2[j-1];
        b[j] = -c1[j]*(c2[j-1] + c2[j]);
        c[j] = c1[j]*c2[j];
    }
}

void DiffusionSystem::get_k(dvec& k)
{
    assert(mathUtils::notnan(splitConst));
    k = splitConst;
    if (grid.leftBC == BoundaryCondition::WallFlux) {
        k[0] += B[0] * (grid.alpha + 1) / hh[0] * wallConst * yInf;
    }
    assert(mathUtils::notnan(k));
}

void DiffusionSystem::resize(size_t N_)
{
    N = N_;
    B.setConstant(N, NaN);
    D.setConstant(N, NaN);
    splitConst.setConstant(N, NaN);
    c1.setConstant(N, 0);
    c2.setConstant(N, 0);
}

void DiffusionSystem::resetSplitConstants()
{
    splitConst.setZero(N);
}
