#include "diffusionSystem.h"
#include <assert.h>

DiffusionSystem::DiffusionSystem()
    : yInf(0)
    , wallConst(0)
{
}

void DiffusionSystem::get_A(sdBandMatrix& A)
{
    // Build the matrix A describing the linear ODE
    SetToZero(A.forSundials());
    size_t N = D.rows();

    assert(mathUtils::notnan(D));
    assert(mathUtils::notnan(B));
    dvector c1(N, 0);
    dvector c2(N, 0);

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
        A(0,0) = -c0;
        A(0,1) = c0;
    } else if (grid.leftBC == BoundaryCondition::WallFlux) {
        jStart = 1;
        double c0 = B[0] * (grid.alpha + 1) / hh[0];
        double d = 0.5 * (D[0]+D[1]);
        A(0,0) = - c0 * (d / hh[0] + wallConst);
        A(0,1) = d * c0 / hh[0];
    } else  { // (leftBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[0] = y[1] in the finite difference formula.
        jStart = 2;
        A(1,1) = -c1[1]*c2[1];
        A(1,2) = c1[1]*c2[1];
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
        A(N-2,N-3) = c1[N-2]*c2[N-3];
        A(N-2,N-2) = -c1[N-2]*c2[N-3];
    }

    // Intermediate points
    for (size_t j=jStart; j<jStop; j++) {
        A(j,j-1) = c1[j]*c2[j-1];
        A(j,j) = -c1[j]*(c2[j-1] + c2[j]);
        A(j,j+1) = c1[j]*c2[j];
    }
}

void DiffusionSystem::get_C(dvec& other_c)
{
    other_c = splitConst;
    if (grid.leftBC == BoundaryCondition::WallFlux) {
        other_c[0] += B[0] * (grid.alpha + 1) / hh[0] * wallConst * yInf;
    }
    assert(mathUtils::notnan(other_c));
}

void DiffusionSystem::resize(int N_)
{
    N = N_;
    B.setConstant(N, NaN);
    D.setConstant(N, NaN);
    splitConst.setConstant(N, NaN);
}

void DiffusionSystem::resetSplitConstants()
{
    splitConst.setZero(N);
}
