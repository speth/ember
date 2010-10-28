#include "diffusionSystem.h"

DiffusionSystem::DiffusionSystem()
{
}

void DiffusionSystem::get_A(sdBandMatrix& A)
{
    // Build the matrix A describing the linear ODE
    SetToZero(A.forSundials());
    int N = D.size();
    dvector c1(N);
    dvector c2(N);
    for (int j=0; j<N-1; j++) {
        c1[j] = 0.5*B[j]/(dlj[j]*r[j]);
        c2[j] = rphalf[j]*(D[j]+D[j+1])/hh[j];
    }

    // Left boundary value
    int jStart;
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        jStart = 1;
    } else if (grid.leftBC == BoundaryCondition::ControlVolume) {
        jStart =  1;
        double centerVol = pow(x[1], grid.alpha+1)/(grid.alpha+1);
        double centerArea = pow(x[1],grid.alpha);
        double c0 = centerArea/centerVol * (D[0]+D[1]) / (hh[0] * B[0]);
        A(0,0) = c0;
        A(0,1) = - c0;
    } else  { // (leftBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[0] = y[1] in the finite difference formula.
        jStart = 2;
        A(1,1) = -c1[1]*c2[1];
        A(1,2) = c1[1]*c2[1];
    }

    // Right boundary value
    int jStop;
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
    for (int j=jStart; j<jStop; j++) {
        A(j,j-1) = c1[j]*c2[j-1];
        A(j,j) = -c1[j]*(c2[j-1] + c2[j]);
        A(j,j+1) = c1[j]*c2[j];
    }

    // Contribution from splitting
    for (int j=0; j<N; j++) {
        A(j,j) += splitLinear[j];
    }
}

void DiffusionSystem::get_C(dvector& other_c)
{
    other_c.assign(splitConst.begin(), splitConst.end());
}

void DiffusionSystem::setGrid(const oneDimGrid& other)
{
    GridBased::setGrid(other);
    jLeft = 0;
    jRight = jj;
}
