#include "diffusionSystem.h"
#include <assert.h>

void DiffusionSystem::get_A(sdBandMatrix& A)
{
    // Build the matrix A describing the linear ODE
    SetToZero(A.forSundials());
    size_t N = D.size();

    assert(N == jRight - jLeft + 1);
    dvector c1(N);
    dvector c2(N);

    for (size_t j=0; j<=jRight-jLeft; j++) {
        c1[j] = 0.5*B[j]/(dlj[j+jLeft]*r[j+jLeft]);
        c2[j] = rphalf[j+jLeft]*(D[j]+D[j+1])/hh[j+jLeft];
    }

    // Left boundary value
    size_t jStart;
    if (jLeft == 0) {
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
    } else {
        // truncated domain -- boundary value is fixed
        jStart = 1;
    }

    // Right boundary value
    size_t jStop;
    if (jRight == jj) {
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
    } else {
        // truncated domain -- boundary value is fixed
        jStop = N-1;
    }

    // Intermediate points
    for (size_t j=jStart; j<jStop; j++) {
        A(j,j-1) = c1[j]*c2[j-1];
        A(j,j) = -c1[j]*(c2[j-1] + c2[j]);
        A(j,j+1) = c1[j]*c2[j];
    }

    // Contribution from splitting
    for (size_t j=0; j<N; j++) {
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

    // Changing the grid resets the region the solver expects to be operating over
    jLeft = 0;
    jRight = jj;
}

void DiffusionSystem::resize(int N_)
{
    N = N_;
    B.resize(N);
    D.resize(N);
    splitConst.resize(N);
    splitLinear.resize(N);
}

void DiffusionSystem::initialize()
{
    splitConst.assign(N, 0);
    splitLinear.assign(N, 0);
}
