#include "diffusionSystem.h"

void DiffusionSystem::get_J(sdBandMatrix& J)
{
    SetToZero(J.forSundials());
    int N = D.size();
    dvector c1(N);
    dvector c2(N);
    dvector c3(N);
    for (int j=0; j<N-1; j++) {
        c1[j] = 0.5/(dlj[j]*r[j]*B[j]);
        c2[j] = rphalf[j]*(D[j]+D[j+1])/hh[j];
    }

    // left boundary value
    int jStart;
    if (leftBC == BoundaryCondition::FixedValue) {
        jStart = 1;
    } else  { // (leftBC == BoundaryCondition::ZeroGradient)
        jStart = 2;
        J(1,1) = -c1[1]*c2[1];
        J(1,2) = c1[1]*c2[1];
    }

    // right boundary value
    int jStop;
    if (rightBC == BoundaryCondition::FixedValue) {
        jStop = N-1;
    } else { // (rightBC == BoundaryCondition::ZeroGradient)
        jStop = N-2;
        J(N-2,N-3) = c1[N-2]*c2[N-3];
        J(N-2,N-2) = -c1[N-2]*c2[N-3];
    }

    // intermediate points
    for (int j=jStart; j<jStop; j++) {
        J(j,j-1) = c1[j]*c2[j-1];
        J(j,j) = -c1[j]*(c2[j-1] + c2[j]);
        J(j,j+1) = c1[j]*c2[j];
    }
}

void DiffusionSystem::get_c(dvector& other_c)
{
    other_c.assign(C.begin(), C.end());
}
