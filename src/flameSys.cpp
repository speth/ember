#include "flameSys.h"
#include <sstream>
#include <vector>
#include <cmath>
#include "debugUtils.h"
#include "mathUtils.h"
#include "dataFile.h"
#include "boost/filesystem.hpp"

using namespace mathUtils;

int flameSys::f(realtype t, sdVector& y, sdVector& ydot, sdVector& res)
{
    tNow = t;
    unrollY(y);
    unrollYdot(ydot);

    int jj = nPoints-1;

    // *********************************
    // *** Physical Property Updates ***
    // *********************************

    // Update the thermodynamic state, and evaluate the thermodynamic, transport and kinetic
    // parameters. Because the transport property evaluations are particularly expensive,
    // these are only done at integrator restart points as part of the initial condition
    // calculations (inGetIC = true) or if specifically requested (options.steadyOnly = false).

    gas.setStateMass(Y,T);

    try {
        updateThermoProperties();
        if (!options.steadyOnly || inGetIC) {
            perfTimerTransportProps.start();
            updateTransportProperties();
            perfTimerTransportProps.stop();
        }
        perfTimerRxnRates.start();
        gas.getReactionRates(wDot);

        perfTimerRxnRates.stop();
    } catch (Cantera::CanteraError) {
        writeStateFile("errorOutput",true);
        throw debugException("Error evaluating thermodynamic properties");
    }

    // *****************************
    // *** Update auxiliary data ***
    // *****************************

    perfTimerResFunc.start();
    for (int j=0; j<=jj; j++) {
        qDot[j] = 0;
        for (int k=0; k<nSpec; k++) {
            qDot[j] -= wDot(k,j)*hk(k,j);
        }
    }

    for (int j=0; j<=jj; j++) {
        double sum = 0;
        for (int k=0; k<nSpec; k++) {
            sum += dYdt(k,j)/W[k];
        }
        drhodt[j] = -rho[j]*(dTdt[j]/T[j] + Wmx[j]*sum);
    }

    double a = strainRate(t);
    double dadt = dStrainRatedt(t);

    if (options.xFlameControl && !inGetIC && !inTestPreconditioner) {
        update_xStag(t, false);
    }

    // *** Calculate diffusion mass fluxes, heat flux, enthalpy flux
    for (int j=0; j<jj; j++) {
        sumcpj[j] = 0;
        for (int k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1)) * ((Y(k,j+1)-Y(k,j))/hh[j]) -
                0.5*(rhoD(k,j)*Y(k,j)/Wmx[j]+Y(k,j+1)*rhoD(k,j+1)/Wmx[j+1])*(Wmx[j+1]-Wmx[j])/hh[j]; 
            jSoret(k,j) = -0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
                * (T[j+1]-T[j])/hh[j];
        }

        qFourier[j] = -0.5*(lambda[j]+lambda[j+1])*(T[j+1]-T[j])/hh[j];
        if (inGetIC) {
            jCorr[j] = 0;
            for (int k=0; k<nSpec; k++) {
                jCorr[j] -= jFick(k,j);
            }
        }
        for (int k=0; k<nSpec; k++) {
            jFick(k,j) += 0.5*(Y(k,j)+Y(k,j+1))*jCorr[j]; // correction to ensure that sum of mass fractions equals 1
            sumcpj[j] += 0.5*(cpSpec(k,j)+cpSpec(k,j+1))/W[k]*(jFick(k,j) + jSoret(k,j));
        }
    }

    // ****************************************
    // *** Left Boundary values for U, T, Y ***
    // ****************************************

    if (grid.leftBoundaryConfig == grid.lbFixedVal) {

        // Zero gradient condition for U
        momentumUnst[0] = rho[0]*dUdt[0];
        momentumProd[0] = rho[0]*U[0]*U[0] - rhou*(dadt + a*a);
        momentumDiff[0] = momentumConv[0] = 0;

        // Fixed values for T, Y
        energyUnst[0] = dTdt[0];
        energyDiff[0] = energyConv[0] = energyProd[0] = 0;

        for (int k=0; k<nSpec; k++) {
            speciesUnst(k,0) = dYdt(k,0);
            speciesDiff(k,0) = speciesConv(k,0) = speciesProd(k,0) = 0;
        }

    } else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {

        // Zero gradient condition for U
        momentumUnst[0] = rho[0]*dUdt[0];
        momentumProd[0] = rho[0]*U[0]*U[0] - rhou*(dadt + a*a);
        momentumDiff[0] = momentumConv[0] = 0;

        // Zero gradient condition for T, Y not at centerline
        energyDiff[0] = (T[1]-T[0])/hh[0];
        energyProd[0] = energyConv[0] = energyProd[0] = 0;

        for (int k=0; k<nSpec; k++) {
            speciesDiff(k,0) = (Y(k,1)-Y(k,0))/hh[0];
            speciesUnst(k,0) = speciesConv(k,0) = speciesProd(k,0) = 0;
        }

    } else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
        // Left boundary corresponds to centerline or symmetry plane
        centerVol = pow(x[1],alpha+1)/(alpha+1);
        centerArea = pow(x[1],alpha);
        double rVzero = (rV[0] > 0) ? rV[0] : 0;

        energyUnst[0] = rho[0]*dTdt[0];
        energyDiff[0] = 2*centerArea/centerVol*qFourier[0]/cp[0];
        energyProd[0] = -qDot[0]/cp[0];
        energyConv[0] = rVzero*(T[0]-Tleft)/centerVol;

        momentumUnst[0] = rho[0]*dUdt[0];
        momentumDiff[0] = -2*centerArea/centerVol*0.5*(mu[0]+mu[1])*(U[1]-U[0])/hh[0];
        momentumProd[0] = (rho[0]*U[0]*U[0] - rhou*(dadt + a*a));
        momentumConv[0] = 0;

        for (int k=0; k<nSpec; k++) {
            speciesUnst(k,0) = rho[0]*dYdt(k,0);
            speciesDiff(k,0) = 2*centerArea/centerVol*(jFick(k,0) + jSoret(k,0));
            speciesProd(k,0) = -wDot(k,0)*W[k];
            speciesConv(k,0) = rVzero*(Y(k,0)-Yleft[k])/centerVol;
        }
    }

    // ***************************************
    // *** Intermediate points for U, T, Y ***
    // ***************************************

    for (int j=1; j<jj; j++) {
        if (options.centeredDifferences) {
            // First derivative for convective terms: centered difference
            dTdx[j] = (T[j-1]*cfm[j] + T[j]*cf[j] + T[j+1]*cfp[j]);
            dUdx[j] = (U[j-1]*cfm[j] + U[j]*cf[j] + U[j+1]*cfp[j]);
            for (int k=0; k<nSpec; k++) {
                dYdx(k,j) = (Y(k,j-1)*cfm[j] + Y(k,j)*cf[j] + Y(k,j+1)*cfp[j]);
            }
        } else {
            // Upwinded convective derivatives:
            if (V[j] > 0) {
                dUdx[j] = (U[j]-U[j-1])/hh[j-1];
                dTdx[j] = (T[j]-T[j-1])/hh[j-1];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j)-Y(k,j-1))/hh[j-1];
                }
            } else {
                dUdx[j] = (U[j+1]-U[j])/hh[j];
                dTdx[j] = (T[j+1]-T[j])/hh[j];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j+1)-Y(k,j))/hh[j];
                }
            }
        }

        // The enthalpy flux term always uses centered differences
        dTdxCen[j] = T[j-1]*cfm[j] + T[j]*cf[j] + T[j+1]*cfp[j];

        // *** Momentum Equation
        momentumUnst[j] = rho[j]*dUdt[j];
        momentumConv[j] = V[j]*dUdx[j];
        momentumDiff[j] = -0.5*( rphalf[j]*(mu[j]+mu[j+1])*(U[j+1]-U[j])/hh[j] -
            rphalf[j-1]*(mu[j-1]+mu[j])*(U[j]-U[j-1])/hh[j-1])/(dlj[j]*r[j]);
        momentumProd[j] = rho[j]*U[j]*U[j] - rhou*(dadt + a*a);

        // *** Energy Equation
        energyUnst[j] = rho[j]*dTdt[j];
        energyConv[j] = V[j]*dTdx[j];
        energyDiff[j] = 0.5*(sumcpj[j]+sumcpj[j-1])*dTdxCen[j]/cp[j] +
            (rphalf[j]*qFourier[j] - rphalf[j-1]*qFourier[j-1])/(dlj[j]*cp[j]*r[j]);
        energyProd[j] = -qDot[j]/cp[j];

        // *** Species Equations
        for (int k=0; k<nSpec; k++) {
            speciesUnst(k,j) = rho[j]*dYdt(k,j);
            speciesConv(k,j) = V[j]*dYdx(k,j);
            speciesDiff(k,j) = (rphalf[j]*(jFick(k,j)+jSoret(k,j))
                - rphalf[j-1]*(jFick(k,j-1)+jSoret(k,j-1)))/(dlj[j]*r[j]);
            speciesProd(k,j) = -wDot(k,j)*W[k];
        }
    }

    // *****************************
    // *** Right boundary values ***
    // *****************************

    // *** Right boundary values for T and Y
    if (grid.unburnedLeft && !grid.fixedBurnedVal) {
        // zero gradient condition
        energyDiff[jj] = (T[jj]-T[jj-1])/hh[jj-1];
        energyProd[jj] = energyConv[jj] = energyUnst[jj] = 0;

        for (int k=0; k<nSpec; k++) {
            speciesDiff(k,jj) = (Y(k,jj)-Y(k,jj-1))/hh[jj-1];
            speciesProd(k,jj) = speciesConv(k,jj) = speciesUnst(k,jj) = 0;
        }
    } else {
        // Fixed values
        energyUnst[jj] = dTdt[jj];
        energyDiff[jj] = energyProd[jj] = energyConv[jj] = 0;

        for (int k=0; k<nSpec; k++) {
            speciesUnst(k,jj) = dYdt(k,jj);
            speciesProd(k,jj) = speciesConv(k,jj) = speciesDiff(k,jj) = 0;
        }
    }

    // Right boundary condition for U is zero gradient
    momentumUnst[jj] = rho[jj]*dUdt[jj];
    momentumProd[jj] = rho[jj]*U[jj]*U[jj] - rhou*(dadt + a*a);
    momentumDiff[jj] = momentumConv[jj] = 0;

    // ***************************
    // *** Continuity Equation ***
    // ***************************

    // *** Left boundary value
    if (options.xStagControl) {
        // Boundary value for V depends on rStag
        if (alpha == 1) {
            continuityUnst[0] = rV[0] - rVzero;
        } else {
            continuityUnst[0] = rV[0] - rVzero;
        }

    } else {
        // Boundary value for V is fixed
        continuityUnst[0] = dVdt[0];
    }

    continuityRhov[0] = continuityStrain[0] = 0;

    // *** Intermediate points
    for (int j=1; j<=jj; j++) {
        continuityRhov[j] = (rV[j]-rV[j-1])/(hh[j-1]*rphalf[j-1]);
        continuityUnst[j] = drhodt[j];
        continuityStrain[j] = rho[j]*U[j];
    }

    // *** Combine the residual components
    rollResiduals(res);
    perfTimerResFunc.stop();
    return 0;
}

int flameSys::preconditionerSetup(realtype t, sdVector& y, sdVector& ydot,
                                          sdVector& res, realtype c_j)
{
    perfTimerPrecondSetup.start();
    unrollY(y);
    unrollYdot(ydot);
    int j;
    int jj = nPoints-1;
    
    // The constant "800" here has been empirically determined to give
    // good performance for typical test cases. This value can have
    // a substantial impact on the convergence rate of the solver.
    double eps = sqrt(DBL_EPSILON)*800;

    SetToZero(bandedJacobian->forSundials());

    // *** Derivatives of reaction rate with respect to temperature
    dvector TplusdT(nPoints);

    for (j=0; j<nPoints; j++) {
        TplusdT[j] = T[j]*(1+eps);
    }

    dvector dwhdT(nPoints,0);
    Array2D dwdT(nSpec,nPoints);
    Array2D wDot2(nSpec,nPoints);

    perfTimerPrecondSetup.stop();
    perfTimerRxnRates.start();

    gas.setStateMass(Y,TplusdT);
    gas.getReactionRates(wDot2);

    perfTimerRxnRates.stop();
    perfTimerPrecondSetup.resume();

    for (j=0; j<nPoints; j++) {
        for (int k=0; k<nSpec; k++) {
            dwdT(k,j) = (wDot2(k,j)-wDot(k,j))/(TplusdT[j]-T[j]);
            dwhdT[j] += hk(k,j)*dwdT(k,j) + cpSpec(k,j)*wDot(k,j);
        }
    }

    // *** Derivatives of reaction rate with respect to species
    vector<Array2D> dwdY(nPoints, Array2D(nSpec,nSpec));
    Array2D hdwdY(nSpec,nPoints,0);
    Array2D YplusdY(nSpec,nPoints);
    for (int k=0; k<nSpec; k++) {
        YplusdY = Y;
        for (j=0; j<nPoints; j++) {
            YplusdY(k,j) = (abs(Y(k,j)) > eps/2) ? Y(k,j)*(1+eps) : eps;
        }
        perfTimerPrecondSetup.stop();
        perfTimerRxnRates.start();

        gas.setStateMass(YplusdY,T);
        gas.getReactionRates(wDot2);

        perfTimerRxnRates.stop();
        perfTimerPrecondSetup.resume();

        for (j=0; j<nPoints; j++) {
            for (int i=0; i<nSpec; i++) {
                dwdY[j](i,k) = (wDot2(i,j)-wDot(i,j))/(YplusdY(k,j)-Y(k,j));
                hdwdY(k,j) += hk(i,j)*dwdY[j](i,k);
            }
        }
    }
    // ********************************
    // *** J = dF/dy + c_j*dF/dydot ***
    // ********************************

    // *** Left Boundary values for U, T, Y
    j=0;
    if (grid.leftBoundaryConfig == grid.lbFixedVal) {
        // Fixed values for T, Y
        jacB(j,kEnergy,kEnergy) = c_j;

        // Momentum equation is always zero gradient
        jacB(j,kMomentum,kMomentum) = c_j*rho[0] + 2*rho[0]*U[0];
        jacB(j,kMomentum,kEnergy) = -rho[0]/T[0]*(U[0]*U[0] + dUdt[0]);
        for (int k=0; k<nSpec; k++) {
            jacB(j,kMomentum,kSpecies+k) = rho[0]*(W[k]-Wmx[0])/(W[k]*(1-Y(k,0)))*(U[0]*U[0] + dUdt[0]);
        }

        for (int k=0; k<nSpec; k++) {
            jacB(j,kSpecies+k,kSpecies+k) = c_j;
        }

    } else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
        // Zero gradient condition for T, Y, U, not at centerline
        jacB(j,kEnergy,kEnergy) = -1/hh[0];
        jacC(j,kEnergy,kEnergy) = 1/hh[0];

        // Momentum equation is always zero gradient
        jacB(j,kMomentum,kMomentum) = c_j*rho[0] + 2*rho[0]*U[0];
        jacB(j,kMomentum,kEnergy) = -rho[0]/T[0]*(U[0]*U[0] + dUdt[0]);
        for (int k=0; k<nSpec; k++) {
            jacB(j,kMomentum,kSpecies+k) = rho[0]*(W[k]-Wmx[0])/(W[k]*(1-Y(k,0)))*(U[0]*U[0] + dUdt[0]);
        }

        for (int k=0; k<nSpec; k++) {
            jacB(j,kSpecies+k,kSpecies+k) = -1/hh[0];
            jacC(j,kSpecies+k,kSpecies+k) = 1/hh[0];
        }

    } else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
        // Left boundary corresponds to centerline or symmetry plane
        centerVol = pow(x[1],alpha+1)/(alpha+1);
        centerArea = pow(x[1],alpha);
        double rVzero, drVzero;
        if (rV[0] > 0) {
            rVzero = rV[0];
            drVzero = 1;
        } else {
            rVzero = 0;
            drVzero = 0;
        }

        // dEnergy/dT_j
        jacB(j,kEnergy,kEnergy) = c_j*rho[j] - rho[j]/T[j]*dTdt[j] + dwhdT[j]/cp[j]
            + 2*centerArea/centerVol*0.5*(lambda[j]+lambda[j+1])/(hh[j]*cp[j]) + rVzero/centerVol;

        // dEnergy/dT_(j+1)
        jacC(j,kEnergy,kEnergy) = - 2*centerArea/centerVol*0.5*(lambda[j]+lambda[j+1])/(hh[j]*cp[j]);

        // dEnergy/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kEnergy,kSpecies+k) = rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j)))*dTdt[j] + hdwdY(k,j)/cp[j];
        }

        // dEnergy/dV
        jacB(j,kEnergy,kContinuity) = drVzero*(T[0]-Tleft)/centerVol;

        // dMomentum/dT
        jacB(j,kMomentum,kEnergy) = -rho[j]/T[j]*(dUdt[j] + U[j]*U[j]);

        // dMomentum/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kMomentum,kSpecies+k) = rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j))) * (dUdt[j] + U[j]*U[j]);
        }

        // dMomentum/dU_j
        jacB(j,kMomentum,kMomentum) = rho[j]*c_j + 2*rho[j]*U[j]
            + 2*centerArea/centerVol*0.5*(mu[j]+mu[j+1])/hh[j];

        // dMomentum/dU_j+1
        jacC(j,kMomentum,kMomentum) = -2*centerArea/centerVol*0.5*(mu[j]+mu[j+1])/hh[j];

        for (int k=0; k<nSpec; k++) {
            // dSpecies/dT
            jacB(j,kSpecies+k,kEnergy) = - rho[j]/T[j]*dYdt(k,j) - dwdT(k,j)*W[k];

            // dSpecies/dY_j
            for (int i=0; i<nSpec; i++) {
                jacB(j,kSpecies+k,kSpecies+i) = rho[j]*(W[i]-Wmx[j])/(W[i]*(1-Y(i,j)))*dYdt(k,j) - dwdY[j](k,i)*W[k];
            }
            jacB(j,kSpecies+k,kSpecies+k) += rho[j]*c_j
                + 2*centerArea/centerVol*0.5*(rhoD(k,j)+rhoD(k,j+1))/hh[j] + rVzero/centerVol;

            // dSpecies/dY_(j+1)
            jacC(j,kSpecies+k,kSpecies+k) = -2*centerArea/centerVol*0.5*(rhoD(k,j)+rhoD(k,j+1))/hh[j];

            // dSpecies/dV
            jacB(j,kSpecies+k,kContinuity) = drVzero*(Y(k,0)-Yleft[k])/centerVol;
        }

    }

    // *** Intermediate points for U, T, Y
    for (j=1; j<jj; j++) {

        if (options.centeredDifferences) {
            // First derivative for convective terms: centered difference
            dTdx[j] = (T[j-1]*cfm[j] + T[j]*cf[j] + T[j+1]*cfp[j]);
            dUdx[j] = (U[j-1]*cfm[j] + U[j]*cf[j] + U[j+1]*cfp[j]);
            for (int k=0; k<nSpec; k++) {
                dYdx(k,j) = (Y(k,j-1)*cfm[j] + Y(k,j)*cf[j] + Y(k,j+1)*cfp[j]);
            }
        } else {
            // Upwinded convective derivatives:
            if (V[j] > 0) {
                dUdx[j] = (U[j]-U[j-1])/hh[j-1];
                dTdx[j] = (T[j]-T[j-1])/hh[j-1];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j)-Y(k,j-1))/hh[j-1];
                }
            } else {
                dUdx[j] = (U[j+1]-U[j])/hh[j];
                dTdx[j] = (T[j+1]-T[j])/hh[j];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j+1)-Y(k,j))/hh[j];
                }
            }
        }

        dTdxCen[j] = T[j-1]*cfm[j] + T[j]*cf[j] + T[j+1]*cfp[j];

        for (int k=0; k<nSpec; k++) {
            // dSpecies/dY
            for (int i=0; i<nSpec; i++) {
                jacB(j,kSpecies+k,kSpecies+i) = -dwdY[j](k,i)*W[k] +
                    rho[j]*(W[i]-Wmx[j])/(W[i]*(1-Y(i,j)))*dYdt(k,j);
            }

            // dSpecies/dT
            jacB(j,kSpecies+k,kEnergy) = -dwdT(k,j)*W[k] - rho[j]/T[j]*dYdt(k,j);

            // dSpecies/dV
            jacB(j,kSpecies+k,kContinuity) = dYdx(k,j);

            // dSpeciesk/dYk_j
            jacB(j,kSpecies+k,kSpecies+k) += rho[j]*c_j  +
                0.5*rphalf[j]*((rhoD(k,j)+rhoD(k,j+1))/hh[j]+jCorr[j])/(dlj[j]*r[j]) +
                0.5*rphalf[j-1]*((rhoD(k,j-1)+rhoD(k,j))/hh[j-1]-jCorr[j-1])/(dlj[j]*r[j]);

            jacB(j,kSpecies+k,kSpecies+k) += 0.5*rhoD(k,j)/Wmx[j]*(-rphalf[j]*(Wmx[j+1]-Wmx[j])/hh[j] + rphalf[j-1]*(Wmx[j]-Wmx[j-1])/hh[j-1])/(dlj[j]*r[j]); 

            // dSpeciesk/dYk_j-1
            jacA(j,kSpecies+k,kSpecies+k) = -0.5*rphalf[j-1]
                * ((rhoD(k,j-1)+rhoD(k,j))/hh[j-1]+jCorr[j-1])
                / (dlj[j]*r[j]);
                
            jacA(j,kSpecies+k,kSpecies+k) += 0.5*rphalf[j-1]*rhoD(k,j-1)/Wmx[j-1]*(Wmx[j]-Wmx[j-1])/(hh[j]*dlj[j]*r[j]);

            // dSpeciesk/dYk_j+1
            jacC(j,kSpecies+k,kSpecies+k) = -0.5*rphalf[j]
                * ((rhoD(k,j)+rhoD(k,j+1))/hh[j]-jCorr[j])
                / (dlj[j]*r[j]);
                
            jacC(j,kSpecies+k,kSpecies+k) -= 0.5*rphalf[j]*rhoD(k,j+1)/Wmx[j+1]*(Wmx[j+1]-Wmx[j])/(hh[j]*dlj[j]*r[j]);

            if (options.centeredDifferences) {
                jacA(j,kSpecies+k,kSpecies+k) += V[j]*cfm[j];
                jacB(j,kSpecies+k,kSpecies+k) += V[j]*cf[j];
                jacC(j,kSpecies+k,kSpecies+k) += V[j]*cfp[j];
            } else {
                if (V[j] > 0) {
                    jacB(j,kSpecies+k,kSpecies+k) += V[j]/hh[j-1];
                    jacA(j,kSpecies+k,kSpecies+k) -= V[j]/hh[j-1];
                } else {
                    jacB(j,kSpecies+k,kSpecies+k) -= V[j]/hh[j];
                    jacC(j,kSpecies+k,kSpecies+k) += V[j]/hh[j];
                }
            }
        }

        // dEnergy/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kEnergy,kSpecies+k) = hdwdY(k,j)/cp[j] + rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j)))*dTdt[j];

            // Enthalpy flux term
            jacA(j,kEnergy,kSpecies+k) += 0.25*dTdxCen[j]*(cpSpec(k,j)+cpSpec(k,j-1))/W[k]
                *(rhoD(k,j)+rhoD(k,j-1))/cp[j]/hh[j-1];
            jacB(j,kEnergy,kSpecies+k) -= 0.25*dTdxCen[j]*(cpSpec(k,j)+cpSpec(k,j-1))/W[k]
                *(rhoD(k,j)+rhoD(k,j-1))/cp[j]/hh[j-1];
            jacB(j,kEnergy,kSpecies+k) += 0.25*dTdxCen[j]*(cpSpec(k,j)+cpSpec(k,j+1))/W[k]
                *(rhoD(k,j)+rhoD(k,j+1))/cp[j]/hh[j];
            jacC(j,kEnergy,kSpecies+k) -= 0.25*dTdxCen[j]*(cpSpec(k,j)+cpSpec(k,j+1))/W[k]
                *(rhoD(k,j)+rhoD(k,j+1))/cp[j]/hh[j];
        }

        // dEnergy/dT_j
        jacB(j,kEnergy,kEnergy) = rho[j]*c_j - rho[j]/T[j]*dTdt[j] +
            (0.5*rphalf[j]*(lambda[j]+lambda[j+1])/hh[j] +
            0.5*rphalf[j-1]*(lambda[j-1]+lambda[j])/hh[j-1])/(dlj[j]*cp[j]*r[j]) +
            dwhdT[j]/cp[j] + 0.5*(sumcpj[j]+sumcpj[j-1])*cf[j]/cp[j];


        // dEnergy/dT_j-1
        jacA(j,kEnergy,kEnergy) = -0.5*rphalf[j-1]*(lambda[j-1]+lambda[j])/(hh[j-1]*dlj[j]*cp[j]*r[j]) +
            0.5*(sumcpj[j]+sumcpj[j+1])*cfm[j]/cp[j];

        // dEnergy/dT_j+1
        jacC(j,kEnergy,kEnergy) = -0.5*rphalf[j]*(lambda[j]+lambda[j+1])/(hh[j]*dlj[j]*cp[j]*r[j]) +
            0.5*(sumcpj[j]+sumcpj[j-1])*cfp[j]/cp[j];

        if (options.centeredDifferences) {
            jacA(j,kEnergy,kEnergy) += V[j]*cfm[j];
            jacB(j,kEnergy,kEnergy) += V[j]*cf[j];
            jacC(j,kEnergy,kEnergy) += V[j]*cfp[j];
        } else {
            if (V[j] > 0) {
                jacB(j,kEnergy,kEnergy) += V[j]/hh[j-1];
                jacA(j,kEnergy,kEnergy) -= V[j]/hh[j-1];
            } else {
                jacB(j,kEnergy,kEnergy) -= V[j]/hh[j];
                jacC(j,kEnergy,kEnergy) += V[j]/hh[j];
            }
        }

        // dEnergy/dV
        jacB(j,kEnergy,kContinuity) = dTdx[j];

        // dMomentum/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kMomentum,kSpecies+k) =  rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j))) *
                (dUdt[j] + U[j]*U[j]);
        }

        // dMomentum/dT
        jacB(j,kMomentum,kEnergy) = -rho[j]/T[j] * (dUdt[j] + U[j]*U[j]);

        // dMomentum/dU_j
        jacB(j,kMomentum,kMomentum) = rho[j]*c_j + 2*U[j]*rho[j] +
            0.5*(rphalf[j]*(mu[j]+mu[j+1])/hh[j] + rphalf[j]*(mu[j-1]+mu[j])/hh[j-1])/(dlj[j]*r[j]);

        // dMomentum/dU_j-1
        jacA(j,kMomentum,kMomentum) = -0.5*rphalf[j]*(mu[j-1]+mu[j])/(hh[j-1]*dlj[j]*r[j]);

        // dMomentum/dU_j+1
        jacC(j,kMomentum,kMomentum) = -0.5*rphalf[j]*(mu[j]+mu[j+1])/(hh[j]*dlj[j]*r[j]);

        if (options.centeredDifferences) {
            jacA(j,kMomentum,kMomentum) += V[j]*cfm[j];
            jacB(j,kMomentum,kMomentum) += V[j]*cf[j];
            jacC(j,kMomentum,kMomentum) += V[j]*cfp[j];
        } else {
            if (V[j] > 0) {
                jacB(j,kMomentum,kMomentum) += V[j]/hh[j-1];
                jacA(j,kMomentum,kMomentum) -= V[j]/hh[j-1];
            } else {
                jacB(j,kMomentum,kMomentum) -= V[j]/hh[j];
                jacC(j,kMomentum,kMomentum) += V[j]/hh[j];
            }
        }

        // dMomentum/dV
        jacB(j,kMomentum,kContinuity) = dUdx[j];
    }

    // *** Right boundary values for T, Y
    j = jj;
    if (grid.unburnedLeft && !grid.fixedBurnedVal) {
        // zero gradient condition
        jacA(j,kEnergy,kEnergy) = -1/hh[j-1];
        jacB(j,kEnergy,kEnergy) = 1/hh[j-1];
        for (int k=0; k<nSpec; k++) {
            jacA(j,kSpecies+k,kSpecies+k) = -1/hh[j-1];
            jacB(j,kSpecies+k,kSpecies+k) = 1/hh[j-1];
        }
    } else {
        // Fixed values
        jacB(j,kEnergy,kEnergy) = c_j;
        for (int k=0; k<nSpec; k++) {
            jacB(j,kSpecies+k,kSpecies+k) = c_j;
        }
    }

    // *** Right boundary values for U (zero gradient)
    jacB(j,kMomentum,kMomentum) = rho[jj]*c_j + 2*U[jj]*rho[jj];
    jacB(j,kMomentum,kEnergy) = -rho[jj]/T[jj]*(U[jj]*U[jj] + dUdt[jj]);
    for (int k=0; k<nSpec; k++) {
        jacB(j,kMomentum,kSpecies+k) = rho[jj]*(W[k]-Wmx[jj])/(W[k]*(1-Y(k,jj)))*(U[jj]*U[jj] + dUdt[jj]);
    }

    // *** Continuity Equation
    if (options.xStagControl) {
        jacB(0,kContinuity,kContinuity) = pow(x[1],alpha);
    } else {
        jacB(0,kContinuity,kContinuity) = c_j;
    }

    for (j=1; j<=jj; j++) {
        double sumdYkWk = 0;
        for (int k=0; k<nSpec; k++) {
            sumdYkWk += dYdt(k,j)/W[k];
        }

        for (int k=0; k<nSpec; k++) {
            // dContinuity/dY
            jacB(j,kContinuity,kSpecies+k) = rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j)))*(U[j] - dTdt[j]/T[j] - 2*Wmx[j]*sumdYkWk) -
                c_j*rho[j]*Wmx[j]/W[k];
        }

        // dContinuity/dT
        jacB(j,kContinuity,kEnergy) = -rho[j]/T[j]*(c_j + U[j]);

        // dContinuity/dU
        jacB(j,kContinuity,kMomentum) = rho[j];

        // dContinuity/dV
        jacB(j,kContinuity,kContinuity) = r[j]/hh[j-1]/rphalf[j-1];

        // dContinuity/dV_(j-1)
        jacA(j,kContinuity,kContinuity) = - r[j-1]/hh[j-1]/rphalf[j-1];
    }

    perfTimerPrecondSetup.stop();

    // *** Get LU Factorization of the Jacobian
    if (!inTestPreconditioner) {
        perfTimerLU.start();
        int iError = BandGBTRF(bandedJacobian->forSundials(),&pMat[0]);
        perfTimerLU.stop();

        if (iError!=0) {
            cout << "Error in LU factorization: i = " << iError << endl;
            debugParameters::debugJacobian = true;
            return 1;
        } else {
            debugParameters::debugJacobian = false;
        }
    }
    return 0;
}


int flameSys::preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn,
                                          sdVector& resIn, sdVector& rhs,
                                          sdVector& outVec, realtype c_j, realtype delta)
{
    // Calculate J^-1*v
    perfTimerPrecondSolve.start();

    double* x = N_VGetArrayPointer(outVec.forSundials());
    for (int i=0; i<N; i++) {
        x[i] = rhs(i);
    }

    BandGBTRS(bandedJacobian->forSundials(),&pMat[0],x);

    perfTimerPrecondSolve.stop();

    return 0;
}

void flameSys::setup(void)
{
    perfTimerSetup.start();
    int nPointsOld = T.size();

    if (nPoints == nPointsOld) {
        return; // nothing to do
    }

    gas.resize(nPoints);
    nSpec = gas.nSpec;

    nVars = 3+nSpec;
    N = nVars*nPoints;

    V.resize(nPoints);
    U.resize(nPoints);
    T.resize(nPoints);
    Y.resize(nSpec,nPoints);

    dVdt.resize(nPoints,0);
    dUdt.resize(nPoints,0);
    dTdt.resize(nPoints,0);
    dYdt.resize(nSpec,nPoints);

    dUdx.resize(nPoints,0);
    dTdx.resize(nPoints,0);
    dYdx.resize(nSpec,nPoints,0);
    dTdxCen.resize(nPoints,0);

    energyUnst.resize(nPoints,0); energyDiff.resize(nPoints,0);
    energyConv.resize(nPoints,0); energyProd.resize(nPoints,0);
    momentumUnst.resize(nPoints,0); momentumDiff.resize(nPoints,0);
    momentumConv.resize(nPoints,0); momentumProd.resize(nPoints,0);
    speciesUnst.resize(nSpec,nPoints,0); speciesDiff.resize(nSpec,nPoints,0);
    speciesConv.resize(nSpec,nPoints,0); speciesProd.resize(nSpec,nPoints,0);
    continuityUnst.resize(nPoints,0); continuityRhov.resize(nPoints,0);
    continuityStrain.resize(nPoints,0);

    rV.resize(nPoints);
    rho.resize(nPoints);
    drhodt.resize(nPoints);
    Wmx.resize(nPoints);
    W.resize(nSpec);

    mu.resize(nPoints);
    lambda.resize(nPoints);
    cp.resize(nPoints);
    cpSpec.resize(nSpec,nPoints);
    qFourier.resize(nPoints);

    rhoD.resize(nSpec,nPoints);
    Dkt.resize(nSpec,nPoints);
    sumcpj.resize(nPoints);

    jFick.resize(nSpec,nPoints);
    jSoret.resize(nSpec,nPoints);
    jCorr.resize(nPoints);
    wDot.resize(nSpec,nPoints);
    hk.resize(nSpec,nPoints);
    qDot.resize(nPoints);

    delete bandedJacobian;
    bandedJacobian = new sdBandMatrix(N,2*nVars+1,2*nVars+1,4*nVars+2);
    SetToZero(bandedJacobian->forSundials());

    pMat.resize(N);
    grid.jj = nPoints-1;
    grid.updateBoundaryIndices();
    perfTimerSetup.stop();
}

void flameSys::copyOptions(void)
{
    strainRateInitial = options.strainRateInitial;
    strainRateFinal = options.strainRateFinal;
    strainRateT0 = options.strainRateT0;
    strainRateDt = options.strainRateDt;

    tStart = options.tStart;
    tEnd = options.tEnd;

    gas.mechanismFile = options.gasMechanismFile;
    gas.phaseID = options.gasPhaseID;
    gas.pressure = options.pressure;

    nPoints = options.nPoints;
    grid.unburnedLeft = options.unburnedLeft;

    grid.absvtol = options.absvtol;
    grid.rmTol = options.rmTol;
    grid.uniformityTol = options.uniformityTol;
    grid.gridMin = options.gridMin;
    grid.centerGridMin = options.centerGridMin;
    grid.gridMax = options.gridMax;
    grid.dampConst = options.dampConst;
    grid.fixedBurnedVal = options.fixedBurnedVal;
    grid.unburnedLeft = options.unburnedLeft;
    grid.boundaryTol = options.boundaryTol;
    grid.boundaryTolRm = options.boundaryTolRm;
    grid.addPointCount = options.addPointCount;

    alpha = grid.alpha = options.gridAlpha;
    kMomentum = grid.kMomentum = options.kMomentum;
    kContinuity = grid.kContinuity = options.kContinuity;
    kEnergy = grid.kEnergy = options.kEnergy;
    kSpecies = grid.kSpecies = options.kSpecies;
}

void flameSys::generateInitialProfiles(void)
{
    cout << "Generating initial profiles from given fuel and oxidizer compositions." << endl;
    x.resize(nPoints);
    int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.
    int jl = jm - 4;
    int jr = jm + 4;
    Yb.resize(nSpec); Yu.resize(nSpec);
    double a = strainRate(tStart);

    Tu = options.Tu;

    // Reactants
    gas.thermo.setState_TPX(Tu,gas.pressure,&options.reactants[0]);
    rhou = gas.thermo.density();
    gas.thermo.getMassFractions(&Yu[0]);
    gas.thermo.getMassFractions(&Y(0,grid.ju));

    // Products
    gas.thermo.setState_TPX(Tu,gas.pressure,&options.reactants[0]);
    Cantera::equilibrate(gas.thermo,"HP");
    Tb = gas.thermo.temperature();
    rhob = gas.thermo.density();
    gas.thermo.getMassFractions(&Yb[0]);
    gas.thermo.getMassFractions(&Y(0,grid.jb));

    // Diluent in the middle
    gas.thermo.setState_TPY(Tu,gas.pressure,options.oxidizer);
    gas.thermo.getMassFractions(&Y(0,jm));

    if (options.unburnedLeft) {
        rhoLeft = rhou;
        Tleft = Tu;
        Yleft = Yu;
        rhoRight = rhob;
        Tright = Tb;
        Yright = Yb;
    } else {
        rhoLeft = rhob;
        Tleft = Tb;
        Yleft = Yb;
        rhoRight = rhou;
        Tright = Tu;
        Yright = Yu;
    }

    T[0] = Tleft; T[grid.jj] = Tright;
    T[jm] = T[grid.ju];

    xRight = options.xRight;
    if (options.twinFlame || options.curvedFlame) {
        x[0] = 0;
        xLeft = options.centerGridMin;
        for (int j=1; j<nPoints; j++) {
            x[j] = xLeft + (xRight-xLeft)*((double) j)/((double) nPoints-1);
        }
    } else {
        // Uniform initial grid
        xLeft = options.xLeft;
        for (int j=0; j<nPoints; j++) {
            x[j] = xLeft + (xRight-xLeft)*((double) j)/((double) nPoints);
        }

    }

    for (int j=1; j<jl; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0);
        }
        T[j] = T[0];
    }

    for (int j=jl; j<jm; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0) + (Y(k,jm)-Y(k,0))*(x[j]-x[jl])/(x[jm]-x[jl]);
        }
        T[j] = T[0] + (T[jm]-T[0])*(x[j]-x[jl])/(x[jm]-x[jl]);
    }

    for (int j=jm+1; j<jr; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,jm) + (Y(k,grid.jj)-Y(k,jm))*(x[j]-x[jm])/(x[jr]-x[jm]);
        }
        T[j] = T[jm] + (T[grid.jj]-T[jm])*(x[j]-x[jm])/(x[jr]-x[jm]);
    }

    for (int j=jr; j<nPoints; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,grid.jj);
        }
        T[j] = T[grid.jj];
    }

    dvector yTemp(nPoints);
    for (int k=0; k<nSpec; k++) {
        for (int j=0; j<nPoints; j++) {
            yTemp[j] = Y(k,j);
        }

        for (int i=0; i<10; i++) {
            mathUtils::smooth(yTemp);
        }

        for (int j=0; j<nPoints; j++) {
            Y(k,j) = yTemp[j];
        }
    }

    for (int i=0; i<5; i++) {
        mathUtils::smooth(T);
    }

    // Grid and initial profiles of T, U and V
    for (int j=0; j<nPoints; j++) {
        gas.thermo.setState_TPY(T[j],gas.pressure,&Y(0,j));
        rho[j] = gas.thermo.density();
        U[j] = a*sqrt(rhou/rho[j]);
    }

    for (int i=0; i<2; i++) {
        mathUtils::smooth(U);
    }

    if (options.fixedLeftLoc)
    {
        jm = 0;
    }

    V[jm] = 0;
    for (int j=jm+1; j<nPoints; j++) {
        V[j] = V[j-1] - rho[j]*U[j]*(x[j]-x[j-1]);
    }

    for (int j=jm-1; j>=0; j--) {
        V[j] = V[j+1] + rho[j]*U[j]*(x[j+1]-x[j]);
    }

    gas.setStateMass(Y,T);
    gas.getMolecularWeights(W);
}

void flameSys::loadInitialProfiles(void)
{
    std::string inputFilename;
    if (options.useRelativeRestartPath) {
        inputFilename = options.inputDir + "/" + options.restartFile;
    } else {
        inputFilename = options.restartFile;
    }

    cout << "Reading initial condition from " << inputFilename << endl;
    DataFile infile(inputFilename);
    x = infile.readVector("x");

    nPoints = x.size();
    setup();

    U = infile.readVector("U");
    V = infile.readVector("V");
    T = infile.readVector("T");
    Y = infile.readArray2D("Y");
    tStart = infile.readScalar("t");

    if (!options.fileNumberOverride) {
        options.outputFileNumber = (int) infile.readScalar("fileNumber");
    }

    infile.close();

    if (options.overrideTu) {
        T[grid.ju] = options.Tu;
    }
    Tu = T[grid.ju];

    // save the burned gas properties for the case where burned values are not fixed
    double TbSave = T[grid.jb];
    dvector YbSave(nSpec);
    for (int k=0; k<nSpec; k++) {
        YbSave[k] = Y(k,grid.jb);
    }

    if (options.overrideReactants) {
        gas.thermo.setState_TPX(Tu,gas.pressure,&options.reactants[0]);
        gas.thermo.getMassFractions(&Y(0,grid.ju));
        gas.thermo.setState_TPX(Tu,gas.pressure,&options.reactants[0]);
        Cantera::equilibrate(gas.thermo,"HP");
        gas.thermo.getMassFractions(&Y(0,grid.jb));
        T[grid.jb] = gas.thermo.temperature();
    }
    Tb = T[grid.jb];

    gas.setStateMass(Y,T);
    gas.getMolecularWeights(W);
    gas.getDensity(rho);
    rhou = rho[grid.ju];
    rhob = rho[grid.jb];
    Yu.resize(nSpec); Yb.resize(nSpec);
    for (int k=0; k<nSpec; k++) {
        Yu[k] = Y(k,grid.ju);
        Yb[k] = Y(k,grid.jb);
    }

    if (options.unburnedLeft) {
        rhoLeft = rhou;
        Tleft = Tu;
        Yleft = Yu;
        rhoRight = rhob;
        Tright = Tb;
        Yright = Yb;
    } else {
        rhoLeft = rhob;
        Tleft = Tb;
        Yleft = Yb;
        rhoRight = rhou;
        Tright = Tu;
        Yright = Yu;
    }

    if (!options.fixedBurnedVal) {
        T[grid.jb] = TbSave;
        for (int k=0; k<nSpec; k++) {
            Y(k,grid.jb) = YbSave[k];
        }
    }

    V2rV();
    updateLeftBC();

    double controlSignal;
    if (grid.leftBoundaryConfig == grid.lbControlVolume && options.xFlameControl) {
        if (alpha == 0) {
            controlSignal = V[0]/rhoLeft;
        } else {
            double tmp = pow(x[0],2) + 2*rV[0]/rhoLeft;
            controlSignal = sign(tmp)*sqrt(abs(tmp));
        }
        flamePosIntegralError = controlSignal/(options.xFlameProportionalGain*options.xFlameIntegralGain);
    }
}

flameSys::flameSys(void)
    : grid(options)
    , x(grid.x)
    , r(grid.r)
    , rphalf(grid.rphalf)
    , hh(grid.hh)
    , dlj(grid.dlj)
    , cfm(grid.cfm)
    , cf(grid.cf)
    , cfp(grid.cfp)
    , flamePosIntegralError(0)
    , bandedJacobian(NULL)
{
    inGetIC = false;
    inTestPreconditioner = false;
    tNow = 0;
    ICfileNumber = 0;
}

flameSys::~flameSys(void)
{
    delete bandedJacobian;
}

void flameSys::unrollYdot(const sdVector& yDot)
{
    for (int j=0; j<nPoints; j++) {
        dVdt[j] = yDot(nVars*j);
        dUdt[j] = yDot(nVars*j+1);
        dTdt[j] = yDot(nVars*j+2);
        for (int k=0; k<nSpec; k++) {
            dYdt(k,j) = yDot(nVars*j+3+k);
        }
    }
}

void flameSys::rollYdot(sdVector& yDot)
{
    for (int j=0; j<nPoints; j++) {
        yDot(nVars*j) = dVdt[j];
        yDot(nVars*j+1) = dUdt[j];
        yDot(nVars*j+2) = dTdt[j];
        for (int k=0; k<nSpec; k++) {
            yDot(nVars*j+3+k) = dYdt(k,j);
        }
    }
}

void flameSys::unrollY(const sdVector& y)
{
    for (int j=0; j<nPoints; j++) {
        V[j] = y(nVars*j);
        U[j] = y(nVars*j+1);
        T[j] = y(nVars*j+2);
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = y(nVars*j+3+k);
        }
    }

    V2rV();
}

void flameSys::rollY(sdVector& y)
{
    for (int j=0; j<nPoints; j++) {
        y(nVars*j) = V[j];
        y(nVars*j+1) = U[j];
        y(nVars*j+2) = T[j];
        for (int k=0; k<nSpec; k++) {
            y(nVars*j+k+3) = Y(k,j);
        }
    }
}

void flameSys::rollResiduals(sdVector& res)
{
    for (int j=0; j<nPoints; j++) {
        res(nVars*j) = continuityRhov[j] + continuityStrain[j] + continuityUnst[j];
        res(nVars*j+1) = momentumConv[j] + momentumDiff[j] + momentumProd[j] + momentumUnst[j];
        res(nVars*j+2) = energyConv[j] + energyDiff[j] + energyProd[j] + energyUnst[j];
        for (int k=0; k<nSpec; k++) {
            res(nVars*j+k+3) = speciesConv(k,j) + speciesDiff(k,j) + speciesProd(k,j) + speciesUnst(k,j);
        }
    }
}

void flameSys::rollVectorVector(const sdVector& y, const dvector& qdot, vector<dvector>& v)
{
    v.resize(nVars+1);
    for (int i=0; i<nVars; i++) {
        v[i].resize(nPoints);
        for (int j=0; j<nPoints; j++) {
            v[i][j] = y(i+nVars*j);
        }
    }
    v[nVars].resize(nPoints);
    for (int j=0; j<nPoints; j++) {
        v[nVars][j] = qdot[j];
    }
}

void flameSys::unrollVectorVector(const vector<dvector>& v)
{
    for (int j=0; j<nPoints; j++) {
        V[j] = v[kContinuity][j];
        U[j] = v[kMomentum][j];
        T[j] = v[kEnergy][j];
    }

    for (int k=0; k<nSpec; k++) {
        for (int j=0; j<nPoints; j++) {
            Y(k,j) = v[kSpecies+k][j];
        }
    }
}

void flameSys::unrollVectorVectorDot(const vector<dvector>& v)
{
    for (int j=0; j<nPoints; j++) {
        dVdt[j] = v[kContinuity][j];
        dUdt[j] = v[kMomentum][j];
        dTdt[j] = v[kEnergy][j];
    }

    for (int k=0; k<nSpec; k++) {
        for (int j=0; j<nPoints; j++) {
            dYdt(k,j) = v[kSpecies+k][j];
        }
    }
}

void flameSys::updateTransportProperties(void)
{
    gas.getTransportProperties(mu, lambda, rhoD, Dkt);
}

void flameSys::updateThermoProperties(void)
{
    gas.getThermoProperties(rho, Wmx, cp, cpSpec, hk);
}

void flameSys::printForMatlab(ofstream& file, dvector& v, int index, char* name)
{
    file << name << "{" << index << "} = [";
    for (unsigned int i=0; i<v.size()-1; i++)
    {
        file << v[i] << ", ";
    }
    file << v[v.size()-1] << "];" << endl;
}

void flameSys::V2rV(void)
{
    if (alpha == 0) {
        for (int j=0; j<nPoints; j++) {
            rV[j] = V[j];
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            rV[j] = x[j]*V[j];
        }
    }

    // In the case where the centerline is part of the domain,
    // V[0] actually stores rV[0] (since V = Inf there)
    if (x[0] == 0) {
        rV[0] = V[0];
    }
}

void flameSys::rV2V(void)
{
    if (alpha == 0) {
        for (int j=0; j<nPoints; j++) {
            V[j] = rV[j];
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            V[j] = rV[j]/x[j];
        }
    }

    // In the case where the centerline is part of the domain,
    // V[0] actually stores rV[0] (since V = Inf there)
    if (x[0] == 0) {
        V[0] = rV[0];
    }
}

void flameSys::writeStateFile(const std::string fileNameStr, bool errorFile)
{
    std::ostringstream fileName(ostringstream::out);
    bool incrementFileNumber = false;

    if (fileNameStr.length() == 0) {
        // Determine the name of the output file (outXXXXXX.h5)
        incrementFileNumber = true;
        if (errorFile) {
            fileName << options.outputDir << "/error";
        } else {
            fileName << options.outputDir << "/prof";
        }
        fileName.flags(ios_base::right);
        fileName.fill('0');
        fileName.width(6);
        fileName << options.outputFileNumber << ".h5";
    } else {
        fileName << options.outputDir << "/" << fileNameStr << ".h5";
    }
    if (errorFile) {
        cout << "Writing error output file: " << fileName.str() << endl;
    } else {
        cout << "Writing output file: " << fileName.str() << endl;
    }

    // Erase the existing file and create a new one
    if (boost::filesystem::exists(fileName.str())) {
        boost::filesystem::remove(fileName.str());
    }
    DataFile outFile(fileName.str());

    // Write the state data to the output file:
    outFile.writeScalar("t", tNow);
    outFile.writeVector("x", x);
    outFile.writeVector("T", T);
    outFile.writeVector("U", U);
    outFile.writeVector("V", V);
    outFile.writeArray2D("Y", Y);
    outFile.writeScalar("a",strainRate(tNow));
    outFile.writeScalar("dadt",dStrainRatedt(tNow));
    outFile.writeScalar("fileNumber", options.outputFileNumber);

    if (options.outputHeatReleaseRate || errorFile) {
        outFile.writeVector("q",qDot);
        outFile.writeVector("rho", rho);
    }

    if (options.outputTimeDerivatives || errorFile) {
        outFile.writeVector("dUdt", dUdt);
        outFile.writeVector("dTdt", dTdt);
        outFile.writeVector("dVdt", dVdt);
        outFile.writeArray2D("dYdt", dYdt);
        outFile.writeVector("drhodt",drhodt);
    }

    if (options.outputAuxiliaryVariables || errorFile) {
        outFile.writeArray2D("wdot", wDot);
        outFile.writeArray2D("rhoD",rhoD);
        outFile.writeVector("lambda",lambda);
        outFile.writeVector("cp",cp);
        outFile.writeVector("mu",mu);
        outFile.writeVector("Wmx",Wmx);
        outFile.writeVector("W",W);
        outFile.writeArray2D("jFick", jFick);
        outFile.writeArray2D("jSoret", jSoret);
        outFile.writeVector("qFourier",qFourier);
        outFile.writeVector("cfp",grid.cfp);
        outFile.writeVector("cf",grid.cf);
        outFile.writeVector("cfm",grid.cfm);
        outFile.writeVector("hh",hh);
        outFile.writeVector("rphalf",grid.rphalf);
        outFile.writeScalar("Tleft",Tleft);
        outFile.writeVector("Yleft",Yleft);
        outFile.writeVector("sumcpj",sumcpj);
    }

    if (options.outputResidualComponents || errorFile) {
        outFile.writeVector("resEnergyUnst",energyUnst);
        outFile.writeVector("resEnergyDiff",energyDiff);
        outFile.writeVector("resEnergyConv",energyConv);
        outFile.writeVector("resEnergyProd",energyProd);
        outFile.writeVector("resMomentumUnst",momentumUnst);
        outFile.writeVector("resMomentumDiff",momentumDiff);
        outFile.writeVector("resMomentumConv",momentumConv);
        outFile.writeVector("resMomentumProd",momentumProd);
        outFile.writeArray2D("resSpeciesUnst",speciesUnst);
        outFile.writeArray2D("resSpeciesDiff",speciesDiff);
        outFile.writeArray2D("resSpeciesConv",speciesConv);
        outFile.writeArray2D("resSpeciesProd",speciesProd);
        outFile.writeVector("resContinuityUnst",continuityUnst);
        outFile.writeVector("resContinuityRhov",continuityRhov);
        outFile.writeVector("resContinuityStrain",continuityStrain);
    }

    outFile.close();
    if (incrementFileNumber) {
        options.outputFileNumber++;
    }

    if (errorFile && options.stopIfError) {
      cout << "Error outputs remaining until termination: " << options.errorStopCount << endl;
      if (options.errorStopCount-- <= 0) {
        throw debugException("Too many integration failures.");
      }
    }
}

double flameSys::strainRate(const double t)
{
    // Strain rate is at initial value until time strainRateT0,
    // then increases linearly to the final value at time strainRateT0+strainRateDt
    return (t <= strainRateT0) ? strainRateInitial
        :  (t >= strainRateT0+strainRateDt) ? strainRateFinal
        : strainRateInitial + (strainRateFinal-strainRateInitial)*(t-strainRateT0)/strainRateDt;
}

double flameSys::dStrainRatedt(const double t)
{
    return (t>tPrev) ? (strainRate(t)-aPrev)/(t-tPrev) : 0;
}

void flameSys::updateLeftBC(void)
{
    // Boundary condition for left edge of domain
    int prev = grid.leftBoundaryConfig;

    if ((options.twinFlame || options.curvedFlame) && x[0] >= 0.0 && x[0] <= options.centerGridMin) {
        grid.leftBoundaryConfig = grid.lbControlVolume;
    } else if (grid.ju == 0 || (grid.jb == 0 && grid.fixedBurnedVal)) {
        grid.leftBoundaryConfig = grid.lbFixedVal;
    } else {
        grid.leftBoundaryConfig = grid.lbZeroGradNonCenter;
    }

    if (prev != grid.leftBoundaryConfig) {
        cout << "updateLeftBC: BC changed from " << prev << " to " << grid.leftBoundaryConfig << "." << endl;
    }
}

double& flameSys::jacA(const int j, const int k1, const int k2)
{
    return (*bandedJacobian)(j*nVars+k1, (j-1)*nVars+k2);
}

double& flameSys::jacB(const int j, const int k1, const int k2)
{
    return (*bandedJacobian)(j*nVars+k1, j*nVars+k2);
}

double& flameSys::jacC(const int j, const int k1, const int k2)
{
    return (*bandedJacobian)(j*nVars+k1, (j+1)*nVars+k2);
}

double flameSys::getHeatReleaseRate(void)
{
    return mathUtils::integrate(x, qDot);
}

double flameSys::getConsumptionSpeed(void)
{
    double QoverCp = mathUtils::integrate(x,qDot/cp);
    double rhouDeltaT = rhou*(Tb-Tu);
    return QoverCp/rhouDeltaT;
}

double flameSys::getFlamePosition(void)
{
    return mathUtils::trapz(x,x*qDot)/mathUtils::trapz(x,qDot);
}

double flameSys::targetFlamePosition(double t)
{
    return (t <= options.xFlameT0) ? options.xFlameInitial
        :  (t >= options.xFlameT0+options.xFlameDt) ? options.xFlameFinal
        : options.xFlameInitial + (options.xFlameFinal-options.xFlameInitial)*(t-options.xFlameT0)/options.xFlameDt;
}

void flameSys::update_xStag(const double t, const bool updateIntError)
{
    xFlameActual = getFlamePosition();
    xFlameTarget = targetFlamePosition(t);
    if (updateIntError) {
        flamePosIntegralError += (xFlameTarget-xFlameActual)*(t-tFlamePrev);
        tFlamePrev = t;
    }

    // controlSignal is approximately a*xStag
    double controlSignal = options.xFlameProportionalGain *
        ( (xFlameTarget-xFlameActual) + (flamePosIntegralError + (xFlameTarget-xFlameActual)*(t-tFlamePrev))*options.xFlameIntegralGain );

    if (debugParameters::debugFlameRadiusControl) {
        cout << "rFlameControl: " << "rF = " << xFlameActual << "   control = " << controlSignal;
        cout << "   P = " <<  options.xFlameProportionalGain*(xFlameTarget-xFlameActual);
        cout << "   I = " << options.xFlameProportionalGain*flamePosIntegralError*options.xFlameIntegralGain;
        cout << "  dt = " << t-tFlamePrev << endl;
    }

    double a = strainRate(t);
    if (alpha == 1) {
        rVzero = 0.5*rhoLeft*(controlSignal*abs(controlSignal)-a*x[0]*x[0]);
    } else {
        rVzero = rhoLeft*(controlSignal-a*x[0]);
    }

}

void flameSys::printPerformanceStats(void)
{
    cout << endl << " **Performance Stats**:  time  (call count)" << endl;
    printPerfString("         General Setup: ", perfTimerSetup);
    printPerfString("  Preconditioner Setup: ", perfTimerPrecondSetup);
    printPerfString("  Factorizing Jacobian: ", perfTimerLU);
    printPerfString("  Preconditioner Solve: ", perfTimerPrecondSolve);
    printPerfString("   Residual Evaluation: ", perfTimerResFunc);
    printPerfString("        Reaction Rates: ", perfTimerRxnRates);
    printPerfString("  Transport Properties: ", perfTimerTransportProps);
    cout << endl;
}

void flameSys::printPerfString(const std::string& label, const perfTimer& T) const
{
    cout << label << mathUtils::stringify(T.getTime(),6) << " (" << T.getCallCount() << ")" << endl;
}

void flameSys::testPreconditioner(void)
{
    // This function computes a handful of Jacobian elements,
    // both using the analytic Jacobian function and numerically
    // by evaluating (f(y+dy)-f(y))/dy. Since the latter method
    // would take too long normally, we only look at 3 points: j=1, j=jj/2 and j=jj
    // This function is for debugging purposes only, and is not called when
    // the code is running normally.

    int jj = nPoints - 1;

    try {
        updateThermoProperties();
        updateTransportProperties();
        gas.getReactionRates(wDot);
    } catch (Cantera::CanteraError) {
        writeStateFile("errorOutput",true);
        throw debugException("Error evaluating thermodynamic properties");
    }

    inTestPreconditioner = true;
    sdVector y(N), ypdy(N), yDot(N), res1(N), res2(N);
    rollY(y);
    rollYdot(yDot);
    f(tNow, y, yDot, res1);

    // First, evaluate using the usual preconditioner function:
    preconditionerSetup(tNow, y, yDot, res1, 0);

    int j2 = jj/2;
    Array2D firstBlock1(nVars,2*nVars), secondBlock1(nVars,3*nVars), thirdBlock1(nVars,2*nVars);
    for (int k=0; k<nVars; k++) {
        for (int i=0; i<nVars; i++) {
            firstBlock1(k,i) = jacB(0,k,i);
            firstBlock1(k,i+nVars) = jacC(0,k,i);

            secondBlock1(k,i) = jacA(j2,k,i);
            secondBlock1(k,i+nVars) = jacB(j2,k,i);
            secondBlock1(k,i+2*nVars) = jacC(j2,k,i);

            thirdBlock1(k,i) = jacA(jj,k,i);
            thirdBlock1(k,i+nVars) = jacB(jj,k,i);
        }
    }

    // Now by finite differences:
    double eps = sqrt(DBL_EPSILON)*800;
    Array2D firstBlock2(nVars,2*nVars), secondBlock2(nVars,3*nVars), thirdBlock2(nVars,2*nVars);

    for (int i=0; i<2*nVars; i++) {
        for (int k=0; k<N; k++) {
            ypdy(k) = y(k);
        }
        ypdy(i) = (abs(y(i)) > eps/2) ? y(i)*(1+eps) : eps;
        f(tNow, ypdy, yDot, res2);
        for (int k=0; k<nVars; k++) {
            firstBlock2(k,i) = (res2(k)-res1(k))/(ypdy(i)-y(i));
        }
    }

    for (int i=0; i<3*nVars; i++) {
        for (int k=0; k<N; k++) {
            ypdy(k) = y(k);
        }
        int iMod = i + (j2-1)*nVars;
        ypdy(iMod) = (abs(y(iMod)) > eps/2) ? y(iMod)*(1+eps) : eps;
        f(tNow, ypdy, yDot, res2);
        for (int k=0; k<nVars; k++) {
            secondBlock2(k,i) = (res2(k+j2*nVars)-res1(k+j2*nVars))/(ypdy(iMod)-y(iMod));
        }
    }

    for (int i=0; i<2*nVars; i++) {
        for (int k=0; k<N; k++) {
            ypdy(k) = y(k);
        }
        int iMod = i + (jj-1)*nVars;
        ypdy(iMod) = (abs(y(iMod)) > eps/2) ? y(iMod)*(1+eps) : eps;
        f(tNow, ypdy, yDot, res2);
        for (int k=0; k<nVars; k++) {
            thirdBlock2(k,i) = (res2(k+jj*nVars)-res1(k+jj*nVars))/(ypdy(iMod)-y(iMod));
        }
    }

    // And then save the results:
    DataFile jacobianFile(options.outputDir+"/jacobianComparison.h5");
    jacobianFile.writeArray2D("aFirst",firstBlock1);
    jacobianFile.writeArray2D("aMiddle",secondBlock1);
    jacobianFile.writeArray2D("aEnd",thirdBlock1);
    jacobianFile.writeArray2D("nFirst",firstBlock2);
    jacobianFile.writeArray2D("nMiddle",secondBlock2);
    jacobianFile.writeArray2D("nEnd",thirdBlock2);
    jacobianFile.close();

    cout << "Wrote output file: " << options.outputDir+"/jacobianComparison.h5" << endl;
    inTestPreconditioner = false;
}

void flameSys::debugFailedTimestep(const sdVector& y)
{
    sdVector res(N);
    rollResiduals(res);
    dvector nRes(N);
    for (int i=0; i<N; i++) {
        nRes[i] = abs(res(i))/(reltol*abs(y(i))+(*abstol)(i));
    }
    cout << "Largest normalized residual components:" << endl;
    for (int i=0; i<10; i++) {
        int m = maxloc(nRes);
        if (nRes[m] > 0.0) {
            int j = m / nVars;
            int k = m % nVars;
            std::string kStr = (k == 0) ? "V" :
                               (k == 1) ? "U" :
                               (k == 2) ? "T"
                                        : "Y_"+stringify(k-3);
            cout << "    j = " << j << ", residual of " << kStr << " = " << nRes[m] << endl;
            nRes[m] = 0;
        }
    }
}
