#include "flameSys.h"
#include <sstream>
#include <vector>
#include <cmath>
#include "debugUtils.h"
#include "mathUtils.h"
#include "matlabFile.h"
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

    if (options.singleCanteraObject) {
        simpleGas.setStateMass(Y,T);
    } else {
        gas.setStateMass(Y,T);
    }

    try {
        updateThermoProperties();
        if (!options.steadyOnly || inGetIC) {
            perfTimerTransportProps.start();
            updateTransportProperties();
            perfTimerTransportProps.stop();
        }
        perfTimerRxnRates.start();
        if (options.singleCanteraObject) {
            simpleGas.getReactionRates(wDot);
        } else {
            gas.getReactionRates(wDot);
        }
        perfTimerRxnRates.stop();
    } catch (Cantera::CanteraError) {
        cout << "Error evaluating thermodynamic properties" << endl;
        writeStateMatFile("errorOutput",true);
        throw;
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

    if (options.flameRadiusControl && !inGetIC && !inTestPreconditioner) {
        update_rStag(t, false);
    }

    // *** Calculate diffusion mass fluxes, heat flux, enthalpy flux
    for (int j=0; j<jj; j++) {
        sumcpj[j] = 0;
        for (int k=0; k<nSpec; k++) {
            jFick(k,j) = -0.5*(rhoD(k,j)+rhoD(k,j+1))
                * (Y(k,j+1)-Y(k,j))/grid.hh[j];
            jSoret(k,j) = -0.5*(Dkt(k,j)/T[j] + Dkt(k,j+1)/T[j+1])
                * (T[j+1]-T[j])/grid.hh[j];

        }
        qFourier[j] = -0.5*(lambda[j]+lambda[j+1])*(T[j+1]-T[j])/grid.hh[j];
        if (inGetIC) {
            jCorr[j] = 0;
            for (int k=0; k<nSpec; k++) {
                jCorr[j] -= jFick(k,j);
            }
        }
        for (int k=0; k<nSpec; k++) {
            jFick(k,j) += 0.5*(Y(k,j)+Y(k,j+1))*jCorr[j]; // correction to ensure that sum of mass fractions equals 1
            sumcpj[j] += cpSpec(k,j)/W[k]*(jFick(k,j) + jSoret(k,j));
        }
    }

    // ****************************************
    // *** Left Boundary values for U, T, Y ***
    // ****************************************

    if (grid.leftBoundaryConfig == grid.lbFixedVal) {
        // Fixed values for T, Y

        energyUnst[0] = dTdt[0];
        energyDiff[0] = energyConv[0] = energyProd[0] = 0;

        // Momentum equation is always zero gradient on the burned side
        if (grid.unburnedLeft) {
            momentumUnst[0] = rho[0]*dUdt[0];
            momentumProd[0] = rho[0]*U[0]*U[0] - rhou*(dadt + a*a);
            momentumDiff[0] = momentumConv[0] = 0;
        } else {
            momentumDiff[0] = (U[1]-U[0])/grid.hh[0];
            momentumProd[0] = momentumConv[0] = momentumProd[0] = 0;
        }

        for (int k=0; k<nSpec; k++) {
            speciesUnst(k,0) = dYdt(k,0);
            speciesDiff(k,0) = speciesConv(k,0) = speciesProd(k,0) = 0;
        }

    } else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
        // Zero gradient condition for T, Y, U, not at centerline
        energyDiff[0] = (T[1]-T[0])/grid.hh[0];
        momentumDiff[0] = (U[1]-U[0])/grid.hh[0];
        energyProd[0] = energyConv[0] = energyProd[0] = 0;
        momentumProd[0] = momentumConv[0] = momentumProd[0] = 0;

        for (int k=0; k<nSpec; k++) {
            speciesDiff(k,0) = (Y(k,1)-Y(k,0))/grid.hh[0];
            speciesUnst(k,0) = speciesConv(k,0) = speciesProd(k,0) = 0;
        }

    } else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
        // Left boundary corresponds to centerline or symmetry plane
        centerVol = pow(grid.x[0],grid.alpha+1)/(grid.alpha+1);
        centerArea = pow(grid.x[0],grid.alpha);
        double rVzero = (rV[0] > 0) ? rV[0] : 0;

        energyUnst[0] = rho[0]*dTdt[0];
        energyDiff[0] = centerArea/centerVol*qFourier[0]/cp[0];
//        energyProd[0] = -qDot[0]/cp[0]*0.0;
        energyProd[0] = 0;
        energyConv[0] = rVzero*(T[0]-Tleft)/centerVol;

        momentumUnst[0] = rho[0]*dUdt[0];
        momentumDiff[0] = -centerArea/centerVol*0.5*(mu[0]+mu[1])*(U[1]-U[0])/grid.hh[0];
//        momentumProd[0] = (dadt*(rho[0]*U[0]-rhou)/a + (rho[0]*U[0]*U[0]-rhou)*a)*0.0;
        momentumProd[0] = 0;
        momentumConv[0] = rVzero*(U[0]-Uleft)/centerVol;

        for (int k=0; k<nSpec; k++) {
            speciesUnst(k,0) = rho[0]*dYdt(k,0);
            speciesDiff(k,0) = centerArea/centerVol*(jFick(k,0) + jSoret(k,0));
//            speciesProd(k,0) = -wDot(k,0)*W[k]*0.0;
            speciesProd(k,0) = 0;
            speciesConv(k,0) = rVzero*(Y(k,0)-Yleft[k])/centerVol;
        }
    }

    // ***************************************
    // *** Intermediate points for U, T, Y ***
    // ***************************************

    for (int j=1; j<jj; j++) {
        if (options.centeredDifferences) {
            // First derivative for convective terms: centered difference
            dTdx[j] = (T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j]);
            dUdx[j] = (U[j-1]*grid.cfm[j] + U[j]*grid.cf[j] + U[j+1]*grid.cfp[j]);
            for (int k=0; k<nSpec; k++) {
                dYdx(k,j) = (Y(k,j-1)*grid.cfm[j] + Y(k,j)*grid.cf[j] + Y(k,j+1)*grid.cfp[j]);
            }
        } else {
            // Upwinded convective derivatives:
            if (V[j] > 0) {
                dUdx[j] = (U[j]-U[j-1])/grid.hh[j-1];
                dTdx[j] = (T[j]-T[j-1])/grid.hh[j-1];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j)-Y(k,j-1))/grid.hh[j-1];
                }
            } else {
                dUdx[j] = (U[j+1]-U[j])/grid.hh[j];
                dTdx[j] = (T[j+1]-T[j])/grid.hh[j];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j+1)-Y(k,j))/grid.hh[j];
                }
            }
        }

        // The enthalpy flux term always uses centered differences
        dTdxCen[j] = T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j];

        // *** Momentum Equation
        momentumUnst[j] = rho[j]*dUdt[j];
        momentumConv[j] = V[j]*dUdx[j];
        momentumDiff[j] = -0.5*( grid.rphalf[j]*(mu[j]+mu[j+1])*(U[j+1]-U[j])/grid.hh[j] -
            grid.rphalf[j-1]*(mu[j-1]+mu[j])*(U[j]-U[j-1])/grid.hh[j-1])/(grid.dlj[j]*grid.r[j]);
        momentumProd[j] = rho[j]*U[j]*U[j] - rhou*(dadt + a*a);

        // *** Energy Equation
        energyUnst[j] = rho[j]*dTdt[j];
        energyConv[j] = V[j]*dTdx[j];
        energyDiff[j] = 0.5*(sumcpj[j]+sumcpj[j-1])*dTdxCen[j]/cp[j] +
            (grid.rphalf[j]*qFourier[j] - grid.rphalf[j-1]*qFourier[j-1])/(grid.dlj[j]*cp[j]*grid.r[j]);
        energyProd[j] = -qDot[j]/cp[j];

        // *** Species Equations
        for (int k=0; k<nSpec; k++) {
            speciesUnst(k,j) = rho[j]*dYdt(k,j);
            speciesConv(k,j) = V[j]*dYdx(k,j);
            speciesDiff(k,j) = (grid.rphalf[j]*(jFick(k,j)+jSoret(k,j))
                - grid.rphalf[j-1]*(jFick(k,j-1)+jSoret(k,j-1)))/(grid.dlj[j]*grid.r[j]);
            speciesProd(k,j) = -wDot(k,j)*W[k];
        }
    }

    // *****************************
    // *** Right boundary values ***
    // *****************************

    // *** Right boundary values for T and Y
    if (grid.unburnedLeft && !grid.fixedBurnedVal) {
        // zero gradient condition
        energyDiff[jj] = (T[jj]-T[jj-1])/grid.hh[jj-1];
        energyProd[jj] = energyConv[jj] = energyUnst[jj] = 0;

        for (int k=0; k<nSpec; k++) {
            speciesDiff(k,jj) = (Y(k,jj)-Y(k,jj-1))/grid.hh[jj-1];
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

    // *** Right boundary values for U
    if (grid.unburnedLeft) {
        momentumDiff[jj] = (U[jj]-U[jj-1])/grid.hh[jj-1];
        momentumProd[jj] = momentumConv[jj] = momentumUnst[jj] = 0;
    } else {
        momentumUnst[jj] = rho[jj]*dUdt[jj];
        momentumProd[jj] = rho[jj]*U[jj]*U[jj] - rhou*(dadt + a*a);
        momentumDiff[jj] = momentumConv[jj] = 0;
    }

    // ***************************
    // *** Continuity Equation ***
    // ***************************

    // *** Left boundary value
    if (options.stagnationRadiusControl) {
        // Boundary value for V depends on rStag
        if (grid.alpha == 1) {
          continuityUnst[0] = rV[0] - 0.5*rhoLeft*a*(options.rStag*abs(options.rStag)-grid.x[0]*grid.x[0]);
      } else {
          continuityUnst[0] = rV[0] - rhoLeft*a*(options.rStag-grid.x[0]);
      }

    } else {
        // Boundary value for V is fixed
        continuityUnst[0] = dVdt[0];
    }

    continuityRhov[0] = continuityStrain[0] = 0;

    // *** Intermediate points
    for (int j=1; j<=jj; j++) {
        continuityRhov[j] = (rV[j]-rV[j-1])/(grid.hh[j-1]*grid.rphalf[j-1]);
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

    BandZero(bandedJacobian->forSundials());

    double a = strainRate(t);
    double dadt = dStrainRatedt(t);

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

    if (options.singleCanteraObject) {
        simpleGas.setStateMass(Y,TplusdT);
        simpleGas.getReactionRates(wDot2);
    } else {
        gas.setStateMass(Y,TplusdT);
        gas.getReactionRates(wDot2);
    }

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

        if (options.singleCanteraObject) {
            simpleGas.setStateMass(YplusdY,T);
            simpleGas.getReactionRates(wDot2);
        } else {
            gas.setStateMass(YplusdY,T);
            gas.getReactionRates(wDot2);
        }

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

        // Momentum equation is always zero gradient on the burned side
        if (grid.unburnedLeft) {
            jacB(j,kMomentum,kMomentum) = c_j + 2*rho[0]*U[0];
            jacB(j,kMomentum,kEnergy) = -rho[0]/T[0]*(U[0]*U[0] + dUdt[0]);
            for (int k=0; k<nSpec; k++) {
                jacB(j,kMomentum,kSpecies+k) = rho[0]*(W[k]-Wmx[0])/(W[k]*(1-Y(k,0)))*(U[0]*U[0] + dUdt[0]);
            }
        } else {
            jacB(j,kMomentum,kMomentum) = -1/grid.hh[0];
            jacC(j,kMomentum,kMomentum) = 1/grid.hh[0];
        }

        for (int k=0; k<nSpec; k++) {
            jacB(j,kSpecies+k,kSpecies+k) = c_j;
        }

    } else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
        // Zero gradient condition for T, Y, U, not at centerline
        jacB(j,kEnergy,kEnergy) = -1/grid.hh[0];
        jacC(j,kEnergy,kEnergy) = 1/grid.hh[0];

        jacB(j,kMomentum,kMomentum) = -1/grid.hh[0];
        jacC(j,kMomentum,kMomentum) = 1/grid.hh[0];

        for (int k=0; k<nSpec; k++) {
            jacB(j,kSpecies+k,kSpecies+k) = -1/grid.hh[0];
            jacC(j,kSpecies+k,kSpecies+k) = 1/grid.hh[0];
        }

    } else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
        // Left boundary corresponds to centerline or symmetry plane
        centerVol = pow(grid.x[0],grid.alpha+1)/(grid.alpha+1);
        centerArea = pow(grid.x[0],grid.alpha);
        double rVzero, drVzero;
        if (rV[0] > 0) {
            rVzero = rV[0];
            drVzero = grid.x[0];
        } else {
            rVzero = 0;
            drVzero = 0;
        }

        // dEnergy/dT_j
        jacB(j,kEnergy,kEnergy) = c_j*rho[j] - rho[j]/T[j]*dTdt[j] +
            + centerArea/centerVol*0.5*(lambda[j]+lambda[j+1])/(grid.hh[j]*cp[j]) + rVzero/centerVol;

        // dEnergy/dT_(j+1)
        jacC(j,kEnergy,kEnergy) = - centerArea/centerVol*0.5*(lambda[j]+lambda[j+1])/(grid.hh[j]*cp[j]);

        // dEnergy/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kEnergy,kSpecies+k) = rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j)))*dTdt[j];
        }

        // dEnergy/dV
        jacB(j,kEnergy,kContinuity) = drVzero*(T[0]-Tleft)/centerVol;

        // dMomentum/dT
        jacB(j,kMomentum,kEnergy) = -rho[j]/T[j]*dUdt[j];

        // dMomentum/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kMomentum,kSpecies+k) = rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j))) * dUdt[j];
        }

        // dMomentum/dU_j
        jacB(j,kMomentum,kMomentum) = rho[j]*c_j
            + centerArea/centerVol*0.5*(mu[j]+mu[j+1])/grid.hh[j] + rVzero/centerVol;

        // dMomentum/dU_j+1
        jacC(j,kMomentum,kMomentum) = -centerArea/centerVol*0.5*(mu[j]+mu[j+1])/grid.hh[j];

        // dMomentum/dV
        jacB(j,kMomentum,kContinuity) = drVzero*(U[0]-Uleft)/centerVol;

        for (int k=0; k<nSpec; k++) {
            // dSpecies/dT
            jacB(j,kSpecies+k,kEnergy) = - rho[j]/T[j]*dYdt(k,j);

            // dSpecies/dY_j
            for (int i=0; i<nSpec; i++) {
                jacB(j,kSpecies+k,kSpecies+i) = rho[j]*(W[i]-Wmx[j])/(W[i]*(1-Y(i,j)))*dYdt(k,j);
            }
            jacB(j,kSpecies+k,kSpecies+k) += rho[j]*c_j
                + centerArea/centerVol*0.5*(rhoD(k,j)+rhoD(k,j+1))/grid.hh[j] + rVzero/centerVol;

            // dSpecies/dY_(j+1)
            jacC(j,kSpecies+k,kSpecies+k) = -centerArea/centerVol*0.5*(rhoD(k,j)+rhoD(k,j+1))/grid.hh[j];

            // dSpecies/dV
            jacB(j,kSpecies+k,kContinuity) = drVzero*(Y(k,0)-Yleft[k])/centerVol;
        }

    }

    // *** Intermediate points for U, T, Y
    for (j=1; j<jj; j++) {

        if (options.centeredDifferences) {
            // First derivative for convective terms: centered difference
            dTdx[j] = (T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j]);
            dUdx[j] = (U[j-1]*grid.cfm[j] + U[j]*grid.cf[j] + U[j+1]*grid.cfp[j]);
            for (int k=0; k<nSpec; k++) {
                dYdx(k,j) = (Y(k,j-1)*grid.cfm[j] + Y(k,j)*grid.cf[j] + Y(k,j+1)*grid.cfp[j]);
            }
        } else {
            // Upwinded convective derivatives:
            if (V[j] > 0) {
                dUdx[j] = (U[j]-U[j-1])/grid.hh[j-1];
                dTdx[j] = (T[j]-T[j-1])/grid.hh[j-1];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j)-Y(k,j-1))/grid.hh[j-1];
                }
            } else {
                dUdx[j] = (U[j+1]-U[j])/grid.hh[j];
                dTdx[j] = (T[j+1]-T[j])/grid.hh[j];
                for (int k=0; k<nSpec; k++) {
                    dYdx(k,j) = (Y(k,j+1)-Y(k,j))/grid.hh[j];
                }
            }
        }

        dTdxCen[j] = T[j-1]*grid.cfm[j] + T[j]*grid.cf[j] + T[j+1]*grid.cfp[j];

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
                0.5*grid.rphalf[j]*((rhoD(k,j)+rhoD(k,j+1))/grid.hh[j]+jCorr[j])/(grid.dlj[j]*grid.r[j]) +
                0.5*grid.rphalf[j-1]*((rhoD(k,j-1)+rhoD(k,j))/grid.hh[j-1]-jCorr[j-1])/(grid.dlj[j]*grid.r[j]);

            // dSpeciesk/dYk_j-1
            jacA(j,kSpecies+k,kSpecies+k) = -0.5*grid.rphalf[j-1]
                * ((rhoD(k,j-1)+rhoD(k,j))/grid.hh[j-1]+jCorr[j-1])
                / (grid.dlj[j]*grid.r[j]);

            // dSpeciesk/dYk_j+1
            jacC(j,kSpecies+k,kSpecies+k) = -0.5*grid.rphalf[j]
                * ((rhoD(k,j)+rhoD(k,j+1))/grid.hh[j]-jCorr[j])
                / (grid.dlj[j]*grid.r[j]);

            if (options.centeredDifferences) {
                jacA(j,kSpecies+k,kSpecies+k) += V[j]*grid.cfm[j];
                jacB(j,kSpecies+k,kSpecies+k) += V[j]*grid.cf[j];
                jacC(j,kSpecies+k,kSpecies+k) += V[j]*grid.cfp[j];
            } else {
                if (V[j] > 0) {
                    jacB(j,kSpecies+k,kSpecies+k) += V[j]/grid.hh[j-1];
                    jacA(j,kSpecies+k,kSpecies+k) -= V[j]/grid.hh[j-1];
                } else {
                    jacB(j,kSpecies+k,kSpecies+k) -= V[j]/grid.hh[j];
                    jacC(j,kSpecies+k,kSpecies+k) += V[j]/grid.hh[j];
                }
            }
        }

        // dEnergy/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kEnergy,kSpecies+k) = hdwdY(k,j)/cp[j] + rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j)))*dTdt[j];

            // Enthalpy flux term
            jacA(j,kEnergy,kSpecies+k) += 0.25*dTdxCen[j]*cpSpec(k,j)/W[k]*(rhoD(k,j)+rhoD(k,j-1))/cp[j]/grid.hh[j-1];
            jacB(j,kEnergy,kSpecies+k) -= 0.25*dTdxCen[j]*cpSpec(k,j)/W[k]*(rhoD(k,j)+rhoD(k,j-1))/cp[j]/grid.hh[j-1];
            jacB(j,kEnergy,kSpecies+k) += 0.25*dTdxCen[j]*cpSpec(k,j)/W[k]*(rhoD(k,j)+rhoD(k,j+1))/cp[j]/grid.hh[j];
            jacC(j,kEnergy,kSpecies+k) -= 0.25*dTdxCen[j]*cpSpec(k,j)/W[k]*(rhoD(k,j)+rhoD(k,j+1))/cp[j]/grid.hh[j];
        }

        // dEnergy/dT_j
        jacB(j,kEnergy,kEnergy) = rho[j]*c_j - rho[j]/T[j]*dTdt[j] +
            (0.5*grid.rphalf[j]*(lambda[j]+lambda[j+1])/grid.hh[j] +
            0.5*grid.rphalf[j-1]*(lambda[j-1]+lambda[j])/grid.hh[j-1])/(grid.dlj[j]*cp[j]*grid.r[j]) +
            dwhdT[j]/cp[j] + 0.5*(sumcpj[j]+sumcpj[j-1])*grid.cf[j]/cp[j];


        // dEnergy/dT_j-1
        jacA(j,kEnergy,kEnergy) = -0.5*grid.rphalf[j-1]*(lambda[j-1]+lambda[j])/(grid.hh[j-1]*grid.dlj[j]*cp[j]*grid.r[j]) +
            0.5*(sumcpj[j]+sumcpj[j+1])*grid.cfm[j]/cp[j];

        // dEnergy/dT_j+1
        jacC(j,kEnergy,kEnergy) = -0.5*grid.rphalf[j]*(lambda[j]+lambda[j+1])/(grid.hh[j]*grid.dlj[j]*cp[j]*grid.r[j]) +
            0.5*(sumcpj[j]+sumcpj[j-1])*grid.cfp[j]/cp[j];

        if (options.centeredDifferences) {
            jacA(j,kEnergy,kEnergy) += V[j]*grid.cfm[j];
            jacB(j,kEnergy,kEnergy) += V[j]*grid.cf[j];
            jacC(j,kEnergy,kEnergy) += V[j]*grid.cfp[j];
        } else {
            if (V[j] > 0) {
                jacB(j,kEnergy,kEnergy) += V[j]/grid.hh[j-1];
                jacA(j,kEnergy,kEnergy) -= V[j]/grid.hh[j-1];
            } else {
                jacB(j,kEnergy,kEnergy) -= V[j]/grid.hh[j];
                jacC(j,kEnergy,kEnergy) += V[j]/grid.hh[j];
            }
        }

        // dEnergy/dV
        jacB(j,kEnergy,kContinuity) = dTdx[j];

        // dMomentum/dY
        for (int k=0; k<nSpec; k++) {
            jacB(j,kMomentum,kSpecies+k) =  rho[j]*(W[k]-Wmx[j])/(W[k]*(1-Y(k,j))) *
                (dUdt[j] + a*U[j]*U[j] + U[j]*dadt/a);
        }

        // dMomentum/dT
        jacB(j,kMomentum,kEnergy) = -rho[j]/T[j] *
            (dUdt[j] + a*U[j]*U[j] + U[j]*dadt/a);

        // dMomentum/dU_j
        jacB(j,kMomentum,kMomentum) = rho[j]*c_j +
            0.5*(grid.rphalf[j]*(mu[j]+mu[j+1])/grid.hh[j] +
            grid.rphalf[j]*(mu[j-1]+mu[j])/grid.hh[j-1])/(grid.dlj[j]*grid.r[j]) +
            2*U[j]*rho[j];

        // dMomentum/dU_j-1
        jacA(j,kMomentum,kMomentum) = -0.5*grid.rphalf[j]*(mu[j-1]+mu[j])/(grid.hh[j-1]*grid.dlj[j]*grid.r[j]);

        // dMomentum/dU_j+1
        jacC(j,kMomentum,kMomentum) = -0.5*grid.rphalf[j]*(mu[j]+mu[j+1])/(grid.hh[j]*grid.dlj[j]*grid.r[j]);

        if (options.centeredDifferences) {
            jacA(j,kMomentum,kMomentum) += V[j]*grid.cfm[j];
            jacB(j,kMomentum,kMomentum) += V[j]*grid.cf[j];
            jacC(j,kMomentum,kMomentum) += V[j]*grid.cfp[j];
        } else {
            if (V[j] > 0) {
                jacB(j,kMomentum,kMomentum) += V[j]/grid.hh[j-1];
                jacA(j,kMomentum,kMomentum) -= V[j]/grid.hh[j-1];
            } else {
                jacB(j,kMomentum,kMomentum) -= V[j]/grid.hh[j];
                jacC(j,kMomentum,kMomentum) += V[j]/grid.hh[j];
            }
        }

        // dMomentum/dV
        jacB(j,kMomentum,kContinuity) = dUdx[j];
    }

    // *** Right boundary values for T, Y
    j = jj;
    if (grid.unburnedLeft && !grid.fixedBurnedVal) {
        // zero gradient condition
        jacA(j,kEnergy,kEnergy) = -1/grid.hh[j-1];
        jacB(j,kEnergy,kEnergy) = 1/grid.hh[j-1];
        for (int k=0; k<nSpec; k++) {
            jacA(j,kSpecies+k,kSpecies+k) = -1/grid.hh[j-1];
            jacB(j,kSpecies+k,kSpecies+k) = 1/grid.hh[j-1];
        }
    } else {
        // Fixed values
        jacB(j,kEnergy,kEnergy) = c_j;
        for (int k=0; k<nSpec; k++) {
            jacB(j,kSpecies+k,kSpecies+k) = c_j;
        }
    }

    // *** Right boundary values for U
    if (grid.unburnedLeft) {
        // zero gradient
        jacA(j,kMomentum,kMomentum) = -1/grid.hh[j-1];
        jacB(j,kMomentum,kMomentum) = 1/grid.hh[j-1];
    } else {
        // fixed value
        jacB(j,kMomentum,kMomentum) = rho[jj]*c_j + 2*U[jj]*rho[jj];
        jacB(j,kMomentum,kEnergy) = -rho[jj]/T[jj]*(U[jj]*U[jj] + dUdt[jj]);
        for (int k=0; k<nSpec; k++) {
            jacB(j,kMomentum,kSpecies+k) = rho[jj]*(W[k]-Wmx[jj])/(W[k]*(1-Y(k,jj)))*(U[jj]*U[jj] + dUdt[jj]);
        }
    }

    // *** Continuity Equation
    if (options.stagnationRadiusControl) {
        jacB(0,kContinuity,kContinuity) = pow(grid.x[0],grid.alpha);
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
        jacB(j,kContinuity,kContinuity) = grid.r[j]/grid.hh[j-1]/grid.rphalf[j-1];

        // dContinuity/dV_(j-1)
        jacA(j,kContinuity,kContinuity) = - grid.r[j-1]/grid.hh[j-1]/grid.rphalf[j-1];
    }

    perfTimerPrecondSetup.stop();

    // *** Get LU Factorization of the Jacobian
    if (!inTestPreconditioner) {
        perfTimerLU.start();
        long int iError = BandGBTRF(bandedJacobian->forSundials(),&pMat[0]);
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

    if (options.singleCanteraObject) {
        simpleGas.resize(nPoints);
        nSpec = simpleGas.nSpec;
    } else {
        gas.resize(nPoints);
        nSpec = gas.nSpec;
    }

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

    X.resize(nSpec,nPoints);
    dXdx.resize(nSpec,nPoints);
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
    BandZero(bandedJacobian->forSundials());

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

    if (options.singleCanteraObject) {
        simpleGas.mechanismFile = options.gasMechanismFile;
        simpleGas.phaseID = options.gasPhaseID;
        simpleGas.pressure = options.pressure;
    } else {
        gas.mechanismFile = options.gasMechanismFile;
        gas.phaseID = options.gasPhaseID;
        gas.pressure = options.pressure;
    }
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

    grid.alpha = options.gridAlpha;
    grid.kMomentum = options.kMomentum;
    grid.kContinuity = options.kContinuity;
    grid.kEnergy = options.kEnergy;
    grid.kSpecies = options.kSpecies;

    kMomentum = options.kMomentum;
    kContinuity = options.kContinuity;
    kEnergy = options.kEnergy;
    kSpecies = options.kSpecies;
}

void flameSys::generateInitialProfiles(void)
{
    cout << "Generating initial profiles from given fuel and oxidizer compositions." << endl;
    grid.x.resize(nPoints);
    int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.
    int jl = jm - 4;
    int jr = jm + 4;
    Yb.resize(nSpec); Yu.resize(nSpec);
    double a = strainRate(tStart);

    Tu = options.Tu;
    if (options.singleCanteraObject) {
        // Reactants
        simpleGas.thermo.setState_TPX(Tu,simpleGas.pressure,&options.reactants[0]);
        rhou = simpleGas.thermo.density();
        simpleGas.thermo.getMassFractions(&Yu[0]);
        simpleGas.thermo.getMassFractions(&Y(0,grid.ju));

        // Products
        simpleGas.thermo.setState_TPX(Tu,simpleGas.pressure,&options.reactants[0]);
        Cantera::equilibrate(simpleGas.thermo,"HP");
        Tb = simpleGas.thermo.temperature();
        rhob = simpleGas.thermo.density();
        simpleGas.thermo.getMassFractions(&Yb[0]);
        simpleGas.thermo.getMassFractions(&Y(0,grid.jb));

        // Diluent in the middle
        simpleGas.thermo.setState_TPY(Tu,simpleGas.pressure,options.oxidizer);
        simpleGas.thermo.getMassFractions(&Y(0,jm));

    } else {
        // Reactants
        gas[grid.ju].setState_TPX(Tu,gas.pressure,&options.reactants[0]);
        rhou = gas[grid.ju].density();
        gas[grid.ju].getMassFractions(&Yu[0]);
        gas[grid.ju].getMassFractions(&Y(0,grid.ju));

        // Products
        gas[grid.jb].setState_TPX(Tu,gas.pressure,&options.reactants[0]);
        Cantera::equilibrate(gas[grid.jb],"HP");
        Tb = gas[grid.jb].temperature();
        rhob = gas[grid.jb].density();
        gas[grid.jb].getMassFractions(&Yb[0]);
        gas[grid.jb].getMassFractions(&Y(0,grid.jb));

        // Diluent in the center to delay ignition
        gas[jm].setState_TPY(Tu,gas.pressure,options.oxidizer);
        gas[jm].getMassFractions(&Y(0,jm));

    }
    Ub = a*sqrt(rhou/rhob);

    if (options.unburnedLeft) {
        rhoLeft = rhou;
        Tleft = Tu;
        Yleft = Yu;
        Uleft = a;
        rhoRight = rhob;
        Tright = Tb;
        Yright = Yb;
        Uright = Ub;
    } else {
        rhoLeft = rhob;
        Tleft = Tb;
        Yleft = Yb;
        Uleft = Ub;
        rhoRight = rhou;
        Tright = Tu;
        Yright = Yu;
        Uright = a;
    }

    T[0] = Tleft; T[grid.jj] = Tright;
    T[jm] = T[grid.ju];

    xLeft = (options.twinFlame || options.curvedFlame) ? max(options.xLeft, options.centerGridMin)
                                                       : options.xLeft;
    xRight = options.xRight;

    // Uniform initial grid
    for (int j=0; j<nPoints; j++) {
        grid.x[j] = xLeft + (xRight-xLeft)*((double) j)/((double) nPoints);
    }

    for (int j=1; j<jl; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0);
        }
        T[j] = T[0];
    }

    for (int j=jl; j<jm; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0) + (Y(k,jm)-Y(k,0))*(grid.x[j]-grid.x[jl])/(grid.x[jm]-grid.x[jl]);
        }
        T[j] = T[0] + (T[jm]-T[0])*(grid.x[j]-grid.x[jl])/(grid.x[jm]-grid.x[jl]);
    }

    for (int j=jm+1; j<jr; j++) {
        for (int k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,jm) + (Y(k,grid.jj)-Y(k,jm))*(grid.x[j]-grid.x[jm])/(grid.x[jr]-grid.x[jm]);
        }
        T[j] = T[jm] + (T[grid.jj]-T[jm])*(grid.x[j]-grid.x[jm])/(grid.x[jr]-grid.x[jm]);
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
    if (options.singleCanteraObject) {
        for (int j=0; j<nPoints; j++) {
            simpleGas.thermo.setState_TPY(T[j],simpleGas.pressure,&Y(0,j));
            rho[j] = simpleGas.thermo.density();
            U[j] = a*sqrt(rhou/rho[j]);
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            gas[j].setState_TPY(T[j],gas.pressure,&Y(0,j));
            rho[j] = gas[j].density();
            U[j] = a*sqrt(rhou/rho[j]);
        }
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
        V[j] = V[j-1] - rho[j]*U[j]*(grid.x[j]-grid.x[j-1]);
    }

    for (int j=jm-1; j>=0; j--) {
        V[j] = V[j+1] + rho[j]*U[j]*(grid.x[j+1]-grid.x[j]);
    }

    if (options.singleCanteraObject) {
        simpleGas.setStateMass(Y,T);
        simpleGas.getMoleFractions(X);
        simpleGas.getMolecularWeights(W);
    } else {
        gas.setStateMass(Y,T);
        gas.getMoleFractions(X);
        gas.getMolecularWeights(W);
    }

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
    matlabFile infile(inputFilename);
    grid.x = infile.readVector("x");

    nPoints = grid.x.size();
    setup();

    U = infile.readVector("U");
    V = infile.readVector("V");
    T = infile.readVector("T");
    Y = infile.readArray2D("Y");
    tStart = infile.readScalar("t");
    double a = strainRate(tStart);
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
        if (options.singleCanteraObject) {
            simpleGas.thermo.setState_TPX(Tu,simpleGas.pressure,&options.reactants[0]);
            simpleGas.thermo.getMassFractions(&Y(0,grid.ju));
            simpleGas.thermo.setState_TPX(Tu,simpleGas.pressure,&options.reactants[0]);
            Cantera::equilibrate(simpleGas.thermo,"HP");
            simpleGas.thermo.getMassFractions(&Y(0,grid.jb));
            T[grid.jb] = simpleGas.thermo.temperature();
        } else {
            gas.thermo(grid.ju).setState_TPX(Tu,gas.pressure,&options.reactants[0]);
            gas.thermo(grid.ju).getMassFractions(&Y(0,grid.ju));
            gas.thermo(grid.jb).setState_TPX(Tu,gas.pressure,&options.reactants[0]);
            Cantera::equilibrate(gas.thermo(grid.jb),"HP");
            gas.thermo(grid.jb).getMassFractions(&Y(0,grid.jb));
            T[grid.jb] = gas.thermo(grid.jb).temperature();
        }
    }
    Tb = T[grid.jb];

    if (options.singleCanteraObject) {
        simpleGas.setStateMass(Y,T);
        simpleGas.getMolecularWeights(W);
        simpleGas.getDensity(rho);
        rhou = rho[grid.ju];
        rhob = rho[grid.jb];
        Yu.resize(nSpec); Yb.resize(nSpec);
        for (int k=0; k<nSpec; k++) {
            Yu[k] = Y(k,grid.ju);
            Yb[k] = Y(k,grid.jb);
        }
    } else {
        gas.setStateMass(Y,T);
        gas.getMolecularWeights(W);
        rhou = gas.thermo(grid.ju).density();
        rhob = gas.thermo(grid.jb).density();
        Yu.resize(nSpec); Yb.resize(nSpec);
        gas[grid.ju].getMassFractions(&Yu[0]);
        gas[grid.jb].getMassFractions(&Yb[0]);
    }
    Ub = a*sqrt(rhou/rhob);

    if (options.unburnedLeft) {
        rhoLeft = rhou;
        Tleft = Tu;
        Yleft = Yu;
        Uleft = a;
        rhoRight = rhob;
        Tright = Tb;
        Yright = Yb;
        Uright = Ub;

    } else {
        rhoLeft = rhob;
        Tleft = Tb;
        Yleft = Yb;
        Uleft = Ub;
        rhoRight = rhou;
        Tright = Tu;
        Yright = Yu;
        Uright = a;
    }

    if (!options.fixedBurnedVal) {
        T[grid.jb] = TbSave;
        for (int k=0; k<nSpec; k++) {
            Y(k,grid.jb) = YbSave[k];
        }
    }

    V2rV();
    updateLeftBC();

    if (grid.leftBoundaryConfig == grid.lbControlVolume && options.flameRadiusControl) {
        if (grid.alpha == 0) {
            options.rStag = V[0]/(rhoLeft*options.strainRateInitial);
        } else {
            double tmp = pow(grid.x[0],2) + 2*rV[0]/(rhoLeft*options.strainRateInitial);
            options.rStag = sign(tmp)*sqrt(abs(tmp));
        }
        flamePosIntegralError = options.rStag/(options.rFlameProportionalGain*options.rFlameIntegralGain);
    }
}

flameSys::flameSys(void)
    : grid(options)
    , bandedJacobian(NULL)
    , flamePosIntegralError(0)
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
        V[j] = v[grid.kContinuity][j];
        U[j] = v[grid.kMomentum][j];
        T[j] = v[grid.kEnergy][j];
    }

    for (int k=0; k<nSpec; k++) {
        for (int j=0; j<nPoints; j++) {
            Y(k,j) = v[grid.kSpecies+k][j];
        }
    }
}

void flameSys::unrollVectorVectorDot(const vector<dvector>& v)
{
    for (int j=0; j<nPoints; j++) {
        dVdt[j] = v[grid.kContinuity][j];
        dUdt[j] = v[grid.kMomentum][j];
        dTdt[j] = v[grid.kEnergy][j];
    }

    for (int k=0; k<nSpec; k++) {
        for (int j=0; j<nPoints; j++) {
            dYdt(k,j) = v[grid.kSpecies+k][j];
        }
    }
}

void flameSys::updateTransportProperties(void)
{
    if (options.singleCanteraObject) {
        simpleGas.getTransportProperties(mu, lambda, rhoD, Dkt);
    } else {
        gas.getViscosity(mu);
        gas.getThermalConductivity(lambda);
        gas.getWeightedDiffusionCoefficients(rhoD);
        gas.getThermalDiffusionCoefficients(Dkt);
    }
}

void flameSys::updateThermoProperties(void)
{
    if (options.singleCanteraObject) {
        simpleGas.getThermoProperties(rho, Wmx, cp, cpSpec, hk);
    } else {
        gas.getDensity(rho);
        gas.getSpecificHeatCapacity(cp);
        gas.getSpecificHeatCapacities(cpSpec);
        gas.getEnthalpies(hk);
        gas.getMixtureMolecularWeight(Wmx);
    }
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

int flameSys::getInitialCondition(double t, sdVector& y, sdVector& ydot, std::vector<bool>& algebraic)
{
    if (debugParameters::debugCalcIC) {
        std::ostringstream fileName(ostringstream::out);
        // Determine the name of the output file (ICaxxxxxx.mat)

        fileName << options.outputDir << "/ICa";
        fileName.flags(ios_base::right);
        fileName.fill('0');
        fileName.width(6);
        fileName << ICfileNumber << ".mat";

        cout << "Writing IC file: " << fileName.str() << endl;
        //cout << "DBG: grid.leftBoundaryConfig" << grid.leftBoundaryConfig << endl;

        // Erase the existing file and create a new one
        if (boost::filesystem::exists(fileName.str())) {
            boost::filesystem::remove(fileName.str());
        }

        dvector resVec(N);
        dvector yVec(N);
        dvector ydotVec(N);
        for (int i=0; i<N; i++) {
            yVec[i] = y(i);
            ydotVec[i] = ydot(i);
        }

        matlabFile outfile(fileName.str());

        outfile.writeScalar("rhoLeft",rhoLeft);
        outfile.writeScalar("rStag",options.rStag);
        outfile.writeScalar("rFlameActual",rFlameActual);
        outfile.writeScalar("rFlameTarget",rFlameTarget);
        outfile.writeScalar("flamePosIntegralError",flamePosIntegralError);
        outfile.writeScalar("tFlamePrev",tFlamePrev);
        outfile.writeScalar("t",t);
        outfile.writeScalar("Pgain",options.rFlameProportionalGain);
        outfile.writeScalar("Igain",options.rFlameIntegralGain);

        outfile.writeVector("x",grid.x);
        outfile.writeVector("y",yVec);
        outfile.writeVector("ydot",ydotVec);
        outfile.writeVector("eUnst",energyUnst);
        outfile.writeVector("eDiff",energyDiff);
        outfile.writeVector("eConv",energyConv);
        outfile.writeVector("eProd",energyProd);
        outfile.writeVector("mUnst",momentumUnst);
        outfile.writeVector("mDiff",momentumDiff);
        outfile.writeVector("mConv",momentumConv);
        outfile.writeVector("mProd",momentumProd);
        outfile.writeVector("cUnst",continuityUnst);
        outfile.writeVector("cStrain",continuityStrain);
        outfile.writeVector("cRhov",continuityRhov);
        outfile.writeArray2D("yUnst",speciesUnst);
        outfile.writeArray2D("yDiff",speciesDiff);
        outfile.writeArray2D("yConv",speciesConv);
        outfile.writeArray2D("yProd",speciesProd);
        outfile.writeArray2D("Y",Y);
        outfile.writeVector("T",T);
        outfile.writeVector("U",U);
        outfile.writeVector("V",V);
        outfile.writeVector("rV",rV);
        outfile.writeArray2D("dYdt",dYdt);
        outfile.writeVector("dTdt",dTdt);
        outfile.writeVector("dUdt",dUdt);

        outfile.writeVector("qDot",qDot);
        outfile.writeArray2D("wDot",wDot);
        outfile.writeArray2D("hk",hk);
        outfile.close();
    }

    inGetIC = true;
    sdVector res(N);

    // Continuity equation: Left boundary value
    double a = strainRate(t);
    double dadt = dStrainRatedt(t);

    if (options.stagnationRadiusControl) {
      if (grid.alpha == 1) {
          rV[0] = 0.5*rhoLeft*a*(options.rStag*abs(options.rStag)-grid.x[0]*grid.x[0]);
      } else {
          rV[0] = rhoLeft*a*(options.rStag-grid.x[0]);
      }
    } else {
        dVdt[0] = 0;
    }
    rV2V();
    rollY(y);
    f(t,y,ydot,res);

    int jj = nPoints-1;

    // Left boundary values
    if (grid.leftBoundaryConfig == grid.lbFixedVal) {
        // Fixed values for T, Y
        dTdt[0] = 0;
        for (int k=0; k<nSpec; k++) {
            dYdt(k,0) = 0;
        }

        if (grid.unburnedLeft) {
            dUdt[0] = -U[0]*U[0] + rhou*(a*a + dadt)/rho[0]; // "equilibrium" value on reactants side
        } else {
            U[0] = U[1]; // zero gradient on products side
        }

    } else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
        T[0] = T[1];
        U[0] = U[1];
        for (int k=0; k<nSpec; k++) {
            Y(k,0) = Y(k,1);
        }
    } else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
        dTdt[0] = -(energyDiff[0]+energyProd[0]+energyConv[0])/rho[0];
        dUdt[0] = -(momentumDiff[0]+momentumProd[0]+momentumConv[0])/rho[0];
        for (int k=0; k<nSpec; k++) {
            dYdt(k,0) = -(speciesDiff(k,0)+speciesProd(k,0)+speciesConv(k,0))/rho[0];
        }
    }

    // Right boundary values
    if (grid.unburnedLeft && !grid.fixedBurnedVal) {
        T[jj] = T[jj-1];
        for (int k=0; k<nSpec; k++) {
            Y(k,jj) = Y(k,jj-1);
        }
    } else {
        dTdt[jj] = 0;
        for (int k=0; k<nSpec; k++) {
            dYdt(k,jj) = 0;
        }
    }

    if (grid.unburnedLeft) {
        U[jj] = U[jj-1];
    } else {
        dUdt[jj] = - U[jj]*U[jj] + rhou*(a*a + dadt)/rho[jj];
    }

    dTdx[jj] = (T[jj]-T[jj-1])/grid.hh[jj-1];
    for (int k=0; k<nSpec; k++) {
        dYdx(k,jj) = (Y(k,jj)-Y(k,jj-1))/grid.hh[jj-1];
    }

    for (int j=1; j<jj; j++) {
        double S = -(energyDiff[j]+energyProd[j])/T[j];
        double B = dTdx[j]/T[j];
        for (int k=0; k<nSpec; k++) {
            S -= (speciesDiff(k,j)+speciesProd(k,j))*Wmx[j]/W[k];
            B += dYdx(k,j)*Wmx[j]/W[k];
        }
        rV[j] = (S - continuityStrain[j] + rV[j-1]/(grid.hh[j-1]*grid.rphalf[j-1]))
            / (1/(grid.hh[j-1]*grid.rphalf[j-1]) + B/grid.r[j]);
    }

    rV2V();
    rollY(y);
    rollYdot(ydot);

    f(t,y,ydot,res);

    // Intermediate points
    for (int j=1; j<jj; j++) {
        dTdt[j] = -(energyDiff[j]+energyProd[j]+energyConv[j])/rho[j];
        dUdt[j] = -(momentumDiff[j]+momentumProd[j]+momentumConv[j])/rho[j];
        for (int k=0; k<nSpec; k++) {
            dYdt(k,j) = -(speciesDiff(k,j)+speciesProd(k,j)+speciesConv(k,j))/rho[j];
        }
    }

    // rV[jj]
    if (grid.unburnedLeft && !grid.fixedBurnedVal) {
        // Zero gradient condition for T,Y
        double sum = 0;
        for (int k=0; k<nSpec; k++) {
            sum += dYdt(k,jj-1)/W[k];
        }
        drhodt[jj-1] = -rho[jj-1]*(dTdt[jj-1]/T[jj-1] + Wmx[jj-1]*sum);
        drhodt[jj] = drhodt[jj-1];
    } else {
        // Fixed values for T,Y
        drhodt[jj] = 0;
    }
    rV[jj] = rV[jj-1] - (rho[jj]*U[jj] + drhodt[jj])*grid.hh[jj-1]*grid.rphalf[jj-1];

    rV2V();
    rollY(y);
    rollYdot(ydot);

    sdVector resTemp(N);
    f(t,y,ydot,resTemp);
    inGetIC = false;

    double resnorm = 0;
    for (int i=0; i<N; i++) {
        resnorm += resTemp(i)*resTemp(i);
    }
    resnorm = sqrt(resnorm);

    cout << "Residual norm after IC calculation: " << resnorm << endl;

    if (debugParameters::debugCalcIC) {
        std::ostringstream fileName(ostringstream::out);
        // Determine the name of the output file (ICbxxxxxx.mat)

        fileName << options.outputDir << "/ICb";
        fileName.flags(ios_base::right);
        fileName.fill('0');
        fileName.width(6);
        fileName << ICfileNumber++ << ".mat";

        cout << "Writing IC file: " << fileName.str() << endl;

        // Erase the existing file and create a new one
        if (boost::filesystem::exists(fileName.str())) {
            boost::filesystem::remove(fileName.str());
        }

        dvector resVec(N);
        dvector yVec(N);
        dvector ydotVec(N);
        for (int i=0; i<N; i++) {
            resVec[i] = resTemp(i);
            yVec[i] = y(i);
            ydotVec[i] = ydot(i);
        }

        matlabFile outfile(fileName.str());

        //writeStateMatFile("errorOutput",true);
        outfile.writeScalar("rhoLeft",rhoLeft);
        outfile.writeScalar("a",a);
        outfile.writeScalar("rStag",options.rStag);
        outfile.writeScalar("rFlameActual",rFlameActual);
        outfile.writeScalar("rFlameTarget",rFlameTarget);
        outfile.writeScalar("flamePosIntegralError",flamePosIntegralError);
        outfile.writeScalar("tFlamePrev",tFlamePrev);
        outfile.writeScalar("t",t);
        outfile.writeScalar("Pgain",options.rFlameProportionalGain);
        outfile.writeScalar("Igain",options.rFlameIntegralGain);

        outfile.writeVector("x",grid.x);
        outfile.writeVector("y",yVec);
        outfile.writeVector("ydot",ydotVec);
        outfile.writeVector("res",resVec);
        outfile.writeVector("eUnst",energyUnst);
        outfile.writeVector("eDiff",energyDiff);
        outfile.writeVector("eConv",energyConv);
        outfile.writeVector("eProd",energyProd);
        outfile.writeVector("mUnst",momentumUnst);
        outfile.writeVector("mDiff",momentumDiff);
        outfile.writeVector("mConv",momentumConv);
        outfile.writeVector("mProd",momentumProd);
        outfile.writeVector("cUnst",continuityUnst);
        outfile.writeVector("cStrain",continuityStrain);
        outfile.writeVector("cRhov",continuityRhov);
        outfile.writeArray2D("yUnst",speciesUnst);
        outfile.writeArray2D("yDiff",speciesDiff);
        outfile.writeArray2D("yConv",speciesConv);
        outfile.writeArray2D("yProd",speciesProd);
        outfile.writeArray2D("Y",Y);
        outfile.writeVector("T",T);
        outfile.writeVector("U",U);
        outfile.writeVector("V",V);
        outfile.writeVector("rV",rV);
        outfile.writeArray2D("dYdt",dYdt);
        outfile.writeVector("dTdt",dTdt);
        cout << "DBG: dTdt[0] = " << dTdt[0] << endl;
        outfile.writeVector("dUdt",dUdt);
        outfile.writeVector("qDot",qDot);
        outfile.writeArray2D("wDot",wDot);
        outfile.writeArray2D("hk",hk);

        outfile.close();
    }

    if (resnorm > 1.0e0) {
        // cout << "IC calculation failed!." << endl;
        if (resnorm > 1.0e3) {
            return 100;
        } else {
            return 1;
        }
    }

    return 0;
}

void flameSys::V2rV(void)
{
    if (grid.alpha == 0) {
        for (int j=0; j<nPoints; j++) {
            rV[j] = V[j];
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            rV[j] = grid.x[j]*V[j];
        }
    }

    // In the case where the centerline is part of the domain,
    // V[0] actually stores rV[0] (since V = Inf there)
    if (grid.x[0] == 0) {
        rV[0] = V[0];
    }
}

void flameSys::rV2V(void)
{
    if (grid.alpha == 0) {
        for (int j=0; j<nPoints; j++) {
            V[j] = rV[j];
        }
    } else {
        for (int j=0; j<nPoints; j++) {
            V[j] = rV[j]/grid.x[j];
        }
    }

    // In the case where the centerline is part of the domain,
    // V[0] actually stores rV[0] (since V = Inf there)
    if (grid.x[0] == 0) {
        V[0] = rV[0];
    }
}

void flameSys::writeStateMatFile(const std::string fileNameStr, bool errorFile)
{
    std::ostringstream fileName(ostringstream::out);
    bool incrementFileNumber = false;

    if (fileNameStr.length() == 0) {
        // Determine the name of the output file (outXXXXXX.mat)
        incrementFileNumber = true;
        if (errorFile) {
            fileName << options.outputDir << "/error";
        } else {
            fileName << options.outputDir << "/prof";
        }
        fileName.flags(ios_base::right);
        fileName.fill('0');
        fileName.width(6);
        fileName << options.outputFileNumber << ".mat";
    } else {
        fileName << options.outputDir << "/" << fileNameStr << ".mat";
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
    matlabFile outFile(fileName.str());

    // Write the state data to the output file:
    outFile.writeScalar("t", tNow);
    outFile.writeVector("x", grid.x);
    outFile.writeVector("T", T);
    outFile.writeVector("U", U);
    outFile.writeVector("V", V);
    outFile.writeArray2D("Y", Y);
    outFile.writeScalar("a",strainRate(tNow));
    outFile.writeScalar("dadt",dStrainRatedt(tNow));
    outFile.writeScalar("fileNumber", options.outputFileNumber);

    if (options.outputHeatReleaseRate || errorFile) {
        outFile.writeVector("q",qDot);
    }

    if (options.outputTimeDerivatives || errorFile) {
        outFile.writeVector("dUdt", dUdt);
        outFile.writeVector("dTdt", dTdt);
        outFile.writeVector("dVdt", dVdt);
        outFile.writeArray2D("dYdt", dYdt);
        outFile.writeVector("drhodt",drhodt);
    }

    if (options.outputAuxiliaryVariables || errorFile) {
        outFile.writeVector("rho", rho);
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
        outFile.writeVector("csp",grid.cfp);
        outFile.writeVector("cs",grid.cf);
        outFile.writeVector("csm",grid.cfm);
        outFile.writeVector("hh",grid.hh);
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
        throw;
      }
    }
}

double flameSys::strainRate(const double t)
{
    // Strain rate is at initial value until time strainRateT0,
    // then increases linearly to the final value at time strainRateT0+strainRateDt
    return (t <= strainRateT0) ? strainRateInitial
        :  (t >= strainRateT0+strainRateDt) ? strainRateFinal
        : strainRateInitial + (strainRateFinal-strainRateInitial)*(t-tStart)/strainRateDt;
}

double flameSys::dStrainRatedt(const double t)
{
    return (t <= strainRateT0) ? 0
        :  (t >= strainRateT0+strainRateDt) ? 0
        : (strainRateFinal-strainRateInitial)/strainRateDt;
}

void flameSys::updateAlgebraicComponents(void)
{
    int jj = nPoints-1;
    algebraic.resize(N);

    algebraic[0] = options.stagnationRadiusControl; // continuity

    if (grid.leftBoundaryConfig == grid.lbFixedVal) {
        algebraic[1] = !grid.unburnedLeft; // momentum
        algebraic[2] = false; // energy
        for (int k=0; k<nSpec; k++) {
            algebraic[3+k] = false; // species
        }
    } else if (grid.leftBoundaryConfig == grid.lbZeroGradNonCenter) {
        algebraic[1] = true; // momentum
        algebraic[2] = true; // energy
        for (int k=0; k<nSpec; k++) {
            algebraic[3+k] = true; // species
        }
    } else if (grid.leftBoundaryConfig == grid.lbControlVolume) {
        algebraic[1] = false; // momentum
        algebraic[2] = false; // energy
        for (int k=0; k<nSpec; k++) {
            algebraic[3+k] = false; // species
        }
    }

    for (int j=1; j<jj; j++) {
        algebraic[nVars*j] = true; // continuity
        algebraic[nVars*j+1] = false; // momentum
        algebraic[nVars*j+2] = false; // energy
        for (int k=0; k<nSpec; k++) {
            algebraic[nVars*j+3+k] = false; //species
        }
    }

    algebraic[nVars*jj] = true; // continuity
    if (grid.unburnedLeft && !grid.fixedBurnedVal) {
        algebraic[nVars*jj+1] = true; // momentum;
        algebraic[nVars*jj+2] = true; // energy
        for (int k=0; k<nSpec; k++) {
            algebraic[nVars*jj+3+k] = true; // species
        }
    } else {
        algebraic[nVars*jj+1] = false; // momentum;
        algebraic[nVars*jj+2] = false; // energy
        for (int k=0; k<nSpec; k++) {
            algebraic[nVars*jj+3+k] = false; // species
        }
    }

    algebraic[nVars*jj+1] = grid.unburnedLeft;

}

void flameSys::updateConstraints(sdVector& constraints)
{
    for (int j=0; j<nPoints; j++) {
        constraints(nVars*j+kContinuity) = 0; // V: no constraint
        constraints(nVars*j+kMomentum) = 0; // U: no constraint
        constraints(nVars*j+kEnergy) = 1.0; // T >= 0
        for (int k=0; k<nSpec; k++) {
            constraints(nVars*j+kSpecies+k) = 1.0; // Y >= 0
        }
    }
}

void flameSys::updateLeftBC(void)
{
    // Boundary condition for left edge of domain
    int prev = grid.leftBoundaryConfig;

    if ((options.twinFlame || options.curvedFlame) && grid.x[0] >= 0.0 && grid.x[0] <= options.centerGridMin) {
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
    if (options.twinFlame || options.curvedFlame) {
        dvector qDotFull(nPoints+1), xFull(nPoints+1);
        qDotFull[0] = qDot[0]; xFull[0] = 0;
        for (int j=0; j<nPoints; j++) {
            qDotFull[j+1] = qDot[j];
            xFull[j+1] = grid.x[j];
        }
        return mathUtils::trapz(xFull,qDotFull);
    } else {
        return mathUtils::integrate(grid.x, qDot);
    }
}

double flameSys::getConsumptionSpeed(void)
{
    double QoverCp;
    if (options.twinFlame || options.curvedFlame) {
        dvector qcp(nPoints+1), xFull(nPoints+1);
        qcp[0] = qDot[0]/cp[0]; xFull[0] = 0;
        for (int j=0; j<nPoints; j++) {
            qcp[j+1] = qDot[j]/cp[j];
            xFull[j+1] = grid.x[j];
        }
        QoverCp = mathUtils::integrate(xFull,qcp);
    } else {
        QoverCp = mathUtils::integrate(grid.x,qDot/cp);
    }
    double rhouDeltaT = rho[grid.ju]*(T[grid.jb]-T[grid.ju]);
    return QoverCp/rhouDeltaT;
}

double flameSys::getFlamePosition(void)
{
    if (options.twinFlame || options.curvedFlame) {
        dvector qDotFull(nPoints+1), xFull(nPoints+1);
        qDotFull[0] = qDot[0]; xFull[0] = 0;
        for (int j=0; j<nPoints; j++) {
            qDotFull[j+1] = qDot[j];
            xFull[j+1] = grid.x[j];
        }
        return mathUtils::trapz(xFull,xFull*qDotFull)/mathUtils::trapz(xFull,qDotFull);
    } else {
        return mathUtils::trapz(grid.x,grid.x*qDot)/mathUtils::trapz(grid.x,qDot);
    }
}

double flameSys::targetFlamePosition(double t)
{
    return (t <= options.rFlameT0) ? options.rFlameInitial
        :  (t >= options.rFlameT0+options.rFlameDt) ? options.rFlameFinal
        : options.rFlameInitial + (options.rFlameFinal-options.rFlameInitial)*(t-options.rFlameT0)/options.rFlameDt;
}

void flameSys::update_rStag(const double t, const bool updateIntError)
{
    rFlameActual = getFlamePosition();
    rFlameTarget = targetFlamePosition(t);
    if (updateIntError) {
        flamePosIntegralError += (rFlameTarget-rFlameActual)*(t-tFlamePrev);
        tFlamePrev = t;
    }

    options.rStag = options.rFlameProportionalGain *
        ( (rFlameTarget-rFlameActual) + (flamePosIntegralError + (rFlameTarget-rFlameActual)*(t-tFlamePrev))*options.rFlameIntegralGain );

    if (debugParameters::debugFlameRadiusControl) {
        cout << "rFlameControl: " << "rF = " << rFlameActual << "   rStag = " << options.rStag;
        cout << "   P = " <<  options.rFlameProportionalGain*(rFlameTarget-rFlameActual);
        cout << "   I = " << options.rFlameProportionalGain*flamePosIntegralError*options.rFlameIntegralGain;
        cout << "  dt = " << t-tFlamePrev << endl;
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

    int jj = nPoints - 1;

    try {
        updateThermoProperties();
        updateTransportProperties();
        if (options.singleCanteraObject) {
            simpleGas.getReactionRates(wDot);
        } else {
            gas.getReactionRates(wDot);
        }
    } catch (Cantera::CanteraError) {
        cout << "Error evaluating thermodynamic properties" << endl;
        writeStateMatFile("errorOutput",true);
        throw;
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
    matlabFile jacobianFile(options.outputDir+"/jacobianComparison.mat");
    jacobianFile.writeArray2D("aFirst",firstBlock1);
    jacobianFile.writeArray2D("aMiddle",secondBlock1);
    jacobianFile.writeArray2D("aEnd",thirdBlock1);
    jacobianFile.writeArray2D("nFirst",firstBlock2);
    jacobianFile.writeArray2D("nMiddle",secondBlock2);
    jacobianFile.writeArray2D("nEnd",thirdBlock2);
    jacobianFile.close();

    cout << "Wrote output file: " << options.outputDir+"/jacobianComparison.mat" << endl;
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
