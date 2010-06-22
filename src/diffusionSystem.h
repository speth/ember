#pragma once

#include "mathUtils.h"
#include "sundialsUtils.h"
#include "chemistry0d.h"

class DiffusionSystem : public sdODE, public GridBased
{
    // This is a system representing diffusion of a single solution component
private:
    // Jacobian data
    sdBandMatrix* bandedJacobian;
    vector<int> pMat;

    // Sundials solver parameters
    sdVector* abstol;
    double reltol;

    // The sum of the terms held constant for this system
    dvector constantTerm;
};

class SpeciesDiffusionSystem : public DiffusionSystem
{
    // This system represents the diffusion of a single species at all grid points
public:
    int f(realtype t, sdVector& y, sdVector& ydot);

    int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                            sdVector& resIn, realtype c_j);

    int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                            sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);

    // Diffusion coefficients
    dvector rhoD; // density-weighted, mixture-averaged diffusion coefficients [kg/m*s] (= rho*Dkm)
    dvector Dkt; // thermal diffusion coefficients [kg/m*s]

    // Diffusion mass fluxes
    dvector jFick; // Normal diffusion (Fick's Law) [kg/m^2*s]
    dvector jSoret; // Soret Effect diffusion [kg/m^2*s]
    dvector dYkdx; // upwinded
};


class TemperatureDiffusionSystem : public DiffusionSystem
{
    // This system represents the diffusion Temperature at all grid points
public:
    int f(realtype t, sdVector& y, sdVector& ydot);

    int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                            sdVector& resIn, realtype c_j);

    int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                            sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);

private:
    dvector sumcpj; // for enthalpy flux term [W/m^2*K]
    dvector lambda; // thermal conductivity [W/m*K]
    dvector cp; // specific heat capacity (average) [J/kg*K]
    Array2D cpSpec; // species specific heat capacity [J/mol*K]
    dvector W; // species molecular weights [kg/kmol]
    dvector qFourier; // heat flux [W/m^2]
    dvector dTdx; // upwinded
    dvector dTdxCen; // centered difference (for enthalpy flux term)
};


class MomentumDiffusionSystem : public DiffusionSystem
{
    // This system represents the diffusion of tangential momentum at all grid points
public:
    int f(realtype t, sdVector& y, sdVector& ydot);

    int preconditionerSetup(realtype t, sdVector& yIn, sdVector& ydotIn,
                            sdVector& resIn, realtype c_j);

    int preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, sdVector& resIn,
                            sdVector& rhs, sdVector& outVec, realtype c_j, realtype delta);


private:
    dvector mu; // viscosity [Pa*s]
    dvector dUdx; // upwinded
};

