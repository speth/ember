#include "profileGenerator.h"
#include "chemistry0d.h"
#include "strainFunction.h"
#include "dataFile.h"

void ProfileGenerator::setOptions(configOptions& _options)
{
    options = _options;
}

void ProfileGenerator::generateProfile(void)
{
    cout << "Generating initial profiles from given fuel and oxidizer compositions." << endl;

    // Set up a CanteraGas object to use for calculating the initial profiles
    CanteraGas gas;
    gas.mechanismFile = options.gasMechanismFile;
    gas.phaseID = options.gasPhaseID;
    gas.pressure = options.pressure;
    gas.initialize(false);
    size_t nSpec = gas.nSpec;

    // Create a uniform initial grid
    x.resize(options.nPoints);
    double xRight, xLeft;
    xRight = options.xRight;
    if (options.twinFlame || options.curvedFlame) {
        x[0] = 0;
        xLeft = options.centerGridMin;
    } else {
        xLeft = options.xLeft;
        x[0] = xLeft;
    }
    // Uniform initial grid
    double dx = (xRight-xLeft)/((double) nPoints-1);
    for (size_t j=1; j<nPoints; j++) {
        x[j] = xLeft + j * dx;
    }

    grid.unburnedLeft = options.unburnedLeft;
    grid.updateBoundaryIndices();
    grid.updateValues();

    int jm = (grid.ju+grid.jb)/2; // midpoint of the profiles.
    int jl = jm - 4;
    int jr = jm + 4;

    Yb.resize(nSpec);
    Yu.resize(nSpec);

    double a = options.strainRateInitial;

    Tu = options.Tu;

    // Reactants
    dvector reactants = calculateReactantMixture();
    gas.thermo.setState_TPX(Tu, gas.pressure, &reactants[0]);
    rhou = gas.thermo.density();
    gas.thermo.getMassFractions(&Yu[0]);
    gas.thermo.getMassFractions(&Y(0,grid.ju));

    // Products
    gas.thermo.setState_TPX(Tu, gas.pressure, &reactants[0]);
    Cantera::equilibrate(gas.thermo,"HP");
    Tb = gas.thermo.temperature();
    rhob = gas.thermo.density();
    gas.thermo.getMassFractions(&Yb[0]);
    gas.thermo.getMassFractions(&Y(0,grid.jb));

    // Diluent in the middle
    gas.thermo.setState_TPY(Tu, gas.pressure, options.oxidizer);
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


    for (int j=1; j<jl; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0);
        }
        T[j] = T[0];
    }

    for (int j=jl; j<jm; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,0) + (Y(k,jm)-Y(k,0))*(x[j]-x[jl])/(x[jm]-x[jl]);
        }
        T[j] = T[0] + (T[jm]-T[0])*(x[j]-x[jl])/(x[jm]-x[jl]);
    }

    for (int j=jm+1; j<jr; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,jm) + (Y(k,grid.jj)-Y(k,jm))*(x[j]-x[jm])/(x[jr]-x[jm]);
        }
        T[j] = T[jm] + (T[grid.jj]-T[jm])*(x[j]-x[jm])/(x[jr]-x[jm]);
    }

    for (size_t j=jr; j<nPoints; j++) {
        for (size_t k=0; k<nSpec; k++) {
            Y(k,j) = Y(k,grid.jj);
        }
        T[j] = T[grid.jj];
    }

    dvector yTemp(nPoints);
    for (size_t k=0; k<nSpec; k++) {
        for (size_t j=0; j<nPoints; j++) {
            yTemp[j] = Y(k,j);
        }

        for (size_t i=0; i<10; i++) {
            mathUtils::smooth(yTemp);
        }

        for (size_t j=0; j<nPoints; j++) {
            Y(k,j) = yTemp[j];
        }
    }

    for (size_t i=0; i<5; i++) {
        mathUtils::smooth(T);
    }

    // Grid and initial profiles of T, U and V
    dvector rho(nPoints);
    for (size_t j=0; j<nPoints; j++) {
        gas.setStateMass(&Y(0,j), T[j]);
        rho[j] = gas.getDensity();
        U[j] = a*sqrt(rhou/rho[j]);
    }

    for (int i=0; i<2; i++) {
        mathUtils::smooth(U);
    }

    if (options.fixedLeftLoc)
    {
        jm = 0;
    }

    dvector V(nPoints);
    V[jm] = 0;
    for (size_t j=jm+1; j<nPoints; j++) {
        V[j] = V[j-1] - rho[j]*U[j]*(x[j]-x[j-1]);
    }

    for (size_t j=jm-1; j>=0; j--) {
        V[j] = V[j+1] + rho[j]*U[j]*(x[j+1]-x[j]);
    }
}

void ProfileGenerator::loadProfile(void)
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

    U = infile.readVector("U");
    T = infile.readVector("T");
    Y = infile.readArray2D("Y");
    t = infile.readScalar("t");

    if (!options.fileNumberOverride) {
        options.outputFileNumber = (int) infile.readScalar("fileNumber");
    }

    infile.close();

    if (options.overrideTu) {
        T[grid.ju] = options.Tu;
    }
    Tu = T[grid.ju];

    CanteraGas gas;
    gas.mechanismFile = options.gasMechanismFile;
    gas.phaseID = options.gasPhaseID;
    gas.pressure = options.pressure;
    gas.initialize(false);
    size_t nSpec = gas.nSpec;

    // save the burned gas properties for the case where burned values are not fixed
    double TbSave = T[grid.jb];
    dvector YbSave(nSpec);
    for (size_t k=0; k<nSpec; k++) {
        YbSave[k] = Y(k,grid.jb);
    }

    if (options.overrideReactants) {
        dvector reactants = calculateReactantMixture();
        gas.thermo.setState_TPX(Tu,gas.pressure, &reactants[0]);
        gas.thermo.getMassFractions(&Y(0,grid.ju));
        gas.thermo.setState_TPX(Tu,gas.pressure, &reactants[0]);
        Cantera::equilibrate(gas.thermo,"HP");
        gas.thermo.getMassFractions(&Y(0,grid.jb));
        T[grid.jb] = gas.thermo.temperature();
    }
    Tb = T[grid.jb];

    gas.setStateMass(&Y(0,grid.ju), T[grid.ju]);
    rhou = gas.getDensity();
    gas.setStateMass(&Y(0,grid.jb), T[grid.ju]);
    rhob = gas.getDensity();
    Yu.resize(nSpec); Yb.resize(nSpec);
    for (size_t k=0; k<nSpec; k++) {
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
        for (size_t k=0; k<nSpec; k++) {
            Y(k,grid.jb) = YbSave[k];
        }
    }

    // TODO: Figure out where this actually should be implemented
//    updateLeftBC();
//
//    double controlSignal;
//    if (grid.leftBoundaryConfig == grid.lbControlVolume && options.xFlameControl) {
//        if (alpha == 0) {
//            controlSignal = V[0]/rhoLeft;
//        } else {
//            double tmp = pow(x[0],2) + 2*rV[0]/rhoLeft;
//            controlSignal = sign(tmp)*sqrt(abs(tmp));
//        }
//        flamePosIntegralError = controlSignal/(options.xFlameProportionalGain*options.xFlameIntegralGain);
//    }
}

dvector ProfileGenerator::calculateReactantMixture(void)
{
    // Calculate the composition of the reactant mixture from compositions of
    // the fuel and oxidizer mixtures and the equivalence ratio.

    Cantera_CXX::IdealGasMix fuel(options.gasMechanismFile,options.gasPhaseID);
    Cantera_CXX::IdealGasMix oxidizer(options.gasMechanismFile,options.gasPhaseID);

    fuel.setState_TPX(options.Tu, options.pressure, options.fuel);
    oxidizer.setState_TPX(options.Tu, options.pressure, options.oxidizer);


    double Cf(0), Hf(0), Of(0); // moles of C/H/O in fuel
    double Co(0), Ho(0), Oo(0); // moles of C/H/O in oxidizer

    int nSpec = fuel.nSpecies();
    int mC = fuel.elementIndex("C");
    int mO = fuel.elementIndex("O");
    int mH = fuel.elementIndex("H");

    dvector Xf(nSpec), Xo(nSpec), Xr(nSpec);
    fuel.getMoleFractions(&Xf[0]);
    oxidizer.getMoleFractions(&Xo[0]);
    dvector a(fuel.nElements());
    for (int k=0; k<nSpec; k++) {
        fuel.getAtoms(k,&a[0]);
        Cf += a[mC]*Xf[k];
        Co += a[mC]*Xo[k];
        Hf += a[mH]*Xf[k];
        Ho += a[mH]*Xo[k];
        Of += a[mO]*Xf[k];
        Oo += a[mO]*Xo[k];
    }
    double stoichAirFuelRatio = -(Of-2*Cf-Hf/2)/(Oo-2*Co-Ho/2);
    Xr = Xf*options.equivalenceRatio + stoichAirFuelRatio*Xo;
    Xr /= mathUtils::sum(Xr);

    return Xr;
}
