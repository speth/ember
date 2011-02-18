#include "test_qssintegrator.h"

#include <boost/format.hpp>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using boost::format;

void QSSTestProblem::odefun(double t, const dvector& y, dvector& q, dvector& d)
{
    // csdfe(y, q, d, t)

    // description:
    // derivative function evaluator(gsub) for an atmospheric chemical
    // relaxation test problem involving cesium and cesium ions. format-
    // ion and loss rates are calculated for this set of "stiff ordinary
    // differential equations" that was suggested by by d. edelson of
    // bell laboratories.

    // argument list definitions:
    // y(i)         r*4  current values of the functions plus the     i/o
    //                   extra data at the end of the array that may be
    //                   passed back and forth between "csdfe" and the
    //                   main program. locations in y(i) which represent
    //                   the functions being advanced should not be
    //                   tampered with here.
    // q(i)         r*4  total formation rates.                         i
    // d(i)         r*4  total loss rates.                              i
    // t            r*4  the value of the independent variable.         i

    // utilize local storage for varibles.
    double o2m = y[0];
    double csp = y[1];
    double cs = y[2];
    double cso2 = y[3];
    double o2 = y[4];
    double n2 = y[5];

    // calculate electron density for local use and transmission back to
    // the main program via y(7). however in this case this value should
    // not be trusted since "chemeq" will not call the "gsub" with the
    // latest function values after the final step has converged. y(7)
    // will be one iteration behind in this case. y(7) and y(6) are
    // examples tho, of how data may be transfered between the "gsub" and
    // the main program.
    double ne = std::max(csp - o2m, 0.0);
//    y[6] = ne;

    // calculate reaction rates.
    double cr1 = 5.00e-08*o2m*csp;
    double cr2 = 1.00e-12*csp*ne;
    double cr3 = 3.24e-03*cs;
    double cr4 = 4.00e-01*o2m;
    double cr5 = 1.00e-31*o2*cs*(cs + cso2 + n2 + o2);
    double cr6 = 1.24e-30*o2*o2*ne;
    double cr7 = 1.00e-31*o2*n2*ne;

    // calculate total formation rates (c(i)) and total loss rates (d(i))
    // for each species.

    // o2m
    q[0] = cr6 + cr7;
    d[0] = cr1 + cr4;

    // cs+
    q[1] = cr3;
    d[1] = cr1 + cr2;

    // cs
    q[2] = cr1 + cr2;
    d[2] = cr3 + cr5;

    // cso2
    q[3] = cr5;
    //      q(4) = q(4) - 1.00e-31*o2*cs*cso2
    //      d(4) = - 1.00e-31*o2*cs*cso2

    // o2
    q[4] = cr1 + cr4;
    d[4] = cr5 + cr6 + cr7;
}

int main(int argc, char** argv) {
    // This is the dirver program for the seven-species cesium
    // mechanism test problem.  The code integrates the system
    // MXCASE times using different values of the chemeq2 variable
    // epsmin (set by passing an entry from array EPS through
    // CHEMSP before each integration).

    QSSTestProblem qssSolver;

    // PROGRAM SPECIFICATIONS.
    dvector Y(10);
    dvector YF(10);
    dvector YMIN(10, 1e-20);
    dvector YI(10);
    dvector epsil(10);

    vector<std::string> SPSYM(7);
    SPSYM[0] = "O2-";
    SPSYM[1] = "CS+";
    SPSYM[2] = "CS";
    SPSYM[3] = "CSO2";
    SPSYM[4] = "O2";
    SPSYM[5] = "N2";
    SPSYM[6] = "NE";

    // For this example, the external subroutine that calculates the
    // source terms is called CSDFE.
    int MXCASE = 9;

    double EPS[15] = {0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005,
                      0.00001, 0.000005, 0.000001, 5e-7, 1e-7, 5e-8, 1e-8};

//1000     FORMAT('CASE NO. ', I5, '    PARAMETERS;', /,
//   .    ' CONVERGENCE PARAMETER EPS = ', 1PE10.3, /,
//   .    ' INNER LOOP LENGTH;', I5)
//1001     FORMAT(/, '    SPECIE    Y - INITAL      Y - FINAL ',
//   .    '  Y - SOLUTION   REL ERR')
//1002     FORMAT(5X, A4, 1P3E15.6, E10.3)
//1003     FORMAT(/, ' T - INITIAL = (', 1PE10.3, ') T - FINAL = (',
//   .    E10.3, ')')
//1004     FORMAT(/' INTEGRATION STATISTICS;')
//1005     FORMAT(' CPU TIME USED FOR INTEGRATION;', 1PE10.3,
//   .    ' SEC.,  CPU TIME NORMALIZED;', I8)
//1006     FORMAT(' SUM OF THE RELATIVE ERRORS SQUARED; ', 1PE10.3)
//1007     FORMAT(/)

    //     Note that the timing routines included may not work on
    //     all systems.  Extra timing options are included as comments.

    // INITIALIZE CONTROL PARAMETERS.

    //     INLP allows the user to subdivide the interval over which
    //     each test is run.  For INLP=1, CHEMEQ2 is sent the full
    //     interval TF-TI (specified below) as the global timestep.
    int INLP = 1;

    // For this particular test, the electron number density is not
    // integrated.  The other five reacting species are integrated,
    // and the electron density is found through charge conservation.
    // This calculation is done within CSDFE.  Therefore, NA = 5 is
    // the number of equations that are integrated, but NS = 7 is the
    // number of species.  Species to be integrated must be placed in
    // first NA positions within the Y array.  CHEMEQ2 only works with
    // these first NA entries since NA is passed in the argument list
    // below, but all NS values are available to and used by CSDFE.

    int NS = 7;
    // int NA = 5;

    // "TI"  -  INITIAL TIME,   "TF"  -  FINAL TIME.
    double TI = 0.0;
    double TF = 1000.0;
    double DELTAT = (TF - TI)/INLP;

    // STORE INITIAL(TI = 0.0) AND FINALT(F = 1000.0) VALUES.
    // O2-
    YI[0] = 5.200e+02;
    YF[0] = 2.59139492061e+04;

    // CS+
    YI[1] = 6.200e+02;
    YF[1] = 7.55718460300e+04;

    // CS
    YI[2] = 1.000e+12;
    YF[2] = 1.53194051722e+03;

    // CSO2
    YI[3] = 0;
    YF[3] = 9.99999923516e+11;

    // O2
    YI[4] = 3.600e+14;
    YF[4] = 3.59000000051e+14;

    // N2
    YI[5] = 1.400e+15;
    YF[5] = 1.40000000000e+15;

    // NE
    YI[6] = 1.000e+02;
    YF[6] = 4.96578968239e+04;

    // LOOP OVER THE TEST CASES.
    for (int ICASE=0; ICASE<MXCASE; ICASE++) {
        cout << ICASE << ", " << EPS[ICASE] << ", " << INLP << endl;
        qssSolver.epsmin = EPS[ICASE];
        qssSolver.itermax = 5;
        qssSolver.ymin = YMIN;

        // RESET "Y" TO INITIAL VALUES "YI".
        for (int i=0; i<NS; i++) {
            Y[i] = YI[i];
        }

        // INNER LOOP TO DETERMINE OVERHEAD OR RELATIVE STARTING EFFECIENCY
        // OF ITEGRATION SCHEME BEING TESTED.
        for (int istep=0; istep<INLP; istep++) {
            // CALL INTEGRATOR.
            // CALL CHEMEQ2(DELTAT, CSDFE, NA, Y)
            qssSolver.initialize(Y, TI);
            qssSolver.integrateToTime(DELTAT);
        }

        Y = qssSolver.y;
        // Calculate final electron density from densities of other charges species
        Y[6] = Y[1] - Y[0];

        // CALCULATE RELATIVE ERROR.
        double sum = 0.0;
        for (int i=0; i<NS; i++) {
            epsil[i] = std::abs(Y[i] - YF[i])/std::min(Y[i] , YF[i]);
            sum += epsil[i]*epsil[i];
        }

        // Root-mean-square error is calculated using ns-1 (rather than ns)
        // since N2 is inert.

        sum = sqrt(sum/(NS-1));
        // PRINT RESULTS.

        cout << format("t - tInitial = %f ;   t - tFinal = %f\n") % TI % TF;
        cout << "Species        Y-Inital     Y-Final    Y-Solution     Rel. Error\n";
        for (int i=0; i<NS; i++) {
            cout << format(" %6s   %012.6e    %012.6e    %012.6e    %012.6e\n") %
                    SPSYM[i] % YI[i] % YF[i] % Y[i] % epsil[i];
        }
        cout << "Integration Statistics:" << endl;
//        WRITE(LO, 1006) SUM
//        WRITE(LO,  1005) CPUT, TNORM
        cout << format("Sum of the relative errors squared: %12.6e") % sum << endl;
//        WRITE(*,699)  EPS(ICASE),
//   &         CPUT,
//   &         TNORM
//    // &         INT(CPUT*1024. + .5)
//   &         ,sum
//699      format(1x,25HEPS, time, ticks, error: ,E7.1,2x,e10.4,2x,
//   &         I5,2x,e10.4)
//        WRITE(LO, 1007)
//        CALL CHEMCT(TF)
    }

    return 0;
}
