#include "qssintegrator.h"

#include <iostream>

using std::cout;
using std::endl;
using std::abs;
using mathUtils::sign;

QSSIntegrator::QSSIntegrator(size_t Neq)
    : y(Neq)
    , N(Neq)
{
    t = 0;
    epsmax = 10;
    epsmin = 1e-2;
    dtmin = 1e-15;
    tstart = 0;
    itermax = 1;

    gcount = 0;
    rcount = 0;

    ymin.resize(N, 1e-20);
    ym1.resize(N);
    ym2.resize(N);
    scrarray.resize(N);
    stabilityCheck = false;

    q.resize(N, 0);
    d.resize(N, 0);
    rtaus.resize(N);
    y1.resize(N);

    ys.resize(N);
    rtau.resize(N);
    qs.resize(N);
}

int QSSIntegrator::integrateToTime(double dtg)
{
    int iter;

    static const double tfd = 1.000008;
    double ts;

    double eps = 1e-10;
    // double rswitch = 5.965900; // rswitch for 4-4 pade: 5.9659

    // Initialize the control parameters.
    double tn = 0.0;

    // Store and limit to 'ymin' the initial values.
    q.assign(N, 0);
    d.assign(N, 0);
    for (size_t i=0; i<N; i++) {
        y[i] = std::max(y[i], ymin[i]);
    }

    // Evaluate the derivatives of the initial values.
    odefun(tn + tstart, y, q, d);
    gcount += 1;

    // Estimate the initial stepsize.

    // strongly increasing functions(q >>> d assumed here) use a step-
    // size estimate proportional to the step needed for the function to
    // reach equilibrium where as functions decreasing or in equilibrium
    // use a stepsize estimate directly proportional to the character-
    // istic stepsize of the function. convergence of the integration
    // scheme is likely since the smallest estimate is chosen for the
    // initial stepsize.

    double scrtch = 1.0e-25;

    for (size_t i=0; i<N; i++) {
        double absq = abs(q[i]);
        double scr2 = abs(1/y[i]) * sign(0.1*epsmin*absq - d[i]);
        double scr1 = scr2 * d[i];
        scrtch = std::max(std::max(scr1, -abs(absq-d[i])*scr2), scrtch);
    }

    double sqreps = 0.5;
    double dt = std::min(sqreps/scrtch, dtg);

    bool newStep = true;

    while (true) {
        if (newStep) {
        // The starting values are stored.
            ts = tn;

            for (size_t i=0; i<N; i++) {
                rtau[i] = dt*d[i]/y[i];
                ys[i] = y[i];
                qs[i] = q[i];
                rtaus[i] = rtau[i];
            }
        }

        // Find the predictor terms.
        for (size_t i=0; i<N; i++) {
            // prediction
            double rtaui = rtau[i];

            // note that one of two approximations for alpha is chosen:
            // 1) Pade b for all rtaui (see supporting memo report)
            //      or
            // 2) Pade a for rtaui<=rswitch,
            //    linear approximation for rtaui > rswitch
            //    (again, see supporting NRL memo report (Mott et al., 2000) )

            // Option 1): Pade b

           double alpha = (180+rtaui*(60+rtaui*(11+rtaui))) / (360 + rtaui*(60 + rtaui*(12 + rtaui)));

           // Option 2): Pade a or linear

           //c         if(rtaui.le.rswitch) then
           //c            alpha = (840.+rtaui*(140.+rtaui*(20.+rtaui)))
           //c     &           / (1680. + 40. * rtaui*rtaui)
           //c         else
           //c            alpha = 1.-1./rtaui
           //c         end if

           scrarray[i] = (q[i]-d[i])/(1.0 + alpha*rtaui);
        }

        iter = 1;
        while (iter < itermax) {
            // limit decreasing functions to their minimum values.
            if (stabilityCheck) {
                for (size_t i=0; i<N; i++) {
                    ym2[i] = ym1[i];
                    ym1[i] = y[i];
                }
            }

            for (size_t i=0; i<N; i++) {
                y[i] = std::max(ys[i] + dt*scrarray[i], ymin[i]);
            }

            if (iter == 1) {
                // The first corrector step advances the time (tentatively) and
                // saves the initial predictor value as y1 for the timestep check later
                tn = ts + dt;
                for (size_t i=0; i<N; i++) {
                    y1[i] = y[i];
                }
            }

            // Evaluate the derivatives for the corrector.
            odefun(tn + tstart, y, q, d);
            gcount += 1;
            eps = 1.0e-10;

            for (size_t i=0; i<N; i++) {
               double rtaub = 0.5*(rtaus[i]+dt*d[i]/y[i]);

               // Same options for calculating alpha as in predictor:

               // Option 1): Pade b
               double alpha = (180.+rtaub*(60.+rtaub*(11.+rtaub))) / (360. + rtaub*(60. + rtaub*(12. + rtaub)));

               //c  Option 2):  Pade a or linear
               //c
               //c             if(rtaub.le.rswitch) then
               //c                alpha = (840.+rtaub*(140.+rtaub*(20.+rtaub)))
               //c     &               / (1680. + 40.*rtaub*rtaub)
               //c             else
               //c                alpha = 1.-1./rtaub
               //c             end if

               double qt = qs[i]*(1.0 - alpha) + q[i]*alpha;
               double pb = rtaub/dt;
               scrarray[i] = (qt - ys[i]*pb) / (1.0 + alpha*rtaub);
            }
            iter += 1;
        }

        // Calculate new f, check for convergence, and limit decreasing
        // functions. the order of the operations in this loop is important.
        for (size_t i=0; i<N; i++) {
            double scr2 = std::max(ys[i] + dt*scrarray[i], 0.0);
            double scr1 = abs(scr2 - y1[i]);
            y[i] = std::max(scr2, ymin[i]);

            if (0.25*(ys[i] + y[i]) > ymin[i]) {
               scr1 = scr1/y[i];
               eps = std::max(.5*(scr1+ std::min(abs(q[i]-d[i])/(q[i]+d[i]+1e-30), scr1)),eps);
            }
        }

        if (stabilityCheck) {
            for (size_t i=0; i<N; i++) {
                ym2[i] = ym1[i];
                ym1[i] = y[i];
            }
        }

        eps /= epsmin;

        // Print out dianostics if stepsize becomes too small.
        // write(*,*) tn, dt, dtmin, tn, dtmin+1.0e-4*tn
        if (dt <= dtmin + 1e-16*tn) {
            //write(lo, 1003) dt, tn, dtmin
            for (size_t i=0; i<N; i++) {
                double dtc = epsmin*y[i]/(abs(q[i]-d[i]) + 1.0e-30);
               //write(lo, 1004) q(i), d(i), y(i), rtau(i), dtc, q(i)-d(i),ys(i), y0(i), ymin(i)
            }

            //1003     format('1    chemeq error;   stepsize too small ! ! !', /,
            //   1    '     dt = ', 1pe10.3, ' tn = ', d25.15,
            //   2    ' dtmin = ',e10.3, //, 11x, 'q', 10x, 'd', 10x, 'y',
            //   3         8x, 'rtau', 8x, 'dtc', 7x, 'q - d',7x, 'ys',
            //   4         9x, 'y0', 8x, 'ymin')
            //1004     format(5x, 1p12e11.3)
            dt = dtg - ts;
            dt = std::min(dtmin, abs(dt));

            // return, indicating an error
            return 1;
        }

        // Check for convergence.

        // The following section is used for the stability check
        double stab = 0;
        if (stabilityCheck) {
            stab = 0.01;
            if (itermax >= 3) {
                for (size_t i=0; i<N; i++) {
                    stab = std::max(stab, abs(y[i]-ym1[i])/(abs(ym1[i]-ym2[i])+1e-20*y[i]));
                }
            }
        }

        if (eps <= epsmax && stab <= 1.0) {
            // Valid step.  Return if dtg has been reached.
            if (dtg <= tn*tfd) {
                return 0;
            }
        } else {
            // Invalid step; reset tn to ts
            tn = ts;
        }

        // perform stepsize modifications.
        // estimate sqrt(eps) by newton iteration.

        double rteps = 0.5*(eps + 1.0);
        rteps = 0.5*(rteps + eps/rteps);
        rteps = 0.5*(rteps + eps/rteps);

        double dto = dt;
        dt = std::min(dt*(1.0/rteps + 0.005), tfd*(dtg - tn));
        if (stabilityCheck) {
            dt = std::min(dt, dto/(stab+.001));
        }

        // begin new step if previous step converged.
        if (eps > epsmax || stab > 1) {
            rcount += 1;

            // After an unsuccessful step the initial timescales don't
            // change, but dt does, requiring rtaus to be scaled by the
            // ratio of the new and old timesteps.
            dto = dt/dto;

            for (size_t i=0; i<N; i++) {
                rtaus[i] = rtaus[i] * dto;
            }

            // After an unsuccessful step, do not recalculate the initial source terms
            newStep = false;
        } else {
            // Successful step; get the source terms for the next step
            odefun(t, y, q, d);
            gcount = gcount + 1;
            newStep = true;
        }
    }
}
