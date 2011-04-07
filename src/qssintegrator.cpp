#include "qssintegrator.h"
#include "debugUtils.h"

#include <iostream>

using std::abs;
using mathUtils::sign;

QSSIntegrator::QSSIntegrator()
{
    N = 0;
    epsmax = 20;
    epsmin = 1e-2;
    dtmin = 1e-15;
    dtmax = 1e-6;
    tstart = 0;
    itermax = 2;
    tfd = 1.000008;
    abstol = 1e-8;

    stabilityCheck = true;

    firstStep = true;
    dt = 0;
    ts = 0;
    gcount = 0;
    rcount = 0;
    tn = 0;

    debug = false;
}

void QSSIntegrator::initialize(size_t N_)
{
    N = N_;
    y.resize(N);

    q.assign(N, 0);
    d.assign(N, 0);
    rtaus.resize(N);
    y1.resize(N);

    ys.resize(N);
    rtau.resize(N);
    qs.resize(N);

    ymin.resize(N, 1e-20);
    enforce_ymin.resize(N, true);
    ym1.resize(N);
    ym2.resize(N);
    scratch.resize(N);
}

void QSSIntegrator::setState(const dvector& yIn, double tstart_)
{
    assert(yIn.size() == N);

    // Store and limit to 'ymin' the initial values.
    for (size_t i=0; i<N; i++) {
        if (enforce_ymin[i]) {
            y[i] = std::max(yIn[i], ymin[i]);
        }
    }

    gcount = 0;
    rcount = 0;
    tstart = tstart_;
    tn = 0.0;
    firstStep = true;
}

void QSSIntegrator::getInitialStepSize(double tf)
{
    // Estimate the initial stepsize.

    // strongly increasing functions(q >>> d assumed here) use a step-
    // size estimate proportional to the step needed for the function to
    // reach equilibrium where as functions decreasing or in equilibrium
    // use a stepsize estimate directly proportional to the character-
    // istic stepsize of the function. convergence of the integration
    // scheme is likely since the smallest estimate is chosen for the
    // initial stepsize.

    firstStep = false;
    double scrtch = 1.0e-25;

    for (size_t i=0; i<N; i++) {
        if (abs(y[i]) > abstol) {
            double absq = abs(q[i]);
            double scr2 = abs(1/y[i]) * sign(0.1*epsmin*absq - d[i]);
            double scr1 = scr2 * d[i];
            scrtch = std::max(std::max(scr1, -abs(absq-d[i])*scr2), scrtch);
        }
    }

    double sqreps = 0.5;
    dt = std::min(sqreps/scrtch, tf);
    dt = std::min(dt, dtmax);
}

int QSSIntegrator::integrateToTime(double tf)
{
    while (tfd*tn < tf) {
        int ret = integrateOneStep(tf);
        if (ret) {
            // Integration failed
            return ret;
        }
    }
    return 0;
}

int QSSIntegrator::integrateOneStep(double tf) {
    // Evaluate the derivatives at the initial state.
    odefun(tn + tstart, y, q, d);
    gcount += 1;

    if (firstStep) {
        getInitialStepSize(tf);
    }

    // Store starting values
    ts = tn;
    for (size_t i=0; i<N; i++) {
        if (enforce_ymin[i]) {
            assert(d[i] == 0);
            rtau[i] = dt*d[i]/y[i];
        } else {
            rtau[i] = 0;
        }
        ys[i] = y[i];
        qs[i] = q[i];
        rtaus[i] = rtau[i];
    }

    // double rswitch = 5.965900; // rswitch for 4-4 pade: 5.9659

    // Repeat integration until a successful timestep has been taken
    while (true) {
        // Find the predictor terms.
        for (size_t i=0; i<N; i++) {
            // prediction
            double rtaui = rtau[i];
            double alpha = (180+rtaui*(60+rtaui*(11+rtaui))) /
                           (360 + rtaui*(60 + rtaui*(12 + rtaui)));
           scratch[i] = (q[i]-d[i])/(1.0 + alpha*rtaui);
        }

        double eps = 1e-10;
        for (int iter=0; iter<itermax; iter++) {
            // limit decreasing functions to their minimum values.
            if (stabilityCheck) {
                ym2 = ym1;
                ym1 = y;
            }

            for (size_t i=0; i<N; i++) {
                y[i] = std::max(ys[i] + dt*scratch[i], ymin[i]);
            }

            if (iter == 0) {
                // The first corrector step advances the time (tentatively) and
                // saves the initial predictor value as y1 for the timestep check later
                tn = ts + dt;
                y1 = y;
            }

            // Evaluate the derivatives for the corrector.
            odefun(tn + tstart, y, q, d, true);
            gcount += 1;
            eps = 1.0e-10;

            double rtaub;
            for (size_t i=0; i<N; i++) {
                if (enforce_ymin[i]) {
                    rtaub = 0.5*(rtaus[i]+dt*d[i]/y[i]);
                } else {
                    rtaub = 0;
                }
                double alpha = (180.+rtaub*(60.+rtaub*(11.+rtaub))) /
                               (360. + rtaub*(60. + rtaub*(12. + rtaub)));
                double qt = qs[i]*(1.0 - alpha) + q[i]*alpha;
                scratch[i] = (qt - ys[i]*rtaub/dt) / (1.0 + alpha*rtaub);
            }
        }

        // Calculate new f, check for convergence, and limit decreasing
        // functions. the order of the operations in this loop is important.
        for (size_t i=0; i<N; i++) {
            double scr2 = std::max(ys[i] + dt*scratch[i], 0.0);
            double scr1 = abs(scr2 - y1[i]);
            y[i] = std::max(scr2, ymin[i]);

            if (abs(y[i]) > abstol && 0.25*(ys[i] + y[i]) > ymin[i]) {
               scr1 = scr1/y[i];
               eps = std::max(.5*(scr1+ std::min(abs(q[i]-d[i])/(q[i]+d[i]+1e-30), scr1)),eps);
            }
        }

        if (stabilityCheck) {
            ym2 = ym1;
            ym1 = y;
        }

        eps /= epsmin;

        if (dt <= dtmin + 1e-16*tn) {
            // Integration failed with a timestep that is too small.
            logFile.write("QSSIntegrator failed: timestep too small.");
            logFile.write(format("dt = %e, tn = %e, dtmin = %e") % dt % tn % dtmin);
            logFile.write("     i    q[i]     d[i]     y[i]    rtau[i]    dtc    q[i]-d[i]    ys[i]    ymin[i]");

            for (size_t i=0; i<N; i++) {
                double dtc = epsmin*y[i]/(abs(q[i]-d[i]) + 1.0e-30);
                logFile.write(format("%3i %9e %9e %9e %9e %9e %9e %9e %9e") % i % q[i] %
                        d[i] % y[i] % rtau[i] % dtc % (q[i]-d[i]) % ys[i] % ymin[i]);
            }

            // Return, indicating an error.
            return 1;
        }

        // Check for convergence.

        // The following section is used for the stability check.
        double stab = 0;
        if (stabilityCheck && itermax >= 3) {
            stab = 0.01;
            for (size_t i=0; i<N; i++) {
                if (abs(y[i]) > abstol) {
                    stab = std::max(stab, abs(y[i]-ym1[i])/(abs(ym1[i]-ym2[i])+1e-20*y[i]));
                }
            }
        }

        if (eps <= epsmax && stab <= 1.0) {
            // Successful step. Return if tf has been reached.
            if (tf <= tn*tfd) {
                return 0;
            }
        } else {
            // Unsuccessful step; reset tn to ts
            tn = ts;
        }

        // Perform stepsize modifications.

        // estimate sqrt(eps) by newton iteration.
        double rteps = 0.5*(eps + 1.0);
        rteps = 0.5*(rteps + eps/rteps);
        rteps = 0.5*(rteps + eps/rteps);

        double dto = dt;
        dt = std::min(dt*(1.0/rteps + 0.005), tfd*(tf - tn));
        dt = std::min(dt, dtmax);
        if (stabilityCheck) {
            dt = std::min(dt, dto/(stab+.001));
        }

        if (eps > epsmax || stab > 1) {
            // After an unsuccessful step, the initial timescales don't
            // change, but dt does, requiring rtaus to be scaled by the
            // ratio of the new and old timesteps.
            rcount += 1;
            dto = dt/dto;

            for (size_t i=0; i<N; i++) {
                rtaus[i] *= dto;
            }
        } else {
            // Successful step
            return 0;
        }
    }
}
