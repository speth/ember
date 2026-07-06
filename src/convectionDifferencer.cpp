#include "convectionDifferencer.h"

#include <cassert>

void ConvectionDifferencer::setScheme(Scheme s)
{
    scheme = s;
}

void ConvectionDifferencer::resize(size_t nPoints)
{
    sigma.resize(nPoints);
}

void ConvectionDifferencer::computeDerivatives(const dvec& y, const dvec& v,
                                               const OneDimGrid& grid, dvec& dydx)
{
    const size_t jj = grid.jj;
    const dvec& hh = grid.hh;

    if (scheme == Scheme::FirstOrderUpwind) {
        // Legacy 1st-order upwind, reproduced bit-identically from the original
        // convection loops: forward difference where the advecting velocity is
        // negative or at the left boundary node, backward difference otherwise.
        for (size_t j = 0; j < jj; j++) {
            if (v[j] < 0 || j == 0) {
                dydx[j] = (y[j+1] - y[j]) / hh[j];
            } else {
                dydx[j] = (y[j] - y[j-1]) / hh[j-1];
            }
        }
        return;
    }

    // *** SecondOrderLimited ***
    const size_t nPoints = grid.nPoints;
    assert(nPoints >= 3);
    if (static_cast<size_t>(sigma.size()) != nPoints) {
        sigma.resize(nPoints);
    }
    const dvec& dlj = grid.dlj;

    // Pass 1: van Albada limited node slopes. The outermost nodes fall back to
    // a zero slope, making the outermost face reconstructions locally 1st-order
    // (these are small-gradient far-field regions).
    sigma[0] = 0.0;
    sigma[jj] = 0.0;
    for (size_t j = 1; j < jj; j++) {
        double sm = (y[j] - y[j-1]) / hh[j-1]; // one-sided slope from the left
        double sp = (y[j+1] - y[j]) / hh[j];   // one-sided slope from the right
        if (sm * sp > 0) {
            // van Albada limiter. Returns exactly the common slope when
            // sm == sp (so linear data yields the exact derivative). The
            // sm*sp > 0 guard both clips at extrema and avoids the 0/0 form on
            // flat data (sm == sp == 0), without an epsilon that would break
            // the linear-exactness property.
            sigma[j] = (sm * sp) * (sm + sp) / (sm * sm + sp * sp);
        } else {
            sigma[j] = 0.0;
        }
    }

    // Pass 2: upwind-biased face reconstruction. The advective derivative at
    // node j is the difference of the reconstructed face values divided by the
    // node spacing dlj[j] = (hh[j-1] + hh[j])/2. Each face value is
    // reconstructed from the upwind side selected by sign(v[j]). dydx[0] is
    // unused by the callers (node 0 uses a boundary closure).
    dydx[0] = 0.0;
    for (size_t j = 1; j < jj; j++) {
        double yFaceR; // reconstruction at the right face x_{j+1/2}
        double yFaceL; // reconstruction at the left face  x_{j-1/2}
        if (v[j] < 0) {
            yFaceR = y[j+1] - sigma[j+1] * hh[j] / 2;
            yFaceL = y[j]   - sigma[j]   * hh[j-1] / 2;
        } else {
            yFaceR = y[j]   + sigma[j]   * hh[j] / 2;
            yFaceL = y[j-1] + sigma[j-1] * hh[j-1] / 2;
        }
        dydx[j] = (yFaceR - yFaceL) / dlj[j];
    }
}
