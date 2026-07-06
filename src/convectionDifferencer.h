#pragma once

#include "mathUtils.h"
#include "grid.h"

//! Computes upwind-biased advective derivatives \f$ \partial y / \partial x \f$
//! on a nonuniform grid for the convection sub-problem.
//!
//! Two schemes are provided:
//!   - #Scheme::FirstOrderUpwind reproduces the legacy 1st-order upwind
//!     stencil bit-identically (forward difference where the advecting velocity
//!     is negative or at the left boundary node, backward difference
//!     otherwise).
//!   - #Scheme::SecondOrderLimited applies a van Albada limited, upwind-biased
//!     face reconstruction that is 2nd-order in smooth regions regardless of
//!     grid nonuniformity and falls back to local 1st-order upwind at extrema.
//!
//! The kernel computes derivatives at nodes `0 .. jj-1`. The right-boundary
//! node `jj` uses a boundary-condition-specific closure and is handled by the
//! callers.
class ConvectionDifferencer
{
public:
    //! Selects the discretization used by computeDerivatives().
    enum class Scheme {
        FirstOrderUpwind, //!< Legacy 1st-order upwind (bit-identical)
        SecondOrderLimited //!< van Albada limited ~2nd-order upwind-biased
    };

    //! Select the active discretization scheme.
    void setScheme(Scheme s);

    //! Size the internal limiter scratch vector. This is the intended
    //! pre-allocation entry point: call it with the current number of grid
    //! points before computeDerivatives() to avoid a reallocation on first
    //! use. computeDerivatives() also self-sizes the scratch vector
    //! defensively if it is ever found to be the wrong size (e.g. after a
    //! regrid where resize() was not called again), so calling it is not
    //! strictly required for correctness, only for avoiding that fallback.
    void resize(size_t nPoints);

    //! Compute the advective derivative \f$ dy/dx \f$ at nodes `0 .. jj-1`.
    //!
    //! The upwind branch is chosen per node from the sign of the advecting
    //! velocity `v[j]` (forward difference where `v[j] < 0`). The right-boundary
    //! node `jj` is not written; callers apply their own boundary closure there.
    //!
    //! @param y    Field values at each grid point.
    //! @param v    Advecting velocity at each grid point (only its sign is used
    //!             for branch selection).
    //! @param grid Grid providing the spacing arrays `hh` and `dlj`.
    //! @param dydx Output derivative, written for indices `0 .. jj-1`.
    void computeDerivatives(const dvec& y, const dvec& v,
                            const OneDimGrid& grid, dvec& dydx);

private:
    Scheme scheme = Scheme::SecondOrderLimited;
    dvec sigma; //!< limited node slopes (van Albada), one per grid point
};
