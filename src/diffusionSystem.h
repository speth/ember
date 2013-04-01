#pragma once

#include "mathUtils.h"
#include "grid.h"
#include "integrator.h"

//! System representing diffusion of a single solution component.
//! Represents an ODE in one of the following forms:
//!     \f[ \dot{y} = B \frac{d}{dx}\left(D \frac{dy}{dx}\right) + C \f]
//!     \f[ \dot{y} = \frac{B}{r} \frac{d}{dr}\left(r D \frac{dy}{dr}\right) + C \f]
//! The ODE in this form may be written as a linear system:
//!     \f[ \dot{y} = Ay + k \f]
//! where the entries of the matrix *A* are determined by the prefactor *B*,
//! the diffusion coefficients *D*, the boundary conditions and the finite
//! difference formulas used.
class DiffusionSystem : public TridiagonalODE, public GridBased
{
public:
    DiffusionSystem();

    //! Build the matrix *A* describing the linear %ODE. `a[j]`, `b[j]`, and
    //! `c[j]` are respectivley the subdiagonal, diagonal, and superdiagonal
    //! elements of row `j`.
    void get_A(dvec& a, dvec& b, dvec& c);

    void get_k(dvec& k);
    void resize(size_t N);

    //! Set the splitting constant *C* to zero.
    void resetSplitConstants();

    dvec B; //!< scaling factor for the spatial derivative
    dvec D; //!< diffusion coefficient or equivalent
    dvec splitConst; //!< Balancing constant introduced by the splitting method
    size_t N; //!< Number of points in the region to be solved

    double yInf; //!< reference value for wall flux boundary condition
    double wallConst; //!< wall flux proportionality constant

private:
    // constants used to compute the matrix coefficients
    dvec c1;
    dvec c2;
};
