#pragma once

#include "readConfig.h"

//! Boundary condition options for each component at the right and left
//! boundaries of the domain.
namespace BoundaryCondition {
    enum BC {
     //! Fix the value at the boundary by setting dy/dt = 0.
     FixedValue,

     //! Zero gradient condition by setting y[j] = y[j-1]
     ZeroGradient,

     //! Treat the boundary value as representing the average value within a
     //! control volume extending from x[0] to x[1].
     ControlVolume,

     //! Flux of y at the boundary driven by difference between yLeft and y[0]
     WallFlux,

     //! Outflow boundary condition
     Floating
    };
}

//! Representation of an adaptive, non-uniform one-dimensional grid.
class OneDimGrid
{
public:
    OneDimGrid();

    dvec x; //!< The coordinates of the grid points [m].

    //! Ratio of convective to diffusive coefficients (e.g.\ nu/v). Used to set
    //! an upper bound on grid spacing to control numerical diffusion.
    dvec dampVal;

    int alpha; //!< Domain curvature parameter. 0: planar / disc, 1: cylindrical
    int beta; //!< Domain curvature parameter. 1: planar / cylindrical; 2: disc
    size_t ju; //!< Index corresponding to the unburned mixture
    size_t jb; //!< Index corresponding to the burned mixture

    // **** Parameters for controlling internal grid points ****

    //! Absolute tolerance for point insertion. Components with ranges smaller
    //! than this value are not considered during grid adaptation.
    double absvtol;

    //! Local-error tolerance for grid adaptation. The estimated local
    //! representation error of each adapted component must be smaller than
    //! errTol times the component's range. Set from ConfigOptions in
    //! setOptions().
    double errTol;

    //! Spatial order p of the active convection scheme (1: firstOrderUpwind,
    //! 2: secondOrderLimited). The adaptation error estimate scales as
    //! h^(p+1) * |d^(p+1)v/dx^(p+1)|. Set in setOptions().
    int errorOrder;

    //! Leading coefficient C_p of the local error estimate
    //! E = C_p * h^(p+1) * |d^(p+1)v|. Initial values are the classical
    //! interpolation-error bounds (1/8 for p=1, 1/15 for p=2); the p=2
    //! value is calibrated so matched errTol gives matched accuracy across
    //! convection schemes. Set in setOptions().
    double errCoeff;

    //! Relative grid point removal tolerance. Grid points can be removed if
    //! all criteria are satisfied to this multiplier on the insertion
    //! tolerances.
    double rmTol;

    //! Maximum ratio of adjacent grid point separation distances.
    double uniformityTol;

    double gridMin; //!< Minimum grid point separation
    double gridMax; //!< Maximum grid point separation

    //! Relative allowable numerical damping. Compared with #dampVal.
    double dampConst;

    //!< Minimum separation of center grid points for curved flames.
    double centerGridMin;

    // **** Parameters for controlling exterior grid points ****

    //! `true` if the composition and temperature of the burned gas should be
    //! held constant.
    bool fixedBurnedVal;

    //! `true` if the unburned mixture is on the left side of the domain.
    bool unburnedLeft;

    //! `true` if the coordinate of the leftmost grid point is fixed. Prevents
    //! addition or removal of grid points at the left boundary.
    bool fixedLeftLoc;

    //! `true` if the left boundary condition is a symmetry condition.
    bool twinFlame;

    //! `true` if the flame is curved (corresponding to #alpha = 1)
    bool cylindricalFlame;

    //! `true` if the flame is opposed axisymmetric jets (corresponding to #beta = 2)
    bool discFlame;

    //! Relative tolerance used to extend the domain in order to satisfy zero-
    //! gradient boundary conditions.
    double boundaryTol;

    //! Relative tolerance used to remove boundary grid points satisfying
    //! zero-gradient boundary conditions.
    double boundaryTolRm;

    //! Maximum extent of the domain downstream of the flame for unstrained flames.
    //! Expressed as a multiple of the flame thickness.
    double unstrainedDownstreamWidth;

    size_t addPointCount; //!< number of points to add when regridding

    // **** Derived mesh parameters (calculated by updateValues()) ****
    dvec cfm; //!< Coefficient for y[j-1] in first centered difference
    dvec cf; //!< Coefficient for y[j] in first centered difference
    dvec cfp; //!< Coefficient for y[j+1] in first centered difference
    dvec dlj; //!< Average of left and right grid spacing
    dvec hh; //!< Grid spacing between x[j] and x[j+1]
    dvec rphalf; //!< "radius" at x[j+1/2]. Equal to 1 for planar geometries.
    dvec r; //!< "radius" at x[j]. Equal to 1 for planar geometries.
    size_t nVars; //!< number of variables at each grid point

    //! Only the first nAdapt variables at each point are used for adaptation
    size_t nAdapt;
    size_t nPoints; //!< number of grid point
    size_t jj; //!< index of last grid point (`== nPoints-1`)

    //! `true` if the grid changed on the last call to adapt() or regrid()
    bool updated;

    //! Update the grid based on the solution state `y` to satisfy all
    //! refinement criteria for internal grid points. Also adjusts the solution
    //! to fit on the new grid. The first #nAdapt elements of `y` are vector
    //! containing the values of a state variable at each grid point.
    //! Subsequent elements are arrays (such as the time derivative of a state
    //! variable) which needs to be updated to fit on the new grid but are not
    //! required to meet adaptation tolerances.
    //!
    //! This function tries to remove unnecessary grid points located in the
    //! regions of small gradients, insert new grid points in the regions of
    //! large gradients, and, at the same time maintain relative uniformity of
    //! the grid.
    //!
    //! *Adaptation algorithm*
    //!
    //! Insertion is performed first. For each adapted component f with
    //! range(f) >= #absvtol, the local representation error of interval j
    //! is estimated as
    //!     E = C_p * hh[j]^(p+1) * max|d^(p+1) f|
    //! where p (#errorOrder) is the spatial order of the convection scheme,
    //! C_p is #errCoeff, and the derivative magnitude is taken as the
    //! maximum of the estimates (computeErrorWeights()) at the nodes within
    //! one interval of j. A point is inserted in interval j if
    //!     E > errTol * range(f)
    //! or if the damping, maximum-spacing, or uniformity criteria require
    //! one.
    //!
    //! Removal is considered next: point j is removed only if, for every
    //! component, the estimate E evaluated for the merged interval
    //! (hh[j]+hh[j-1], derivative maximum over nodes within two of j) stays
    //! below rmTol * errTol * range(f), and the damping, maximum-spacing,
    //! and uniformity criteria all permit it. Since merging roughly doubles
    //! h and E scales as h^(p+1), removal has strong built-in hysteresis.
    void adapt(vector<dvector>& y);

    //! Add and remove points at the boundaries of the domain to satisfy
    //! boundary tolerances. The structure of `y` is the same as for adapt().
    void regrid(vector<dvector>& y);

    //! Add and remove points at the boundaries of the domain to satisfy
    //! boundary tolerances for unstrained flames. The structure of `y` is the
    //! same as for adapt().
    void regridUnstrained(vector<dvector>& y, dvec& qdot);

    //! Set adaptation and regridding tolerances specified in a ConfigOptions
    //! object.
    void setOptions(const ConfigOptions& options);

    //! Estimate the nodal magnitude of the (#errorOrder + 1)-th derivative
    //! of v, the weight in the local error estimate used by adapt().
    //! Repeated application of the nonuniform first-derivative stencils;
    //! end nodes copy the nearest interior estimate. updateValues() must
    //! have been called first.
    void computeErrorWeights(const dvector& v, dvector& W) const;

    //! Recompute derived mesh parameters after the grid has changed.
    void updateValues(void);

    //! Set the size of the grid after making changes to the grid.
    void setSize(const size_t N);

    //! Update values of #ju and #jb after regridding.
    void updateBoundaryIndices(void);

    BoundaryCondition::BC leftBC; //!< Boundary condition applied at `j = 0`
    BoundaryCondition::BC rightBC; //!< Boundary condition applied at `j = jj`

private:
    //! Remove the point `jRemove` and update the solution vector `y`.
    void removePoint(int jRemove, vector<dvector>& y);

    //! Add a point to the right of `jInsert` and update the solution vector
    //! `y`.
    void addPoint(int jInsert, vector<dvector>& y);

    //! Remove any unnecessary points from the right side of the domain and
    //! update the solution vector `y`.
    bool removeRight(vector<dvector>& y);

    //! Remove any unnecessary points from the left side of the domain and
    //! update the solution vector `y`.
    bool removeLeft(vector<dvector>& y);

    //! Add points to the right side of the domain if necessary and update the
    //! solution vector `y`.
    bool addRight(vector<dvector>& y);

    //! Add points to the left side of the domain if necessary and update the
    //! solution vector `y`.
    bool addLeft(vector<dvector>& y);

    //! Add points to the right side of an unstrained flame if necessary and
    //! update the solution vector `y`. `q` is the local heat release rate.
    bool addRightUnstrained(vector<dvector>& y, dvec& q);

    //! Remove unnecessary points from the right side of an unstrained flame
    //! and update the solution vector `y`. `q` is the local heat release
    //! rate.
    bool removeRightUnstrained(vector<dvector>& y, dvec& q);
};

//! A "mix-in" class used by classes that need frequent access to grid
//! parameters.
class GridBased
{
public:
    GridBased();
    virtual ~GridBased() {}

    //! Copy the specified grid to this object
    virtual void setGrid(const OneDimGrid& grid);

    OneDimGrid grid; //!< the actual grid

protected:
    // local names for some things that are part of the grid
    dvec& x; //!< The coordinates of the grid points [m].
    dvec& r; //!< "radius" at x[j]. Equal to 1 for planar geometries.
    dvec& rphalf; //!< "radius" at x[j+1/2]. Equal to 1 for planar geometries.
    dvec& hh; //!< Grid spacing between x[j] and x[j+1]
    dvec& dlj; //!< Average of left and right grid spacing
    dvec& cfm; //!< Coefficient for y[j-1] in first centered difference
    dvec& cf; //!< Coefficient for y[j] in first centered difference
    dvec& cfp; //!< Coefficient for y[j+1] in first centered difference
    int& alpha; //!< curved grid exponent. alpha = 1 for curved flames, 0 for planar flames and axisymmetric jet flames.
    int& beta; //!< curved grid exponent. beta = 2 for axisymmetric jet flames, 1 for planar flames and curved flames.

    size_t& nPoints; //!< number of grid point
    size_t& jj; //!< index of last grid point (`== nPoints-1`)
};
