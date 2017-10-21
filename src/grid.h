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

    //! Default relative solution variable tolerance for point insertion.
    //! Set from ConfigOptions in setOptions().
    double vtol_in;

    //! Default derivative tolerance for point insertion.
    //! Set from ConfigOptions in setOptions().
    double dvtol_in;

    //! Relative solution variable tolerance for point insertion for each
    //! solution component. Length #nVars.
    dvec vtol;

    //! solution variable derivative tolerance for point insertion for each
    //! solution component. Length #nVars.
    dvec dvtol;

    //! Absolute tolerance for point insertion. Components with ranges smaller
    //! than this value are not considered during grid adaptation.
    double absvtol;

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
    //! The insertion of the grid points is performed first. For each
    //! component of the solution vector:
    //!
    //! 1. Find its range and the range of its derivative.
    //! 2. Apply four criteria and and find where insertions are
    //!    needed. the criteria for a component f(j) are:
    //!      - `|f[j+1]-f[j]| < vtol*range(f)`
    //!      - `|dfdy[j+1]-dfdy[j]| < dvtol*range(dfdy)`
    //!      - `1/uniformityTol < hh[j]/hh[j-1] < uniformityTol`
    //! 3. If any of these criteria is not satisfied, a grid point j
    //!    is inserted.
    //!
    //! Next, the unnecessary grid points are removed, and the algorithm is
    //! applied in reverse. If the criteria:
    //!   - `|f[j]-f[j-1]| > rmTol*vtol*range(f)`
    //!   - `|dfdy[j]-dfdy[j-1]| > rmTol*dvtol*range(dfdy)`
    //!   - `hh[j]+hh[j-1] < uniformityTol*hh[j-2]`
    //!   - `hh[j]+hh[j-1] < uniformityTol*hh[j+1]`
    //!
    //! are satisfied for all components at a point, it is removed.
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
