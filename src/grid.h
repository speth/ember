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

class OneDimGrid
{
public:
    OneDimGrid();

    dvec x; // The grid points
    dvec dampVal; // ratio of convective to diffusive coefficients (e.g. nu/v)

    int alpha; // domain curvature parameter. 0: planar, 1: cylindrical
    size_t ju, jb; // indices of burned / unburned boundaries

    // Parameters for controlling internal grid points:
    double vtol_in;
    double dvtol_in;
    dvec vtol; // relative solution variable tolerance for point insertion
    dvec dvtol; // global derivative solution variable tolerance for point insertion
    double absvtol; // absolute tolerance (ignore components with range smaller than this)
    double rmTol; // relative grid point removal tolerance
    double uniformityTol; // maximum ratio of adjacent grid point separation distances
    double gridMin, gridMax; // minimum and maximum grid point separation
    double dampConst; // relative allowable numerical damping
    double centerGridMin; // minimum separation of center grid points for curved flames

    // Parameters for controlling exterior grid points:
    bool fixedBurnedVal;
    bool unburnedLeft;
    bool fixedLeftLoc;
    bool twinFlame;
    bool curvedFlame;

    double boundaryTol;
    double boundaryTolRm;
    double unstrainedDownstreamWidth;
    size_t addPointCount; // number of points to add when regridding


    // Derived mesh parameters (calculated by updateValues)
    dvec cfm, cf, cfp; // first centered difference
    dvec dlj;
    dvec hh;
    dvec rphalf;
    dvec r;
    size_t nVars; // number of variables at each grid point
    size_t nAdapt; // Only the first nAdapt variables at each point are used for adaptation
    size_t nPoints; // number of grid point
    size_t jj; // index of last grid point ( = nPoints-1)

    bool updated; // true if the grid changed on the last call to adapt or regrid

    // Update the grid based on the solutionState, and adjust it to fit
    // on the new grid. Each element of the solutionState is a vector
    // containing the value of a solution variable at each grid point.
    // Return value is true if the grid has been modified
    // adapt inserts and removes points within the problem domain
    // regrid inserts and removes points at the ends of the problem domain
    void adapt(vector<dvector>& y);
    void regrid(vector<dvector>& y);
    void regridUnstrained(vector<dvector>& y, dvec& qdot);
    void setOptions(const ConfigOptions& options);
    void updateValues(void);
    void setSize(const size_t N);

    void updateBoundaryIndices(void);

    BoundaryCondition::BC leftBC;
    BoundaryCondition::BC rightBC;

private:
    void removePoint(int jRemove, vector<dvector>& y);
    void addPoint(int jInsert, vector<dvector>& y);

    bool removeRight(vector<dvector>& y);
    bool removeLeft(vector<dvector>& y);
    bool addRight(vector<dvector>& y);
    bool addLeft(vector<dvector>& y);
    bool addRightUnstrained(vector<dvector>& y, dvec& q);
    bool removeRightUnstrained(vector<dvector>& y, dvec& q);
};


class GridBased
{
public:
    GridBased();
    virtual ~GridBased() {}

    virtual void setGrid(const OneDimGrid& grid);

    // the grid:
    OneDimGrid grid;

protected:

    // local names for some things that are part of the grid:
    dvec& x;
    dvec& r;
    dvec& rphalf;
    dvec& hh;
    dvec& dlj;
    dvec& cfm;
    dvec& cf;
    dvec& cfp;
    int& alpha; //!< curved grid exponent. alpha = 1 for curved flames, 0 for planar flames.

    size_t& nPoints;
    size_t& jj;
};
