#include "grid.h"
#include "debugUtils.h"

#include <assert.h>

using namespace mathUtils;
using std::max; using std::min;

OneDimGrid::OneDimGrid()
    : updated(true)
    , leftBC(BoundaryCondition::FixedValue)
    , rightBC(BoundaryCondition::FixedValue)
{
}

void OneDimGrid::setOptions(const ConfigOptions& options)
{
    vtol_in = options.vtol;
    dvtol_in = options.dvtol;
    absvtol = options.absvtol;
    rmTol = options.rmTol;
    uniformityTol = options.uniformityTol;
    gridMin = options.gridMin;
    gridMax = options.gridMax;
    dampConst = options.dampConst;
    centerGridMin = options.centerGridMin;

    fixedBurnedVal = options.fixedBurnedVal;
    unburnedLeft = options.unburnedLeft;
    fixedLeftLoc = options.fixedLeftLoc;
    twinFlame = options.twinFlame;
    cylindricalFlame = options.cylindricalFlame;
    discFlame = options.discFlame;

    boundaryTol = options.boundaryTol;
    boundaryTolRm = options.boundaryTolRm;
    unstrainedDownstreamWidth = options.unstrainedDownstreamWidth;
    addPointCount = options.addPointCount;
    alpha = (options.cylindricalFlame) ? 1 : 0;
    beta = (options.discFlame) ? 2 : 1;
}

void OneDimGrid::updateValues()
{
    // Derived mesh sizes
    hh.resize(jj);
    cfm.resize(jj);
    cf.resize(jj);
    cfp.resize(jj);

    dlj.resize(jj);
    rphalf.resize(jj);
    r.resize(jj+1);
    dampVal.resize(jj+1);

    for (index_t j=0; j<x.rows()-1; j++) {
        hh[j] = x[j+1]-x[j];
        rphalf[j] =  pow(0.5*(x[j]+x[j+1]),alpha);
    }
    r = x.pow(alpha);

    for (index_t j=1; j<x.rows()-1; j++) {
        cfp[j] = hh[j-1]/(hh[j]*(hh[j]+hh[j-1]));
        cf[j] = (hh[j]-hh[j-1])/(hh[j]*hh[j-1]);
        cfm[j] = -hh[j]/(hh[j-1]*(hh[j]+hh[j-1]));

        dlj[j]= 0.5 * (x[j+1]-x[j-1]);
    }
}

void OneDimGrid::setSize(const size_t N)
{
    nPoints = N;
    jj = N-1;
}

void OneDimGrid::adapt(vector<dvector>& y)
{
    nVars = y.size();
    assert(nAdapt <= nVars);
    assert((dampVal > 0).all());
    setSize(y[0].size());

    vtol.resize(nAdapt);
    dvtol.resize(nAdapt);
    for (size_t k=0; k<nAdapt; k++) {
        vtol[k] = vtol_in;
        dvtol[k] = dvtol_in;
    }

    // Used for informational purposes only
    std::vector<int> insertionIndicies;
    std::vector<int> removalIndices;

    // *** Grid point insertion algorithm

    size_t j = 0;
    dvector dv(jj+1); // dv/dx

    while (j < jj) {
        updateValues();
        dv.resize(jj+1);
        bool insert = false;

        // Consider tolerances for each variable v in the solution y
        for (size_t k=0; k<nAdapt; k++) {

            dvector& v = y[k];
            for (size_t i=1; i<jj; i++) {
                dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
            }

            double vRange = mathUtils::range(v);
            double dvRange = mathUtils::range(dv,1,jj-1);

            if (vRange < absvtol) {
                continue; // Ignore minor species
            }

            // Apply grid point addition criteria:

            // resolution of v
            if (abs(v[j+1]-v[j]) > vtol[k]*vRange) {
                insert = true;
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: v resolution wants grid point"
                        " j = %i, k = %i; |v(j+1)-v(j)|/vrange = %g > %g") %
                        j % k % (abs(v[j+1]-v[j])/vRange) % vtol[k]);
                }
            }

            // resolution of dv
            if (j!=0 && j!=jj-1 && abs(dv[j+1]-dv[j]) > dvtol[k]*dvRange) {
                insert = true;
                if (debugParameters::debugAdapt) {
                    logFile.write(format(
                        "Adapt: dv resolution (global) wants grid point"
                        " j = %i, k = %i; |dv(j+1)-dv(j)|/vrange = %g > % g") %
                        j % k % (abs(dv[j+1]-dv[j])/dvRange) % dvtol[k]);
                }
            }
        }

        // Damping of high-frequency numerical error
        if (hh[j] > dampConst*dampVal[j]) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format(
                    "Adapt: damping criterion wants a grid point"
                    " j = %i; hh[j] = %g > %g") %
                    j % hh[j] % (dampConst*dampVal[j]));
            }
        }

        // Maximum grid size
        if (hh[j] > gridMax) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format(
                    "Adapt: Maximum grid size criterion wants a grid point"
                    " j = %i; hh[j] = %g < %g") %
                    j % hh[j] % gridMax);
            }
        }

        // Left uniformity
        if (j!=0 && hh[j]/hh[j-1] > uniformityTol) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format(
                    "Adapt: left uniformity wants grid point"
                    " j = %i; hh(j)/hh(j-1) = %g > %g") %
                    j % (hh[j]/hh[j-1]) % uniformityTol);
            }
        }

        // Right uniformity
        if (j!=jj-1 && hh[j]/hh[j+1] > uniformityTol) {
            insert = true;
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: right uniformity wants grid point"
                        " j = %i; hh(j)/hh(j+1) = %g > %g") %
                        j % (hh[j]/hh[j+1]) % uniformityTol);
            }
        }

        // Special minimum grid size for flames pinned at x=0
        if (j == 0 && (leftBC == BoundaryCondition::ControlVolume ||
                       leftBC == BoundaryCondition::WallFlux))
        {
            double xLeftMin = min(centerGridMin, 0.02*x[jj]);
            if (hh[j] < 2*xLeftMin) {
                insert = false;
                if (debugParameters::debugAdapt) {
                    logFile.write(format("Adapt: grid point addition canceled"
                        " by minimum center grid size j = %i; hh[j] = %g < %g") %
                        j % hh[j] % (2*xLeftMin));
                }
            }
        }

        // Enforce minimum grid size
        if (insert && hh[j] < 2*gridMin) {
            insert = false;
            if (debugParameters::debugAdapt) {
                logFile.write(format("Adapt: grid point addition canceled"
                    " by minimum grid size j = %i; hh[j] = %g < %g") %
                    j % hh[j] % (2*gridMin));
            }
        }

        if (insert) {
            // Insert a new point
            insertionIndicies.push_back(static_cast<int>(j));
            addPoint(static_cast<int>(j), y);
            updated = true;
            setSize(nPoints+1);
            j+=2;
        } else {
            // No insertion; step to the next point.
            j++;
        }
    }

    if (insertionIndicies.size() != 0) {
        logFile.write("Grid points inserted at j = ", false);
        for (unsigned int i=0; i<insertionIndicies.size()-1; i++) {
            logFile.write(format("%i, ") % insertionIndicies[i], false);
        }
        logFile.write(insertionIndicies[insertionIndicies.size()-1]);
    }

    // *** Grid point removal algorithm

    j = 1;
    while (j < jj) {
        updateValues();

        // Assume removal, then look for a condition which prevents removal
        bool remove = true;

        // Consider tolerances each variable v in the solution y
        for (size_t k=0; k<nAdapt; k++) {

            dvector& v = y[k];
            for (size_t i=1; i<jj; i++) {
                dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
            }

            double vRange = mathUtils::range(v);
            double dvRange = mathUtils::range(dv,1,jj-1);

            if (vRange < absvtol) {
                continue; // Ignore minor species
            }

            // Apply grid point removal criteria:

            // resolution of v
            if (abs(v[j+1]-v[j-1]) > rmTol*vtol[k]*vRange) {
                if (debugParameters::debugAdapt) {
                    logFile.write(format(
                        "Adapt: no removal - v res. j = %i, k = %i;"
                        " |v[j+1]-v[j-1]|/vtrange = %g > %g") %
                        j % k % (abs(v[j+1]-v[j-1])/vRange) % (vtol[k]*rmTol));
                }
                remove = false;
            }

            // resolution of dv
            if (j!=2 && j!=jj-1 && abs(dv[j+1]-dv[j-1]) > rmTol*dvtol[k]*dvRange) {
                if (debugParameters::debugAdapt) {
                    logFile.write(format(
                        "Adapt: no removal - dv res. j = %i, k = %i;"
                        " |dv(j+1)-dv(j-1)|/dvrange = %g > %g") %
                        j % k % (abs(dv[j+1]-dv[j-1])/dvRange) % (dvtol[k]*rmTol));
                }
                remove = false;
            }
        }

        // Damping of high-frequency numerical error
        if (hh[j]+hh[j-1] >= rmTol*dampConst*dampVal[j]) {
            if (debugParameters::debugAdapt) {
                logFile.write(format(
                    "Adapt: no removal - damping criterion. j = %i;"
                    " hh(j)+hh(j-1) = %g > %g") %
                    j % (hh[j]+hh[j-1]) % (dampConst*dampVal[j]));
            }
            remove = false;
        }

        // Enforce maximum grid size
        if (hh[j]+hh[j-1] > gridMax) {
            if (debugParameters::debugAdapt) {
                logFile.write(format(
                    "Adapt: no removal - maximum grid size. j = %i;"
                    " hh(j)+hh(j-1) = %g > %g") %
                    j % (hh[j]+hh[j-1]) % gridMax);
            }
            remove = false;
        }

        // Enforce left uniformity
        if (j>=2 && hh[j]+hh[j-1] > uniformityTol*hh[j-2]) {
            if (debugParameters::debugAdapt) {
                logFile.write(format(
                    "Adapt: no removal - left uniformity. j = %i;"
                    " (hh(j)+hh(j-1))/hh(j-2) = %g > %g") %
                    j % ((hh[j]+hh[j-1])/hh[j-2]) % uniformityTol);
            }
            remove = false;
        }

        // Enforce right uniformity
        if (j<=jj-2 && hh[j]+hh[j-1] > uniformityTol*hh[j+1]) {
            if (debugParameters::debugAdapt) {
                logFile.write(format(
                    "Adapt: no removal - right uniformity. j = %i;"
                    " (hh(j)+hh(j-1))/hh(j+1) = %g > %g") %
                    j % ((hh[j]+hh[j-1])/hh[j+1]) % uniformityTol);
            }
            remove = false;
        }

        // Special fixed grid for flames pinned at x=0
        if (j == 1 && (leftBC == BoundaryCondition::ControlVolume ||
                       leftBC == BoundaryCondition::WallFlux))
        {
            if (debugParameters::debugAdapt) {
                logFile.write("Adapt: no removal - fixed grid near r = 0.");
            }
            remove = false;
        }

        if (remove) {
            removalIndices.push_back(static_cast<int>(j));
            removePoint(static_cast<int>(j), y);
            setSize(nPoints-1);
            updated = true;
        } else {
            j++;
        }
    }

    if (removalIndices.size() != 0) {
        logFile.write("Grid points removed at j = ", false);
        for (unsigned int i=0; i<removalIndices.size()-1; i++) {
            logFile.write(format("%i, ") % removalIndices[i], false);
        }
        logFile.write(removalIndices[removalIndices.size()-1]);
    }

    if (updated) {
        updateValues();
        updateBoundaryIndices();
    }
}

void OneDimGrid::addPoint(int jInsert, vector<dvector>& y)
{
    assert(x.rows() == dampVal.rows());
    dvec::Index N = x.rows();

    double xInsert = 0.5*(x[jInsert+1]+x[jInsert]);

    double val = mathUtils::splines(x, dampVal, xInsert);
    dvec tmp(N + 1);
    tmp << dampVal.head(jInsert + 1), val, dampVal.tail(N - jInsert - 1);
    dampVal = tmp;

    for (dvector& row : y) {
        double yNew = mathUtils::splines(x, Eigen::Map<dvec>(&row[0], row.size()),
                                         xInsert);
        row.insert(row.begin()+jInsert+1, yNew);
    }

    tmp << x.head(jInsert + 1), xInsert, x.tail(N - jInsert - 1);
    x = tmp;
}

void OneDimGrid::removePoint(int jRemove, vector<dvector>& y)
{
    dvec tmp(x.rows() - 1);
    tmp << x.head(jRemove), x.tail(x.rows() - jRemove - 1);
    x = tmp;
    tmp << dampVal.head(jRemove), dampVal.tail(dampVal.rows() - jRemove - 1);
    dampVal = tmp;
    for (dvector& row : y) {
        row.erase(row.begin() + jRemove);
    }
}

bool OneDimGrid::addRight(vector<dvector>& y)
{
    // *** Criteria for addition to right (j==jj) ***

    // Pick the comparison point for the grid flatness criterion,
    // depending on whether a zero-gradient condition is being enforced

    // djMom is for momentum equation, which always uses a special
    // zero-gradient condition.
    size_t djMom = 1;

    // All other variables use fixed or zero gradient conditions
    // depending on fixedBurnedVal
    size_t djOther = (jb==jj && !fixedBurnedVal) ? 2 : 1;

    bool pointAdded = false; // Assume no addition

    // Check flatness of temperature, velocity and species
    // profiles at the boundary.
    for (size_t k=0; k<nAdapt; k++) {
        size_t dj = (k == kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][jj]-y[k][jj-dj])/ymax > boundaryTol && ymax > absvtol ) {
            if (!pointAdded && debugParameters::debugRegrid) {
                logFile.write(format(
                    "Regrid: Flatness of component %i requires"
                    " right addition: %g > %g") %
                    k % (abs(y[k][jj]-y[k][jj-dj])/ymax) % boundaryTol);
            }
            pointAdded = true;
            break;
        }
    }

    if (!pointAdded && debugParameters::debugRegrid) {
        logFile.write("Regrid: no addition on right");
    }

    if (pointAdded) {
        logFile.write("Regrid: Adding points to right side.");

        x.conservativeResize(nPoints + addPointCount);
        dampVal.conservativeResize(nPoints + addPointCount);
        for (size_t i=0; i<addPointCount; i++) {
            x[jj+1] = x[jj] + pow(uniformityTol, 1.0/(1+addPointCount)) * (x[jj]-x[jj-1]);
            if (debugParameters::debugRegrid) {
                logFile.write(format("Regrid: point added at %g (hh = %g)") %
                    x[jj+1] % (x[jj+1]-x[jj]));
            }
            dampVal[jj+1] = dampVal[jj];
            for (size_t k=0; k<nVars; k++) {
                // keep constant boundary value
                y[k].push_back(y[k][jj]);
            }
            setSize(nPoints+1);
        }
        updateBoundaryIndices();
    }

    return pointAdded;
}


bool OneDimGrid::addRightUnstrained(vector<dvector>& y, dvec& qdot)
{
    // *** Criteria for addition to right (j==jj) ***

    // Pick the point to check for the heat release rate criterion,
    // depending on whether a zero-gradient condition is being enforced
    bool pointAdded = false; // Assume no addition

    dvec::Index jqMax;
    double qMax = qdot.maxCoeff(&jqMax);

    size_t jLeft = 0;
    size_t jRight = jj;
    for (size_t j=jqMax; j>=1; j--) {
        if (std::abs(qdot[j]) < 0.01 * qMax) {
            jLeft = j;
            break;
        }
    }

    for (size_t j=jqMax; j<=jj; j++) {
        if (std::abs(qdot[j]) < 0.01 * qMax) {
            jRight = j;
            break;
        }
    }

    double reactionWidth = x[jRight] - x[jLeft];

    if (x[jqMax] + unstrainedDownstreamWidth * reactionWidth > x[jj]) {
        pointAdded = true;
    }

    if (!pointAdded && debugParameters::debugRegrid) {
        logFile.write("Regrid: no addition on right");
    }

    if (pointAdded) {
        logFile.write("Regrid: Adding points to right side.");

        x.conservativeResize(nPoints + addPointCount);
        qdot.conservativeResize(nPoints + addPointCount);
        dampVal.conservativeResize(nPoints + addPointCount);
        for (size_t i=0; i<addPointCount; i++) {
            x[jj+1] = x[jj] + pow(uniformityTol, 1.0/(1+addPointCount)) * (x[jj]-x[jj-1]);
            qdot[jj+1] = qdot[jj];
            dampVal[jj+1] = dampVal[jj];
            for (size_t k=0; k<nVars; k++) {
                // keep constant boundary value
                y[k].push_back(y[k][jj]);
            }
            setSize(nPoints+1);
        }
        updateBoundaryIndices();
    }

    return pointAdded;
}


bool OneDimGrid::addLeft(vector<dvector>& y)
{
    // *** Criteria for addition to the left (j==0) ***

    int djOther = (jb==1 && !fixedBurnedVal) ? 2 : 1;
    int djMom = 1;

    bool pointAdded = false;
    if (!fixedLeftLoc) {
        for (size_t k=0; k<nAdapt; k++) {
            double ymax = maxval(y[k]);
            int dj = (k==kMomentum) ? djMom : djOther;
            if (abs(y[k][dj]-y[k][0])/ymax > boundaryTol && ymax > absvtol) {
                if (!pointAdded && debugParameters::debugRegrid) {
                    logFile.write(format("Regrid: Flatness of component %i"
                            " requires left addition: %g > %g") %
                            k % (abs(y[k][dj]-y[k][0])/ymax) % boundaryTol);
                }
                pointAdded = true;
                break;
            }
        }
    }

    if (fixedLeftLoc &&
        leftBC != BoundaryCondition::ControlVolume &&
        leftBC != BoundaryCondition::WallFlux &&
        x[0] > 0.0)
    {
        if (!pointAdded && debugParameters::debugRegrid) {
            logFile.write(
                "Regrid: Adding point to force left boundary toward x = 0");
        }
        pointAdded = true;
    }

    if (!pointAdded && debugParameters::debugRegrid) {
        logFile.write("Regrid: No addition to left.");
    }

    if (pointAdded) {
        logFile.write("Regrid: Adding points to left side.");
        // Add point to the left.
        for (size_t i=0; i<addPointCount; i++) {
            double xLeft = x[0] - sqrt(uniformityTol) * (x[1]-x[0]);
            if (twinFlame || cylindricalFlame) {
                if (x[0] == 0) {
                    break;
                }

                double xLeftMin = std::min(centerGridMin, 0.005*x[jj]);
                if (xLeft < 0.0) {
                    xLeft = 0.0;
                } else if (xLeft < xLeftMin) {
                    xLeft = xLeftMin;
                }
            }
            if (debugParameters::debugRegrid) {
                logFile.write(format("Regrid: adding point at %i (hh = %g)") %
                    xLeft % (x[0]-xLeft));
            }
            dvec tmp(x.rows() + 1);
            tmp << xLeft, x;
            x = tmp;
            tmp << dampVal[0], dampVal;
            dampVal = tmp;

            for (size_t k=0; k<nVars; k++) {
                // keep constant boundary value
                y[k].insert(y[k].begin(),y[k][0]);
            }
            setSize(nPoints+1);
        }
        updateBoundaryIndices();
    }

    return pointAdded;
}

bool OneDimGrid::removeRight(vector<dvector>& y)
{
    // *** Criteria for removal from the right (j==jj) ***

    // Comparison point for flatness criteria, depending on
    // zero-gradient or fixed value boundary condition
    size_t djMom = 2;
    size_t djOther = (jb==jj && !fixedBurnedVal) ? 3 : 2;

    bool pointRemoved = true; // assume removal
    for (size_t k=0; k<nAdapt; k++) {
        size_t dj = (k==kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][jj]-y[k][jj-dj])/ymax > boundaryTolRm && ymax > absvtol) {
            if (pointRemoved && debugParameters::debugRegrid) {
                logFile.write(format(
                    "Regrid: Right removal prevented by component %i"
                    " flatness: %g > %g") %
                    k % (abs(y[k][jj]-y[k][jj-dj])/ymax) % boundaryTolRm);
            }
            pointRemoved = false;
            break;
        }
    }

    if (jj < 3) {
        pointRemoved = false;
    }

    if (pointRemoved) {
        removePoint(static_cast<int>(jj), y);
        setSize(nPoints-1);
        updateBoundaryIndices();
    }

    return pointRemoved;
}

bool OneDimGrid::removeRightUnstrained(vector<dvector>& y, dvec& qdot)
{
    // *** Criterion for removal from the right (j==jj) ***

    bool pointRemoved = false;
    dvec::Index jqMax;
    double qMax = qdot.maxCoeff(&jqMax);

    size_t jLeft = 0;
    size_t jRight = jj;
    for (size_t j=jqMax; j>=1; j--) {
        if (std::abs(qdot[j]) < 0.01 * qMax) {
            jLeft = j;
            break;
        }
    }

    for (size_t j=jqMax; j<=jj; j++) {
        if (std::abs(qdot[j]) < 0.01 * qMax) {
            jRight = j;
            break;
        }
    }

    double reactionWidth = x[jRight] - x[jLeft];

    if (x[jj] > x[jqMax] + 1.25 * unstrainedDownstreamWidth * reactionWidth) {
        pointRemoved = true;
    }

    if (jj < 3) {
        pointRemoved = false;
    }

    if (pointRemoved) {
        removePoint(static_cast<int>(jj), y);
        qdot.conservativeResize(nPoints-1);
        setSize(nPoints-1);
        updateBoundaryIndices();
    }

    return pointRemoved;
}

bool OneDimGrid::removeLeft(vector<dvector>& y)
{
    // *** Criteria for removal from the left (j==0) ***
    int djMom = 2;
    int djOther = (jb==1 && !fixedBurnedVal) ? 3 : 2;

    // Don't remove points if the location of the left boundary is fixed
    bool pointRemoved = true; // assume removal
    if (fixedLeftLoc) {
        if (pointRemoved && debugParameters::debugRegrid) {
            logFile.write("Regrid: left removal prevented by fixed left boundary");
        }
        pointRemoved = false;
    }

    for (size_t k=0; k<nAdapt; k++) {
        size_t dj = (k==kMomentum) ? djMom : djOther;
        double ymax = maxval(y[k]);
        if (abs(y[k][dj]-y[k][0])/ymax > boundaryTolRm && ymax > absvtol) {
            if (pointRemoved && debugParameters::debugRegrid) {
                logFile.write(format(
                    "Regrid: Left removal prevented component %i"
                    " flatness requirement: %g > %g") %
                    k % (abs(y[k][dj]-y[k][0])/ymax) % boundaryTolRm);
            }
            pointRemoved = false;
        }
    }

    if (jj < 3) {
        pointRemoved = false;
    }

    if (pointRemoved) {
        removePoint(0,y);
        setSize(nPoints-1);
        updateBoundaryIndices();
    }

    return pointRemoved;
}


void OneDimGrid::regrid(vector<dvector>& y)
{
    nVars = y.size();
    assert(nAdapt <= nVars);
    assert((dampVal > 0).all());

    setSize(y[0].size());

    bool rightAddition = addRight(y);
    bool leftAddition = addLeft(y);

    bool rightRemoval = false;
    bool leftRemoval = false;

    int rightRemovalCount = 0;

    if (!rightAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeRight(y);
            if (continueRemoval) {
                rightRemoval = true;
                rightRemovalCount++;
            }
        }
    }
    if (rightRemovalCount > 0) {
        std::string suffix = (rightRemovalCount == 1) ? "" : "s";
        logFile.write(format("Removed %i point%s from the right side.") %
                rightRemovalCount % suffix);
    }

    int leftRemovalCount = 0;
    if (!leftAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeLeft(y);
            if (continueRemoval) {
                leftRemoval = true;
                leftRemovalCount++;
            }
        }
    }

    if (leftRemovalCount > 0) {
        std::string suffix = (leftRemovalCount == 1) ? "" : "s";
        logFile.write(format("Removed %i point%s from the left side.") %
                leftRemovalCount % suffix);
    }

    updated = (updated || leftAddition || rightAddition ||
               leftRemoval || rightRemoval);

    if (updated) {
        updateValues();
        updateBoundaryIndices();
    }
}

void OneDimGrid::regridUnstrained(vector<dvector>& y, dvec& qdot)
{
    nVars = y.size();
    assert(nAdapt <= nVars);
    assert((dampVal > 0).all());

    setSize(y[0].size());

    bool rightAddition = addRightUnstrained(y, qdot);
    bool rightRemoval = false;
    int rightRemovalCount = 0;

    if (!rightAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeRightUnstrained(y, qdot);
            if (continueRemoval) {
                rightRemoval = true;
                rightRemovalCount++;
            }
        }
    }

    if (rightRemovalCount > 0) {
        std::string suffix = (rightRemovalCount == 1) ? "" : "s";
        logFile.write(format("Removed %i point%s from the right side.") %
                rightRemovalCount % suffix);
    }

    bool leftAddition = addLeft(y);
    bool leftRemoval = false;
    int leftRemovalCount = 0;

    if (!leftAddition) {
        bool continueRemoval = true;
        while (continueRemoval) {
            continueRemoval = removeLeft(y);
            if (continueRemoval) {
                leftRemoval = true;
                leftRemovalCount++;
            }
        }
    }

    if (leftRemovalCount > 0) {
        std::string suffix = (leftRemovalCount == 1) ? "" : "s";
        logFile.write(format("Removed %i point%s from the left side.") %
                leftRemovalCount % suffix);
    }

    updated = (updated || leftAddition || rightAddition ||
               leftRemoval || rightRemoval);

    if (updated) {
        updateValues();
        updateBoundaryIndices();
    }
}

void OneDimGrid::updateBoundaryIndices(void) {
    if (unburnedLeft) {
        ju = 0;
        jb = jj;
    } else {
        jb = 0;
        ju = jj;
    }
}

GridBased::GridBased()
    : x(grid.x)
    , r(grid.r)
    , rphalf(grid.rphalf)
    , hh(grid.hh)
    , dlj(grid.dlj)
    , cfm(grid.cfm)
    , cf(grid.cf)
    , cfp(grid.cfp)
    , alpha(grid.alpha)
    , beta(grid.beta)
    , nPoints(grid.nPoints)
    , jj(grid.jj)
{
}

void GridBased::setGrid(const OneDimGrid& other)
{
    grid = other;
}
