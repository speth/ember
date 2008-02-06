#include "grid.h"
#include "debugUtils.h"

using std::cout;
using std::endl;

void oneDimGrid::updateDerivedSizes()
{
	hh.resize(jj);
	cfm.resize(jj);
	cf.resize(jj);
	cfp.resize(jj);
	dlj.resize(jj);
	rphalf.resize(jj);

	for (unsigned int j=0; j<x.size()-1; j++) {
		hh[j] = x[j+1]-x[j];
		rphalf[j] =  pow(0.5*(x[j]+x[j+1]),alpha);
	}
	hh[jj-1] = hh[jj-2];

	for (unsigned int j=1; j<x.size()-1; j++) {
		cfp[j] = hh[j-1]/(hh[j]*(hh[j]+hh[j-1]));
		cf[j] = (hh[j]-hh[j-1])/(hh[j]*hh[j-1]);
		cfm[j] = -hh[j]/(hh[j-1]*(hh[j]+hh[j-1]));
		dlj[j]= 0.5 * (x[j+1]-x[j-1]);
	}
}

bool oneDimGrid::adapt(vector<dvector>& y)
{
	nVars = y.size();
	jj = y[0].size();

//! ************** subroutine : adaptation_x *************************
//! this subroutine takes unadapted solution vector, analyzes it and
//! returns an adapted solution vector. it tries to remove unnecessary
//! grid points located in the regions of small gradients, insert
//! new grid points in the regions of large gradients and, at the same
//! time maintain relative uniformity of the grid. it is a passive
//! adaption procedure.
//!
//! inputs:
//!
//! kmax/kk - max./actual number of species,
//! jmax/jj - maximum/actual mesh point number,
//! jzr- stagnation point index,
//! bottom- a numerical parameter. if max(y(k)) < bottom , y(k) is
//!         not to be adapted,
//! tol/ told - numerical tolerance parameters for insertion criterion
//!             based on the function/ its derivative,
//! aconst- a numerical parameters which determines the non-uniformity
//!          of the grid,
//! tmax- the maximum temperature in kelvins, 
//
//! inputs/ outputs:
//!
//! rhov(j)- mass flux in gm/cm^2-s,
//! hh(j)=x(j+1)-x(j)-  mesh sizes in cm,
//! x(j)- mesh points in cm, 
//! etamd(j)=(x(j+1)-x(j-1))/2 in cm,
//! t(j)- the temperature in kelvins,
//! y2(k,j)- mass fractions,
//! u(j)- non-dimensional x-component of velocity,
//! adapt- a flag that shows whether a change in the solution
//! vector has occured.
//!
//! adaptation algorithm
//! 
//! the insertion of the grid points is performed first.
//! solution vector is defined as (y2(k,j) t(j) rhov(j) u(j)).
//! here is the algorithm :
//! 1. take a component of the solution vector
//! 2. find its max/ min values and the max/min values of its 
//!    derivative.
//! 3. apply four criteria and and find where insertions are 
//!    needed. the criteria for a component f(j)
//!    a. |f(j)-f(j-1)| < tol* |max f(j)- minf(j)|
//!    b. |df(j)/dy -df(j-1)/dy | < tol* |max df(j)/dy - minf(j)/dy|
//!    c. 1/aconst < hh(j)/hh(j-1) < aconst
//! 4. if any of these criteria is not satisfied a grid point j
//!    is inserted in the array insrt(i2)
//! 5. perform actual adaptation of the whole solution vector
//! 6. take the next component and go to 1 until the whole 
//!    solution vector satisfies all the criteria.
//!
//! next, the unnecessary grid points are removed.
//! the algorithm is somewhat reversed. two regions are identified
//! inside the flame: the region of temperature gradient and the 
//! rest. for the former the newly introduced parameter tol2 is 
//! of order of 0.1, for the latter it is about 0.5.
//! after that the following criteria are applied to all components
//!    a. |f(j)-f(j-1)| > tol2* tol* |max f(j)- minf(j)|
//!    b. |df(j)/dy -df(j-1)/dy | > tol2* tol* |max df(j)/dy - minf(j)/dy|
//! if they are satisfied at some point, it's removed.
//! *******************************************************************
//subroutine adaptation_x (kmax,kk,jmax,jj,kn2,jzr, tmax, rhov,&
//     rrhov, hh, x, t, y2, u, y2n, y2nn, tn, tnn, un, unn, rrhovn, rrhovnn,srct,adapt)
//  
//  use define_kind
//  use flags
//  use transport_properties
//  use mixture_properties
//  use input_parameters
//  implicit none
//
//  integer kmax,k,jmax,j,jj,kk,i,ii
//  integer code,jzr,adapt
//  integer kn2, insert_done, insert, remove, ninsert, jinsert(jmax)
//  logical, parameter:: debug_output = .TRUE.
//
//  real(kind=real8):: srct(jmax)
//  real(kind=real8):: y2n(kmax,jmax), y2nn(kmax,jmax)
//  real(kind=real8):: tn(jmax), tnn(jmax), un(jmax), unn(jmax)
//  real(kind=real8):: rrhov(jmax), rrhovn(jmax), rrhovnn(jmax)
//  real(kind=real8):: rhov(jmax)
//
//  real(kind=real8):: x(jmax),r(jmax),t(jmax), xn(jmax)
//  real(kind=real8):: y2(kmax,jmax)
//  real(kind=real8):: cfp(jmax),cf(jmax),cfm(jmax), hh(jmax), dlj(jmax), rphalf(jmax)
//  real(kind=real8):: cfpn(jmax),cfn(jmax),cfmn(jmax), hhn(jmax)
//
//  real(kind=real8):: vt(jmax),minvt,maxvt,mindvt,maxdvt,dvt(jmax)
//  real(kind=real8):: vtrange, dvtrange
//  real(kind=real8):: tmax,hhmn,u(jmax), tol, told
//  real(kind=real8):: dampVal(jmax)
//
//  ! tol2 = 0.67 ! relative grid point removal tolerance
//  ! srctol = 0.75 ! relative tolerance for chemical source terms
//
//! a diagnostic file
//!  open (17,file='adaptation.out', status= 'unknown')

//  if (debug_adapt == 1) then
//     open(20,file='prof-adapt1',status='unknown')
//     open(21,file='spec-adapt1',status='unknown')
//     do j=1,jj
//        write(20,*) x(j), t(j), u(j), rrhov(j)
//        write(21,*) x(j), (y2(k,j), k=1,kk)
//     end do
//     close(20)
//     close(21)
//  end if
//
//!  gridmin = 1.5d-4
//!  gridmax = 0.2
//!  dampconst = 3.0 ! limit grid size based on relative strength of convection and diffusion

	bool adapt = false; // flag if adaption has occured

	//! ------------------------------
	//! GRID POINT INSERTION ALGORITHM
	//! ------------------------------

	int j = 0;
	vector<double> dv(jj);

	while (j < jj-1) {
		updateDerivedSizes();
		dv.resize(jj);
		bool insert = false;

		//dampVal(j) = min(minval(dkm(1:kk,j)), mu(j)/rho(j), lambda(j)/(rho(j)*cpmx(j))) &
		//     / (min(abs(rrhov(j)),abs(rhov(j)))/rho(j))

		// Consider tolerances each variable v in the solution y
		for (int k=0; k<nVars; k++) {

			vector<double>& v = y[k];
			for (int i=1; i<jj-1; i++) {
				dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
			}

			double vRange = mathUtils::range(v);

			if (vRange < absvtol) {
				continue;
			}

			double dvRange = mathUtils::range(dv,1,jj-2);

			// Apply grid point addition criteria:

			// resolution of v
			if (abs(v[j+1]-v[j]) > vtol*vRange) {
				insert = true;
				if (debugParameters::debugAdapt) {
					cout << "Adapt: v resolution wants grid point j = " << j << ", k = " << k;
					cout << " |v(j+1)-v(j)|/vrange = " << abs(v[j+1]-v[j])/vRange << " > " << vtol << endl;
				}
			}

			// resolution of dv
			if (j!=0 && j!=jj-2 && abs(dv[j+1]-dv[j]) > dvtol*dvRange) {
				insert = true;
				if (debugParameters::debugAdapt) {
					cout << "Adapt: dv resolution wants grid point j = " << j << ", k = " << k;
					cout << " |dv(j+1)-dv(j)|/vrange = " << abs(dv[j+1]-dv[j])/dvRange << " > " << dvtol << endl;
				}
			}

		}
		
		// Left uniformity
		if (j!=0 && hh[j]/hh[j-1] > uniformityTol) {
			insert = true;
			if (debugParameters::debugAdapt) { 
				cout << "Adapt: left uniformity wants grid point j = " << j;
				cout << " hh(j)/hh(j-1) = " << hh[j]/hh[j-1] << ' > ' << uniformityTol << endl;
			}
		}

		// Right uniformity
		if (j!=jj && hh[j]/hh[j+1] > uniformityTol) {
			insert = true;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: right uniformity wants grid point j = " << j;
				cout << " hh(j)/hh(j+1) = " << hh[j]/hh[j+1] << " > " << uniformityTol << endl;
			}
		}

		// Maximum grid size
		if (hh[j] > gridMax) {
			insert = true;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: Maximum grid size criterion wants a grid point j = " << j;
				cout << " hh(j) = " << hh[j] << " > " << gridMax << endl;
			}
		}

		// Damping of high-frequency numerical error
		if (hh[j] > dampConst*dampVal[j]) {
			insert = true;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: damping criterion wants a grid point j = " << j;
				cout << " hh(j) = " << hh[j] << " > " << dampConst*dampVal[j] << endl;
			}
		}
		
		// Enforce minimum grid size
		if (hh[j] < 2*gridMin) {
			insert = false;
			if (debugParameters::debugAdapt) {
				cout << "Adapt: grid point addition cancelled by minimum grid size j = " << j;
				cout << " hh(j) = " << hh[j] << " < " << 2*gridMin << endl;
			}
		}

		if (insert) {
			// Insert a new point
			cout << "Adapt: inserting grid point j = " << j << endl;
			addPoint(j, y);
			adapt = true;
			jj++;
		} else {
			// No insertion; step to the next point.
			j++;
		}
				

	}

	// ----------------------------    
	// GRID POINT REMOVAL ALGORITHM
	// ----------------------------

	j = 1;
	while (j < jj-1) {
		updateDerivedSizes();

		// Assume removal, then look for a conditon which prevents removal
		bool remove = true;

		// Consider tolerances each variable v in the solution y
		for (int k=0; k<nVars; k++) {

			vector<double>& v = y[k];
			for (int i=1; i<jj-1; i++) {
				dv[i] = cfp[i]*v[i+1] + cf[i]*v[i] + cfm[i]*v[i-1];
			}

			double vRange = mathUtils::range(v);
			double dvRange = mathUtils::range(dv,1,jj-2);

			if (vRange < absvtol) {
				continue; // ignore minor componenents
			}

			// Apply grid point removal criteria:

			// resolution of v
			if (abs(v[j+1]-v[j-1]) > rmTol*vtol*vRange) {
				if (debugParameters::debugAdapt) {
					cout << "Adapt: no removal - v res. j = " << j << ", k = " << k;
					cout << " |v(j+1)-v(j-1)|/vtrange = " << abs(v[j+1]-v[j-1])/vRange << " > " << vtol*rmTol << endl;
				}
				remove = false;
			}

			if (j!=2 && j!=jj-2 && abs(dv[j+1]-dv[j-1]) > rmTol*dvtol*dvRange) {
				if (debugParameters::debugAdapt) {
					cout << "Adapt: no removal - dv res. j = " << j << ", k = " << k;
					cout << " |dv(j+1)-dv(j-1)|/dvrange = " << abs(dv[j+1]-dv[j-1])/dvRange << " > " << dvtol*rmTol << endl;
				}
				remove = false;
			}
		}

		// Enforce maximum grid size
		if (hh[j]+hh[j-1] > gridMax) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - maximum grid size. j = " << j;
				cout << " hh(j)+hh(j-1) = " << hh[j]+hh[j-1] << " > " << gridMax << endl;
			}
			remove = false;
		}

		// Damping of high-frequency numerical error
		if (hh[j]+hh[j-1] >= rmTol*dampConst*dampVal[j]) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - damping criterion. j = " << j;
				cout << "hh(j) = " << hh[j] << " > " << dampConst*dampVal[j] << endl;
			}
			remove = false;
		}

		// Enforce left uniformity
		if (j>=2 && hh[j]+hh[j-1] > uniformityTol*hh[j-2]) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - left uniformity. j = " << j;
				cout << "(hh(j)+hh(j-1))/hh(j-2) = " << (hh[j]+hh[j-1])/hh[j-2] << " > " << uniformityTol << endl;
			}
			remove = false;
		}

		// Enforce right uniformity
		if (j<=jj-3 && hh[j]+hh[j-1] > uniformityTol*hh[j+1]) {
			if (debugParameters::debugAdapt) {
				cout << "Adapt: no removal - right uniformity. j = " << j;
				cout << "(hh(j)+hh(j-1))/hh(j+1) = " << (hh[j]+hh[j-1])/hh[j+1] << " > " << uniformityTol << endl;
			}
			remove = false;
		}

		if (remove) {
			cout << "Adapt: removing grid node j = " << j << endl;
			removePoint(j, y);
            jj--;
		    adapt = true;
		} else {
			j++;
		}
	}

	updateDerivedSizes();

//! if the mesh has been changed, recalculate the location
//! of the stagnation point
//! keep the stag point where x = 0. this is the "anchor" for 
//! the computational domain.
//! this search routine provides jzr, which maintains the
//! boundary condition on rhov whenver the continuity equation
//! is integrated
//
//  if (adapt /= 0) then
//      jzr = minval(minloc(abs(rrhov(2:jj)))) + 1
//      if (sign(1.0d0,rrhov(jzr)) == sign(1.0d0,rrhov(jzr+1))) then
//         jzr = jzr - 1
//      end if
//   end if

   //if (debug_adapt == 1) then
   //   open(20,file='prof-adapt2',status='unknown')
   //   open(21,file='spec-adapt2',status='unknown')
   //   do j=1,jj
   //      write(20,*) x(j), t(j), u(j), rrhov(j)
   //      write(21,*) x(j), (y2(k,j), k=1,kk)
   //   end do
   //   close(20)
   //   close(21)
   //end if

//!   close(17) ! adaptation.out
	return adapt;

}

void oneDimGrid::addPoint(int jInsert, vector<dvector>& y)
{
  	dvector::iterator iter = x.begin() + jInsert;
	x.insert(iter+1,0.5*(x[jInsert+1]+x[jInsert]));

	iter = dampVal.begin() + jInsert;
	dampVal.insert(iter+1,0.5*(dampVal[jInsert+1]+dampVal[jInsert]));

	vector<dvector>::iterator i;
	for (i=y.begin(); i!=y.end(); i++) {
		iter = i->begin();
		iter += jInsert;
		i->insert(iter+1, 0.5*(*(iter+1) + *(iter)));
	}
}

void oneDimGrid::removePoint(int jRemove, vector<dvector>& y)
{
	x.erase(x.begin() + jRemove);
	dampVal.erase(dampVal.begin() + jRemove);

	vector<dvector>::iterator i;
	for (i=y.begin(); i!=y.end(); i++) {
		i->erase(i->begin() + jRemove);
	}
}

bool oneDimGrid::regrid(std::vector<dvector>& solutionState)
{
	vector<dvector>& y = solutionState;
	//! ********* subroutine : regrid_x *********************************
//! this subroutine compares the thickness of the flame based on the 
//! temperature profile with the distance from the flame to the left
//! and right boundaries of the domain. if the flame is approaching
//! either of these boundaries some new points are added and the 
//! solution vector is re-calculated.
//!
//! inputs:
//!
//! kmax/kk - max./actual number of species,
//! jmax/jj - maximum/actual mesh point number,
//! jzr- stagnation point index,
//! xl/ xr- left/right flame boundaries in cm, 
//! dlt- the flame thickness in cm,
//! hh(j)=x(j+1)-x(j)-  mesh sizes in cm,
//! dlj(j)=(x(j+1)-x(j-1))/2 in cm,
//! tsuji=(1/0)- a configuration parameter
//!
//! inputs/ outputs :
//! x(j)- mesh points in cm, 
//! t(j)- the temperature in kelvins,
//! y2(k,j)- mass fractions,
//! u(j)- non-dimensional x-component of velocity,
//! rhov(j)- mass flux in gm/cm^2-s,
//! rgrd- a flag that shows wether a change in the mesh
//! vector has occured.
//! ****************************************************************
//

  //integer jzr, jmax, jj,j, kmax, kk, k, rgrd

  //real(kind=real8):: x(jmax), t(jmax), y2(kmax,jmax)
  //real(kind=real8):: un(jmax), tn(jmax), rhovn(jmax)
  //real(kind=real8):: unn(jmax), tnn(jmax), rhovnn(jmax)
  //real(kind=real8):: y2n(kmax,jmax), y2nn(kmax,jmax)
  //real(kind=real8):: u(jmax), hh(jmax), dlj(jmax)
  //real(kind=real8):: xl, xr, dlt, hmean, rhov(jmax)
  //real(kind=real8):: det, detl, detr, dxl, dxr
  //real(kind=real8):: rhon(jmax)
  //
  //integer addright, rmright, addleft, rmleft
  //integer dj, dju ! distance to comparison point for smoothness criteria

	bool gridUpdated = false;

	/*

  // ------------------------------------------------
  //       Criteria for addition to right (j=jj)
  // ------------------------------------------------

  if (jb == jj .AND. fixedBurnedValFlag == 0) then
     dj = 2 ! condition at j=jj is zero-gradient
  else
     dj = 1 ! condition at j=jj is fixed value
  end if

  if (jb == jj) then ! dj for momemtum
     dju = 2
  else
     dju = 1
  end if

  addright = 0 ! assume no addition
  ! compare flame thickness with distance between flame and boundary
  if ( abs(x(jj)-xr) < addtol*dlt) then
     if (addright == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid1: flame location wants right addition:', x(jj)-xr, '<', addtol*dlt
     end if
     addright = 1 
  end if

  ! compare flame thickness with distance between stagnation point and boundary
  if ( abs(x(jj)-x(jzr)) < addtol*dlt) then
     if (addright == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid1: stagnation pt location wants right addition:', x(jj)-x(jzr),'<', addtol*dlt
     end if
     addright = 1
  end if

  ! check flatness of temperature, velocity and species profiles at the boundary
  if ( abs(t(jj)-t(jj-dj))/(maxval(t(1:jj))-minval(t(1:jj))) > boundarytol ) then 
     if (addright == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid1: Temperature flatness wants right addition:', &
             abs(t(jj)-t(jj-1))/(maxval(t(1:jj))-minval(t(1:jj))), '>', boundarytol
     end if
     addright = 1
  end if

  if ( abs(u(jj)-u(jj-dju))/(maxval(u(1:jj))-minval(u(1:jj))) > boundarytol ) then
     if (addright == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid1: Velocity flatness wants right addition:', &
             abs(u(jj)-u(jj-dju))/(maxval(u(1:jj))-minval(u(1:jj))), '>', boundarytol
     end if
     addright = 1
  end if

  do k=1,kk
     if ( abs(y2(k,jj)-y2(k,jj-dj))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))) > boundarytol & 
          .AND. maxval(y2(k,1:jj)) > bottom) then
        if (addright == 0 .and. debug_regrid == 1) then
           write(*,*) 'regrid1: species', k, 'flatness wants right addition:', &
                abs(y2(k,jj)-y2(k,jj-dj))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))), '>', boundarytol
        end if
        addright = 1
     end if
  end do

  if (addright == 0 .and. debug_regrid == 1) then
     write(*,*) 'regrid1: No addition on right'
  end if


  ! ----------------------------------------------
  !       Criteria for addition to left (j=1)
  ! ----------------------------------------------
  if (jb == 1 .AND. fixedBurnedValFlag == 0) then
     dj = 2 ! condition at j=jj is zero-gradient
  else
     dj = 1 ! condition at j=jj is fixed value
  end if

  if (jb == 1) then ! dj for momemtum
     dju = 2
  else
     dju = 1
  end if

  addleft = 0 ! assume no addition

  ! compare flame thickness with distance between flame and boundary
  if (abs(xl-x(1)) < addtol*dlt) then 
     if (addleft == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid2: flame location wants left addition:', xl-x(1), '<', addtol*dlt
     end if
     addleft = 1
  end if

  ! compare flame thickness with distance between stagnation point and boundary
  if (abs(x(jzr)-x(1)) < addtol*dlt) then 
     if (addleft == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid2: stagnation pt location wants left addition:', x(jzr)-x(1),'<', addtol*dlt
     end if
     addleft = 1
  end if

  ! check flatness of temperature, velocity and species profiles at the boundary
  if ( abs(t(1+dj)-t(1))/(maxval(t(1:jj))-minval(t(1:jj))) > boundarytol )   then 
     if (addleft == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid2: Temperature flatness wants left addition:', &
             abs(t(1+dj)-t(1))/(maxval(t(1:jj))-minval(t(1:jj))), '>', boundarytol
     end if
     addleft = 1
  end if

  if ( abs(u(1+dju)-u(1))/(maxval(u(1:jj))-minval(u(1:jj))) > boundarytol ) then
     if (addleft == 0 .and. debug_regrid == 1) then
        write(*,*) 'regrid2: Velocity flatness wants left addition:', &
             abs(u(1+dju)-u(1))/(maxval(u(1:jj))-minval(u(1:jj))), '>', boundarytol
     end if
     addleft = 1
  end if

  do k=1,kk
     if ( abs(y2(k,1+dj)-y2(k,1))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))) > boundarytol & 
          .AND. maxval(y2(k,1:jj)) > bottom) then
        if (addleft == 0 .and. debug_regrid == 1) then
           write(*,*) 'regrid2: Species', k, 'flatness wants left addition:', &
                abs(y2(k,1+dj)-y2(k,1))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))), '>', boundarytol
        end if
        addleft = 1
     end if
  end do

  if (fixedLeftLocFlag == 1 .AND. x(1) /= 0) then
     if (addleft == 0 .AND. debug_regrid == 1) then
        write(*,*) 'regrid2: Adding point to force left boundary towards x=0'
     end if
     addleft = 1
  end if

  ! don't add points if the symmetry boundary is already at the edge of the domain
  if (abs(x(1)) <= gridmin) then
     if (addleft == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid2: Cannot add point to left because x(1) = 0.'
     end if
     addleft = 0
  end if

  if (addleft == 0 .and. debug_regrid == 1) then
     write(*,*) 'regrid2: No addition on left'
  end if


  ! ----------------------------------------------
  !     Critera for removal from right (j=jj)
  ! ----------------------------------------------
  if (jb == jj .AND. fixedBurnedValFlag == 0) then
     dj = 3 ! condition at j=jj is zero-gradient
  else
     dj = 2 ! condition at j=jj is fixed value
  end if

  if (jb == jj) then ! dj for momemtum (always zero-gradient on burned side)
     dju = 3
  else
     dju = 2
  end if


  rmright = 1 ! assume removal
  if (abs(x(jj-1)-xr) < rmtol*dlt) then 
     if (rmright == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid3: right removal prevented by flame location:', x(jj-1)-xr, '<', rmtol*dlt
     end if
     rmright = 0
  end if

  if (abs(x(jj-1)-x(jzr)) < rmtol*dlt) then
     if (rmright == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid3: right removal prevented by stagnation point location:', &
             x(jj-1)-x(jzr), '<', rmtol*dlt
     end if
     rmright = 0
  end if

  if ( abs(t(jj)-t(jj-dj))/(maxval(t(1:jj))-minval(t(1:jj))) > btolrm ) then
     if (rmright == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid3: right removal prevented by temperature flatness:', &
             abs(t(jj)-t(jj-dj))/(maxval(t(1:jj))-minval(t(1:jj))), '>', btolrm
     end if
     rmright = 0
     end if

  if ( abs(u(jj)-u(jj-dju))/(maxval(u(1:jj))-minval(u(1:jj))) > btolrm ) then
     if (rmright == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid3: right removal prevented by velocity flatness:', &
             abs(u(jj)-u(jj-dju))/(maxval(u(1:jj))-minval(u(1:jj))), '>', btolrm
     end if
     rmright = 0
  end if

  do k=1,kk
     if ( abs(y2(k,jj)-y2(k,jj-dj))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))) > btolrm & 
          .AND. maxval(y2(k,1:jj)) > bottom) then
        if (rmright == 1 .and. debug_regrid == 1) then
           write(*,*) 'regrid3: right removal prevented by species', k, 'flatness:', &
                abs(y2(k,jj)-y2(k,jj-dj))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))), '>', btolrm
        end if
        rmright = 0
     end if
  end do

  if (rmright == 1 .and. debug_regrid == 1) then
     write(*,*) 'regrid3: Removal on right allowed'
  end if

  ! ---------------------------------------------
  !      Critera for removal from left (j=1)
  ! ---------------------------------------------
  if (jb == 1 .AND. fixedBurnedValFlag == 0) then
     dj = 3 ! condition at j=jj is zero-gradient
  else
     dj = 2 ! condition at j=jj is fixed value
  end if

  if (jb == 1) then ! dj for momemtum
     dju = 3
  else
     dju = 2
  end if

  rmleft = 1 ! assume removal

  ! don't remove points if the location of the left boundary is fixed
  if (fixedLeftLocFlag == 1) then
     if (rmleft == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid4: left removal prevented by fixed left boundary'
     end if
     rmleft = 0
  end if

  if (abs(xl-x(2)) < rmtol*dlt) then 
     if (rmleft == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid4: left removal prevented by flame location:', xl-x(2), '<', rmtol*dlt
     end if
     rmleft = 0
  end if

  if (abs(x(jzr)-x(2)) < rmtol*dlt) then
     if (rmleft == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid4: left removal prevented by stagnation point location:', &
             abs(x(jzr)-x(2)), '<', rmtol*dlt
     end if
     rmleft = 0
  end if

  if ( abs(t(1+dj)-t(1))/(maxval(t(1:jj))-minval(t(1:jj))) > btolrm ) then
     if (rmleft == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid4: left removal prevented by temperature flatness:', &
             abs(t(1+dj)-t(1))/(maxval(t(1:jj))-minval(t(1:jj))), '>', btolrm
     end if
     rmleft = 0
  end if
             
  if ( abs(u(1+dju)-u(1))/(maxval(u(1:jj))-minval(u(1:jj))) > btolrm ) then
     if (rmleft == 1 .and. debug_regrid == 1) then
        write(*,*) 'regrid4: left removal prevented by velocity flatness:', &
             abs(u(1+dju)-u(1))/(maxval(u(1:jj))-minval(u(1:jj))), '>', btolrm
     end if
     rmleft = 0
  end if

  do k=1,kk
     if ( abs(y2(k,1+dj)-y2(k,1))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))) > btolrm & 
          .AND. maxval(y2(k,1:jj)) > bottom) then ! species flatness
        if (rmleft == 1 .and. debug_regrid == 1) then
           write(*,*) 'regrid4: left removal prevented by species', k, 'flatness:', &
                abs(y2(k,1+dj)-y2(k,1))/(maxval(y2(k,1:jj))-minval(y2(k,1:jj))), '>', btolrm
        end if
        rmleft = 0
     end if
  end do

  if (rmleft == 1 .and. debug_regrid == 1) then
     write(*,*) 'regrid4: Removal on left allowed'
  end if


  ! *** Perform the indicated additions and removals ***
  if ( addright == 1) then
     ! ---------------------------------
     !           add to right
     ! ---------------------------------
     rgrd = 1
     write(6,*) 'regrid5: regrid right immediately ! (addition) '

     ! location for the new point
!     x(jj+1) = x(jj)+ max( min(x(jj)-x(jj-1),0.5*dlt) , 0.2*dlt)
     x(jj+1) = x(jj) + (x(jj)-x(jj-1))

     ! linear extrapolation for rhov
     rhov(jj+1) = rhov(jj) + (rhov(jj)-rhov(jj-1)) / (x(jj)-x(jj-1)) * (x(jj+1)-x(jj))
     rhovn(jj+1) = rhovn(jj) + (rhovn(jj)-rhovn(jj-1)) / (x(jj)-x(jj-1)) * (x(jj+1)-x(jj))
     rhovnn(jj+1) = rhovnn(jj) + (rhovnn(jj)-rhovnn(jj-1)) / (x(jj)-x(jj-1)) * (x(jj+1)-x(jj))

     ! everything else is constant
     t(jj+1) = t(jj)
     u(jj+1) = u(jj)
     un(jj+1) = un(jj)
     unn(jj+1) = unn(jj)
     tn(jj+1) = tn(jj)
     tnn(jj+1) = tnn(jj)
     rhovn(jj+1) = rhovn(jj)
     rhovnn(jj+1) = rhovnn(jj)
     rhon(jj+1) = rhon(jj)
     y2(1:kk,jj+1) = y2(1:kk,jj)
     y2n(1:kk,jj+1) = y2n(1:kk,jj)
     y2nn(1:kk,jj+1) = y2nn(1:kk,jj)

     jj=jj+1

  elseif ( rmright == 1 ) then
     ! ---------------------------------
     !       remove from right
     ! ---------------------------------
     rgrd = 1
     write(*,*) 'regrid5: regrid right immediately! (removal)'

     ! move data from the rightmost point over, except for rhov and x
     t(jj-1) = t(jj)
     u(jj-1) = u(jj)
     un(jj-1) = un(jj)
     unn(jj-1) = unn(jj)
     tn(jj-1) = tn(jj)
     tnn(jj-1) = tnn(jj)
     rhovn(jj-1) = rhovn(jj)
     rhovnn(jj-1) = rhovnn(jj)
     rhon(jj-1) = rhon(jj)
     y2(1:kk,jj-1) = y2(1:kk,jj)
     y2n(1:kk,jj-1) = y2n(1:kk,jj)
     y2nn(1:kk,jj-1)=y2nn(1:kk,jj)

     jj = jj-1
  end if

  if ( addleft == 1) then
     ! ---------------------------------
     !           add to left
     ! ---------------------------------
     rgrd = 1
     write(6,*) 'regrid6: regrid left immediately! (addition)'

     ! move everything to the right one index
     ! the values remaining at the first index are correct except for x and rhov
     do j=jj,1,-1
        t(j+1) = t(j)
        rhov(j+1) = rhov(j)
        u(j+1) = u(j)
        x(j+1) = x(j)
        un(j+1) = un(j)
        unn(j+1) = unn(j)
        tn(j+1) = tn(j)
        tnn(j+1) = tnn(j)
        rhovn(j+1) = rhovn(j)
        rhovnn(j+1) = rhovnn(j)
        rhon(j+1) = rhon(j)
        y2(1:kk,j+1)=y2(1:kk,j)
        y2n(1:kk,j+1)=y2n(1:kk,j)
        y2nn(1:kk,j+1)=y2nn(1:kk,j)
     end do

     ! new x(1)
!     x(1)=x(2)-max( min(x(3)-x(2),0.5*dlt) , 0.2*dlt)
     x(1) = x(2) - (x(3)-x(2))
     if (x(2) > 0 .and. x(1) < gridmin) x(1) = 0.0d0

     rhov(1)= rhov(2) - (rhov(3)-rhov(2)) / (x(3)-x(2)) * (x(2)-x(1))
     rhovn(1)= rhovn(2) - (rhovn(3)-rhovn(2)) / (x(3)-x(2)) * (x(2)-x(1))
     rhovnn(1)= rhovnn(2) - (rhovnn(3)-rhovnn(2)) / (x(3)-x(2)) * (x(2)-x(1))
     
     jj = jj+1

  elseif ( rmleft == 1) then
     ! ---------------------------------
     !         remove from left
     ! ---------------------------------
     rgrd = 1
     write(*,*) 'regrid6: regrid left immediately! (removal)'

     ! shift the boundary values of rhov and x
     x(1) = x(2)
     rhov(1) = rhov(2)
     rhovn(1) = rhovn(2)
     rhovnn(1) = rhovnn(2)

     ! shift the indicies, but keep the value at the boundary
     do j=2,jj
        t(j) = t(j+1)
        tn(j) = tn(j+1)
        tnn(j) = tnn(j+1)
        u(j) = u(j+1)
        un(j) = un(j+1)
        unn(j) = unn(j+1)
        rhov(j) = rhov(j+1)
        rhovn(j) = rhovn(j+1)
        rhovnn(j) = rhovnn(j+1)
        rhon(j) = rhon(j+1)
        y2(1:kk,j) = y2(1:kk,j+1)
        y2n(1:kk,j) = y2n(1:kk,j+1)
        y2nn(1:kk,j) = y2nn(1:kk,j+1)
        x(j) = x(j+1)
     end do

     jj = jj-1

  end if

  ! if the domain has been expanded, re-calculate the location
  ! of the stagnation point and mesh deltas.
  if (rgrd == 1) then

     jzr = minval(minloc(abs(rhov(2:jj)))) + 1
     if (sign(1.0d0,rhov(jzr)) == sign(1.0d0,rhov(jzr+1))) then
        jzr = jzr - 1
     end if

     hh(1) = x(2)-x(1)
     do j=2,jj-1
        dlj(j) = 0.5 * (x(j+1)-x(j-1))
        hh(j) = x(j+1)-x(j)
     end do

  end if
	*/
	return gridUpdated;

}

void oneDimGrid::update_jZero(void)
{
      //jzr = minval(minloc(abs(rrhov(2:jj)))) + 1
      //if (sign(1.0d0,rrhov(jzr)) == sign(1.0d0,rrhov(jzr+1))) then
      //   jzr = jzr - 1
      //end if
}