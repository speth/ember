#include "strainedFlame.h"
#include <iostream>
#include <vector>

#include <cmath>

typedef vector<double> dvector;

// All Cantera names are in namespace Cantera. You can either
// reference everything as Cantera::<name>, or include the following
// 'using namespace' line.
//using namespace Cantera;

int main(int argc, char** argv) {
		
	char* python_cmd = getenv("PYTHON_CMD");
	if (!python_cmd) {
		putenv("PYTHON_CMD=python");
	}
    try {
    	strainedFlame();
    }
	catch (Cantera::CanteraError) {
		Cantera::showErrors(cout);
    }
    
	return 0;
}

void strainedFlame() {

    cout << "**** strainedFlame (1dflame Version 2.0) ****" << std::endl;

	// output file:
	ofstream outFile;
	outFile.open("out.m");
	clock_t t1, t2;
	t1 = clock();

    strainedFlameSys theSys;
	theSys.nPoints = 1000;
	//theSys.nPoints = 10;
	theSys.xLeft = -0.05;
	theSys.xRight = 0.05;
	theSys.Tleft = 300;
	theSys.Tright = 600;
	theSys.strainRate = 20;
	theSys.mu = 1.8e-5;
	theSys.lambda = 0.0243;
	theSys.cp = 1.005;
	theSys.rhoLeft = 1.20;
	theSys.tStart = 0;
	theSys.tEnd = 0.2;
	// Initial Conditions for ODE
	theSys.setup();

	// Sundials IDA Solver:
	sundialsIDA theSolver(theSys.N);
	theSolver.reltol = 1e-5;
	theSolver.nRoots = 0;
	theSolver.findRoots = false;

	// Initial condition:
	theSys.rollY(theSolver.y);
	for (int i=0; i<theSys.nPoints; i++) {
		theSolver.abstol(3*i) = 1e-6;
		theSolver.abstol(3*i+1) = 1e-6;
		theSolver.abstol(3*i+2) = 1e-6;
		theSolver.componentId(3*i) = 0.0;
		theSolver.componentId(3*i+1) = 1.0;
		theSolver.componentId(3*i+2) = 1.0;
	}
	theSolver.componentId(0) = 1;

	for (int i=0; i<theSys.N; i++) {
		theSolver.ydot(i) = 0;
	}
	theSolver.t0 = theSys.tStart;

	theSolver.setDAE(&theSys);
	theSolver.initialize();
	

	// In loop, call CVode and print results
	//   Break out of loop when NOUT preset output times have been reached.

	int iout = 0;  
	double dt = theSys.tEnd/100;
	double t = dt;
	int i=1;
	
	int flag;
	double R = Cantera::GasConstant;
	theSys.printForMatlab(outFile, theSys.x, 1, "x");
	theSys.unrollY(theSolver.y);
	theSys.printForMatlab(outFile, theSys.rhov, i, "rhov");
	theSys.printForMatlab(outFile, theSys.T, i, "T");
	theSys.printForMatlab(outFile, theSys.U, i, "U");
	outFile << "t(" << i << ") = " << 0 << ";" << endl;
	while(t <= theSys.tEnd) {
		
		flag = theSolver.integrateToTime(t);

		if (flag == CV_SUCCESS) {
			i++;
			
			theSys.unrollY(theSolver.y);
			theSys.printForMatlab(outFile, theSys.rhov, i, "rhov");
			theSys.printForMatlab(outFile, theSys.T, i, "T");
			theSys.printForMatlab(outFile, theSys.U, i, "U");
			outFile << "t(" << i << ") = " << t << ";" << endl;
			cout << "t(" << i << ") = " << t << ";" << endl;
			t += dt;
		}
	}
	theSolver.printStats();
	t2 = clock();
	cout << "Runtime: " << ((double)(t2-t1))/CLOCKS_PER_SEC << " seconds." << endl;
	
	outFile.close();
 	int blargh = 0;
}

void strainedFlameSys::setup(void)
{
	N = 3*nPoints;
	jacBW = 3;
	jacBWdot = 1;
	nVars = 3;

	T.resize(nPoints);
	U.resize(nPoints);
	rhov.resize(nPoints);

	dTdt.resize(nPoints,0);
	dUdt.resize(nPoints,0);
	drhovdt.resize(nPoints,0);

	resContinuity.resize(nPoints);
	resMomentum.resize(nPoints);
	resEnergy.resize(nPoints);

	rho.resize(nPoints);
	drhodt.resize(nPoints);

	x.resize(nPoints);
	dx.resize(nPoints);

	dFdy.resize(N, vector<double> (N));
	dFdydot.resize(N, vector<double> (N));

	bandedJacobian = new sdBandMatrix(N,nVars+1,nVars+1,2*nVars+2);
	BandZero(bandedJacobian->forSundials());

	pMat.resize(N);

	// Grid and initial profiles of T, U and rhov
	for (int i=0; i<nPoints; i++) {
		x[i] = xLeft + (xRight-xLeft)*((double) i)/((double) nPoints);
		T[i] = Tleft + (Tright-Tleft)*(x[i]-xLeft)/(xRight-xLeft);
		rho[i] = rhoLeft*Tleft/T[i];
		U[i] = sqrt(rho[0]/rho[i]);
	}
	rhov[0] = -strainRate*rho[0]*x[0];
	for (int i=1; i<nPoints; i++) {
		rhov[i] = rhov[i-1] - rho[i]*U[i]*strainRate*(x[i]-x[i-1]);
		dx[i-1] = x[i]-x[i-1];
	}

}

int strainedFlameSys::f(realtype t, sdVector& y, sdVector& ydot, sdVector& res)
{
	unrollY(y);
	unrollYdot(ydot);

	// Update auxillary data:
	for (int i=0; i<nPoints; i++) {
		rho[i] = rhoLeft*Tleft/T[i];
		drhodt[i] = rho[i]/T[i]*dTdt[i];
	}

	int jj = nPoints-1;

	// Left boundary:
	resContinuity[0] = drhovdt[0];
	resEnergy[0] = dTdt[0];
	resMomentum[0] = dUdt[0];

	// Intermediate points (energy & momentum equations)
	for (int i=1; i<jj; i++) {
		resEnergy[i] = rhov[i]*(T[i+1]-T[i-1])/(dx[i-1]+dx[i])
			- lambda/cp*(T[i+1]-2*T[i]+T[i-1])/(dx[i-1]*dx[i])
			+ rho[i]*dTdt[i];

		resMomentum[i] = rhov[i]*(U[i+1]-U[i-1])/(dx[i-1]+dx[i])
			- mu*(U[i+1]-2*U[i]+U[i-1])/(dx[i-1]*dx[i])
			+ rho[i]*dUdt[i]
			+ rho[i]*U[i]*U[i]*strainRate
			- rhoLeft*strainRate;
	}

	// Intermediate points (continuity equation)
	for (int i=1; i<nPoints; i++) {
		resContinuity[i] = drhodt[i] + rho[i]*U[i]*strainRate
			+ (rhov[i]-rhov[i-1])/dx[i-1];
	}

	// Right Boundary: fixed values
	resEnergy[jj] = dTdt[jj];
	resMomentum[jj] = dUdt[jj];

	rollResiduals(res);
	double resNorm = 0;
	double resC = 0;
	double resE = 0;
	double resM = 0;
	for (int i=0; i<N; i++) {
		resNorm += res(i)*res(i);
	}
	for (int i=0; i<nPoints; i++) {
		resC += resContinuity[i]*resContinuity[i];
		resE += resEnergy[i]*resEnergy[i];
		resM += resMomentum[i]*resMomentum[i];
	}
	resNorm = sqrt(resNorm);
	resC = sqrt(resC);
	resE = sqrt(resE);
	resM = sqrt(resM);
	//cout << "resNorm = " << resNorm << endl;
	//cout << "resC = " << resC << endl;
	//cout << "resE = " << resE << endl;
	//cout << "resM = " << resM << endl;

	return 0;
}

int strainedFlameSys::JvProd(realtype t, sdVector& y, sdVector& ydot,
			                 sdVector& res, sdVector& v, sdVector& Jv, realtype c_j)
{
	sdVector resTemp(N);
	sdVector yTemp(N);
	double yNorm = 0;
	double vNorm = 0;
	double errNorm = 0;
	for (int i=0; i<N; i++) {
		yNorm += y(i)*y(i);
		vNorm += v(i)*v(i);
	}
	yNorm = sqrt(yNorm);
	vNorm = sqrt(vNorm);
	double eps = pow((1+yNorm)*(DBL_EPSILON),0.3333)/vNorm;
	for (int i=0; i<N; i++) {
		yTemp(i) = y(i) + v(i)*eps;
	}
	f(t, yTemp, ydot, resTemp);
	sdVector JvAlt(N);
	for (int i=0; i<N; i++) {
		JvAlt(i) = (resTemp(i) - res(i))/eps;
		Jv(i) = 0;
		for (int j=0; j<N; j++)
			Jv(i) += dFdy[i][j]*v(j) + dFdydot[i][j]*c_j*v(j);
	}

	//for (int i=0; i<N; i++) {
	//	errNorm += pow(JvAlt(i) - Jv(i),2);
	//}
	//errNorm = sqrt(errNorm);
	//cout << "errNorm = " << errNorm << endl;
	return 0;
}

int strainedFlameSys::preconditionerSetup(realtype t, sdVector& y, sdVector& ydot, 
										  sdVector& res, realtype c_j)
{
	double eps = sqrt(DBL_EPSILON);
	sdVector yTemp(N);
	sdVector ydotTemp(N);
	sdVector resTemp(N);
	BandZero(bandedJacobian->forSundials());

	// J = dF/dy
	// Banded, upper & lower bandwidths = nVars. 
	int bw = 2*nVars+1;
	for (int s=0; s<bw; s++)
	{
		for (int i=0; i<N; i++) {
			yTemp(i) = y(i);
		}
		// Form the perturbed y vector
		for (int i=s; i<N; i+=bw) {
			
			if (y(i) != 0) {
				yTemp(i) += yTemp(i)*eps;
			} else {
				yTemp(i) = eps;
			}
		}
		// Evaluate the residuals
		f(t, yTemp, ydot, resTemp);
		
		// Fill in the corresponding Jacobian elements
		int iStart = max(0,s-nVars);
		int counter;
		int col;
		//int counter = max(1,nVars-s);
		if (s >= nVars) {
			counter = 0;
			col = -bw;
			
		} else {
			counter = nVars - s;
			col = 0;
		}
		
		for (int i=iStart; i<N; i++) {
			if (counter++ % bw == 0) {
				col += bw;
			}
			if (col+s>=N) {
				continue;
			}
			(*bandedJacobian)(i,s+col) = (resTemp(i)-res(i))/(yTemp(s+col)-y(s+col));
			//outFile << "Jv2(" << i+1 << "," << s+col+1 << ") = " << (resTemp(i)-res(i))/(yTemp(s+col)-y(s+col)) << ";" << endl;
		}
	}

	// J += c_j*dF/dydot
	// Banded, upper & lower bandwidths = nVars-1. 
	bw = 2*nVars-1;
	for (int s=0; s<bw; s++)
	{
		for (int i=0; i<N; i++) {
			yTemp(i) = y(i);
		}
		// Form the perturbed ydot vector
		for (int i=0; i<N; i++) {
			ydotTemp(i) = ydot(i);
		}
		for (int i=s; i<N; i+=bw) {
			// Form the perturbed y vector
			if (ydot(i) != 0) {
				ydotTemp(i) *= 1+eps;
			} else {
				ydotTemp(i) = eps;
			}
		}

		// Evaluate the residuals
		f(t, y, ydotTemp, resTemp);

		
		// Fill in the corresponding Jacobian elements
		int iStart = max(0,s-nVars);
		int counter;
		int col;
		if (s >= nVars) {
			counter = 0;
			col = -bw;
			
		} else {
			counter = nVars - s;
			col = 0;
		}
		
		for (int i=iStart; i<N; i++) {
			if (counter++ % bw == 0) {
				col += bw;
			}
			if (col+s>=N) {
				continue;
			}
			(*bandedJacobian)(i,s+col) += c_j*(resTemp(i)-res(i))/(ydotTemp(s+col)-ydot(s+col));
		}
	}

	BandGBTRF(bandedJacobian->forSundials(),&pMat[0]);
	
 	return 0;
}

int strainedFlameSys::preconditionerSolve(realtype t, sdVector& yIn, sdVector& ydotIn, 
										  sdVector& resIn, sdVector& rhs, 
										  sdVector& outVec, realtype c_j, realtype delta)
{
	//if (c_j != cjSave) {
	//	cout << "Warning: c_j changed unexpectedly!" << endl;
	//}
	vector<double> xVec(N);
	for (int i=0; i<N; i++) {
		xVec[i] = rhs(i);
	}

//	denGETRS(jacMatrix->forSundials(),N,&pMat[0],&xVec[0]);
	BandGBTRS(bandedJacobian->forSundials(),&pMat[0],&xVec[0]);


	for (int i=0; i<N; i++) {
		outVec(i) = xVec[i];
	}

	//// Diagonal preconditioner (?)
	//for (int i=0; i<N; i++) {
	//	outVec(i) = rhs(i)/(dFdy[i][i] + c_j*dFdydot[i][i]);
	//}
	return 0;
}


strainedFlameSys::strainedFlameSys(void) 
{
}

strainedFlameSys::~strainedFlameSys(void)
{
	delete bandedJacobian;
}

void strainedFlameSys::unrollYdot(const sdVector& yDot)
{
	for (int i=0; i<nPoints; i++) {
		drhovdt[i] = yDot(3*i);
		dUdt[i] = yDot(3*i+1);
		dTdt[i] = yDot(3*i+2);
	}
}

void strainedFlameSys::rollYdot(sdVector& yDot)
{
	for (int i=0; i<nPoints; i++) {
		yDot(3*i) = drhovdt[i];
		yDot(3*i+1) = dUdt[i];
		yDot(3*i+2) = dTdt[i];
	}
}

void strainedFlameSys::unrollY(const sdVector& y)
{
	for (int i=0; i<nPoints; i++) {
		rhov[i] = y(3*i);
		U[i] = y(3*i+1);
		T[i] = y(3*i+2);
	}
}

void strainedFlameSys::rollY(sdVector& y)
{
	for (int i=0; i<nPoints; i++) {
		y(3*i) = rhov[i];
		y(3*i+1) = U[i];
		y(3*i+2) = T[i];
	}
}

void strainedFlameSys::rollResiduals(sdVector& res)
{
	for (int i=0; i<nPoints; i++) {
		res(3*i) = resContinuity[i];
		res(3*i+1) = resMomentum[i];
		res(3*i+2) = resEnergy[i];
	}
}

void strainedFlameSys::printForMatlab(ofstream& file, vector<double>& v, int index, char* name)
{
	file << name << "(:," << index << ") = [";
	for (unsigned int i=0; i<v.size()-1; i++)
	{
		file << v[i] << ", ";
	}
	file << v[v.size()-1] << "];" << endl;
}