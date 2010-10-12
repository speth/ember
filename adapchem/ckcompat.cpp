#include "ckcompat.h"
#include <iostream>

ChemkinGas::ChemkinGas(const std::string& filename, bool quiet)
{
    int input_unit = 555;
    f90_open_unformatted_(&input_unit, const_cast<char*>(filename.c_str()), filename.length());

    int output_unit;
    if (quiet) {
        output_unit = 666;
        std::string output = "/dev/null";
        f90_open_formatted_(&output_unit, const_cast<char*>(output.c_str()), output.length());
    } else {
        output_unit = 6; // STDOUT
    }

    int len_iwork;
    int len_rwork;
    int len_cwork;
    int iflag;
    cklen_(&input_unit, &output_unit, &len_iwork, &len_rwork, &len_cwork, &iflag);
    std::cout << len_iwork << "," << len_rwork << "," << len_cwork << std::endl;

    _rwork.resize(len_rwork);
    _iwork.resize(len_iwork);
    _cwork.reserve(len_cwork*16);
    f90_ckinit_(&input_unit, &output_unit,
              &len_iwork, &_iwork[0], &len_rwork, &_rwork[0],
              &len_cwork, &_cwork[0], len_cwork);

    f90_close_(&input_unit);

    _iptr = &_iwork[0];
    _rptr = &_rwork[0];

    int mm, ii, nfit;
    ckindx_(_iptr, _rptr, &mm, &_nSpec, &ii, &nfit);
    _Y.resize(_nSpec);
}

void ChemkinGas::setStatePTY(double P, double T, const double* const Y)
{
    #ifdef CKCOMPAT_USE_MKS
        _P = P*10;
    #else
        _P = P;
    #endif

    _T = T;
    for (int k=0; k<_nSpec; k++) {
        _Y[k] = Y[k];
    }
}

void ChemkinGas::setStateMass(const double* const Y, double T)
{
    _T = T;
    for (int k=0; k<_nSpec; k++) {
        _Y[k] = Y[k];
    }
}

void ChemkinGas::setPressure(double P)
{
    #ifdef CKCOMPAT_USE_MKS
        _P = P*10;
    #else
        _P = P;
    #endif
}


void ChemkinGas::getReactionRates(double* wdot)
{
    ckwyp_(&_P, &_T, &_Y[0], _iptr, _rptr, wdot);

    #ifdef CKCOMPAT_USE_MKS
        for (int k=0; k<_nSpec; k++) {
            wdot[k] *= 1000;
        }
    #endif
}



AdapChem::AdapChem(const std::string& filename, bool quiet)
{
    int input_unit = 555;
    f90_open_unformatted_(&input_unit, const_cast<char*>(filename.c_str()), filename.length());

    int output_unit;
    if (quiet) {
        output_unit = 666;
        std::string output = "/dev/null";
        f90_open_formatted_(&output_unit, const_cast<char*>(output.c_str()), output.length());
    } else {
        output_unit = 6; // STDOUT
    }

    int len_iwork;
    int len_rwork;
    int len_cwork;
    int iflag;
    cklen_(&input_unit, &output_unit, &len_iwork, &len_rwork, &len_cwork, &iflag);
    std::cout << len_iwork << "," << len_rwork << "," << len_cwork << std::endl;

    _rwork.resize(len_rwork);
    _iwork.resize(len_iwork);
    _cwork.reserve(len_cwork*16);
    f90_ckinit_(&input_unit, &output_unit,
              &len_iwork, &_iwork[0], &len_rwork, &_rwork[0],
              &len_cwork, &_cwork[0], len_cwork);

    f90_close_(&input_unit);

    _iptr = &_iwork[0];
    _rptr = &_rwork[0];

    int mm, ii, nfit;
    ckindx_(_iptr, _rptr, &mm, &_nSpec, &ii, &nfit);
    _Y.resize(_nSpec);

    _wdot.resize(_nSpec);

    adapchemini_();
    adpchm2_.icount = 0;
}

AdapChem::~AdapChem()
{
    adapchemend_();
}

double AdapChem::getRelTol() const
{
    return adpchm1_.resid_tol;
}

void AdapChem::incrementStep()
{
    adpchm2_.icount += 1;
}

void AdapChem::setGridSize(int nGrid) {
    _nGrid = nGrid;
    _nGridm1 = _nGrid - 1;

    _iSpecies.resize(_nGrid);
    _keepSpecies.resize(_nSpec);
}

void AdapChem::initializeStep(double P, double T, const double* const Y, int iGrid)
{
    _iGrid = iGrid;
    setStatePTY(P, T, Y, iGrid);
    callAdapChem(-1);
    callAdapChem(0);

    // Store the indices of species that are being kept at this point
    _iSpecies[_iGrid].clear();
    for (int i=0; i<_nSpec; i++) {
        if (_keepSpecies[i]) {
            _iSpecies[iGrid].push_back(i);
        }
    }
}

void AdapChem::initializeStep(const double* const Y, double T, int iGrid)
{
    _iGrid = iGrid;
    setStateMass(Y, T, iGrid);
    callAdapChem(-1);
    callAdapChem(0);

    // Store the indices of species that are being kept at this point
    _iSpecies[_iGrid].clear();
    for (int i=0; i<_nSpec; i++) {
        if (_keepSpecies[i]) {
            _iSpecies[iGrid].push_back(i);
        }
    }
}


void AdapChem::setStatePTY(double P, double T, const double* const Y, int iGrid)
{
    #ifdef CKCOMPAT_USE_MKS
        _P = P*10;
    #else
        _P = P;
    #endif

    _T = T;
    for (int k=0; k<_nSpec; k++) {
        _Y[k] = Y[k];
    }

    _iGrid = iGrid;
}


void AdapChem::setStateMass(const double* const Y, double T, int iGrid)
{
    _T = T;
    for (int k=0; k<_nSpec; k++) {
        _Y[k] = Y[k];
    }

    _iGrid = iGrid;
}


void AdapChem::setPressure(double P)
{
    #ifdef CKCOMPAT_USE_MKS
        _P = P*10;
    #else
        _P = P;
    #endif
}


void AdapChem::getReactionRates(double* wdot)
{
    callAdapChem(1);

    const std::vector<size_t>& iSpec = _iSpecies[_iGrid];
    for (size_t i=0; i<iSpec.size(); i++) {
        #ifdef CKCOMPAT_USE_MKS
            wdot[iSpec[i]] = 1000.0*_wdot[i];
        #else
            wdot[iSpec[i]] = _wdot[i];
        #endif
    }

}

void AdapChem::callAdapChem(int jSolv)
{
    int igridTmp = _iGrid+1; // Fortran is 1-indexed...
    adapchem_(&igridTmp, &_nGrid, &_nGridm1,
              &_P, &_T, &_Y[0],
              _iptr, _rptr,
              &_wdot[0], &_keepSpecies[0], &jSolv);
}
