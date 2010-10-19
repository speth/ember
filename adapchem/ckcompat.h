#pragma once

// Work in MKS units rather than Chemkin's default CGS units
#define CKCOMPAT_USE_MKS

#include <string>
#include <vector>

extern "C"
{
    // Functions defined in helper.f90
    extern void f90_open_unformatted_(int* unit, char* name, int len_name);
    extern void f90_open_formatted_(int* unit, char* name, int len_name);
    extern void f90_close_(int* unit);

    extern void f90_ckinit_(int* in_unit, int* out_unit, int* lenick, int* ickwrk, int* lenrwk, double* rckwrk,
                           int* lencwk, char* cckwrk, int lencck2);

    extern void f90_setrtol_(double* rtol);

    // Functions defined by Chemkin:
    extern void cklen_(int* linc, int* lout, int* li, int* lr, int* lc, int* iflag);
    extern void ckindx_(int* iwork, double* rwork, int* mm, int* kk, int* ii, int* nfit);
    extern void ckwyp_(double* P, double* T, double* Y, int* iwork, double* rwork, double* wdot);

    // Functions defined by AdapChem:
    extern void adapchemini_(void);
    extern void adapchemend_(void);
    extern void adapchem_(int* igrid, int* igridTot, int* igm1,
                          double* P, double* T, double* Y,
                          int* iwork, double* rwork, double* wdot,
                          int* ispecies, int* jsolv);

    extern void adapbox_(int* igridTot, int* nmaxbox, double* Y,
                         double* P, double* T, double* rwork, double* iwork);

    // COMMON blocks defined by AdapChem:
    extern struct {
        double resid_tol;
    } adpchm1_;

    extern struct {
        int icount;
    } adpchm2_;
}


class ChemkinGas
{
public:
    ChemkinGas(const std::string& filename, bool quiet=false);

    void setStatePTY(double P, double T, const double* const Y);
    void setStateMass(const double* const Y, double T);
    void setPressure(double P);

    void getReactionRates(double* wdot);

private:
    // Chemkin's work arrays
    std::vector<double> _rwork;
    std::vector<int> _iwork;
    std::vector<char> _cwork;

    // Pointers declared for convenience
    double* _rptr;
    int* _iptr;

    int _nSpec; // number of species

    // current state: pressure, temperature and mass fraction
    double _P;
    double _T;
    std::vector<double> _Y;
};


class AdapChem
{
public:
    AdapChem(const std::string& filename, bool quiet=false);
    ~AdapChem();

    double getRelTol() const;

    // *** Initialize a new step ***

    // call to mark the beginning of the new step
    void incrementStep();

    // call if grid size changed
    void setGridSize(int nGrid);

    // call (one of these) once for each point to set the initial condition
    void initializeStep(double P, double T, const double* const Y, int iGrid);
    void initializeStep(const double* const Y, double T, int iGrid);

    // *** Functions called during integration
    void setStatePTY(double P, double T, const double* const Y, int iGrid);
    void setStateMass(const double* const Y, double T, int iGrid);
    void setPressure(double P);
    void getReactionRates(double* wdot);

    // Functions for informational purposes:
    size_t getNumSpec(size_t iGrid); // number of species in reduced model at grid point iGrid

private:
    // Chemkin's work arrays
    std::vector<double> _rwork;
    std::vector<int> _iwork;
    std::vector<char> _cwork;

    // Pointers declared for convenience
    double* _rptr;
    int* _iptr;

    int _nSpec; // number of species
    int _nGrid; // number of grid points
    int _nGridm1; // = _nGrid - 1

    std::vector<int> _keepSpecies; // binary vector returned by adapchem (temporary)

    // vector indicating indices of species present in reduced model at each point
    std::vector<std::vector<size_t> > _iSpecies;

    // production rate, but only for species included in model
    std::vector<double> _wdot;

    // current state: pressure, temperature, mass fraction, and grid point index
    double _P;
    double _T;
    std::vector<double> _Y;
    int _iGrid;

    void callAdapChem(int jSolv);
};
