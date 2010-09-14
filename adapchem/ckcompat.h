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
    void setPressure(double P);
    void setStateMass(const double* const Y, double T);

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
    void incrementStep();

    void setStatePTY(double P, double T, const double* const Y);
    void wdot(double* wdot);


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
