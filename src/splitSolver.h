#pragma once

#include "mathUtils.h"

class SplitSolver
{
public:
    SplitSolver() {}
    virtual ~SplitSolver() {}

    void resize(index_t nRows, index_t nCols) {
        deltaConv.resize(nRows, nCols);
        deltaDiff.resize(nRows, nCols);
        deltaProd.resize(nRows, nCols);
        ddtConv.resize(nRows, nCols);
        ddtDiff.resize(nRows, nCols);
        ddtProd.resize(nRows, nCols);
        ddtCross.resize(nRows, nCols);
        Sstart.resize(nRows, nCols);
        state.resize(nRows, nCols);
    }

    // Current state
    dmatrix state;

    // State at the start of the current integrator stage
    dmatrix Sstart;

    // Changes in each state variable for each terms of the governing equations
    dmatrix deltaConv, deltaDiff, deltaProd;

    // Estimated time derivatives for each split term
    dmatrix ddtConv, ddtDiff, ddtProd, ddtCross;
};
