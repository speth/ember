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

protected:
    dmatrix deltaConv, deltaDiff, deltaProd;
    dmatrix ddtConv, ddtDiff, ddtProd, ddtCross;
    dmatrix Sstart, state;
};
