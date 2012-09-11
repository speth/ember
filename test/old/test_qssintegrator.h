#pragma once

#include "../src/qssintegrator.h"

class QSSTestProblem : public QssIntegrator
{
public:
    QSSTestProblem() {}
    void odefun(double t, const dvector& y, dvector& q, dvector& d);
};
