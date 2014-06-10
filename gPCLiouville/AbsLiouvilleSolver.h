//
//  AbsLiouvilleSolver.h
//  gPCLiouville
//
//  Created by 马 征 on 14-6-9.
//  Copyright (c) 2014年 mz. All rights reserved.
//

#ifndef gPCLiouville_AbsLiouvilleSolver_h
#define gPCLiouville_AbsLiouvilleSolver_h

#include "Eigen/Dense"

const double EPS = 1e-8;
typedef double (*TypeFunction)(double x, double v);

template <class T>
class AbsLiouvilleSolver {
protected:
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mesh;
    double xLeft;
    double vBottom;
    double deltaX;
    double deltaV;
    double deltaT;
    int timeSteps;
    TypeFunction leftPotential;
    TypeFunction rightPotential;
public:
    double getX(int i) { return xLeft + i * deltaX; }
    double getV(int j) { return vBottom + j * deltaV; }
    int getMeshIndexOfV(double v) {
        for (int j = 0; j < mesh.cols(); j++) {
            if (getV(j) - v <= EPS && getV(j + 1) - v > EPS) {
                return j;
            }
        }
        return -1;
    }
    T positiveXFluxDiff(int i, int j);
    T negativeXFluxDiff(int i, int j);
    T positiveVfluxDiff(int i, int j);
    T negativeVfluxDiff(int i, int j);
    virtual void solveEquation() = 0;
};
template <class T>
T AbsLiouvilleSolver<T>::positiveXFluxDiff(int i, int j) {
    T uMinus = mesh(i, j), uPlus;
    auto x = getX(i), v = getV(j);
    auto potentialDiff = leftPotential(x - 0.5 * deltaX, v) - rightPotential(x - 0.5 * deltaX, v);
    if (v * v - 2.0 * potentialDiff < -EPS) {
        uPlus = mesh(i, getMeshIndexOfV(-v));
    } else {
        auto vLeft = sqrt(v * v - 2.0 * potentialDiff);
        int k = getMeshIndexOfV(vLeft);
        if (k == -1) {
            uPlus = 0.0;
        } else {
            uPlus = (getV(k + 1) - vLeft) * mesh(i - 1, k) / deltaV + (vLeft - getV(k)) * mesh(i - 1, k + 1) / deltaV;
        }
    }
    return (uMinus - uPlus) / deltaX;
}
template <class T>
T AbsLiouvilleSolver<T>::negativeXFluxDiff(int i, int j) {
    T uMinus, uPlus = mesh(i, j);
    double x = getX(i), v = getV(j);
    double potentialDiff = leftPotential(x + 0.5 * deltaX, v) - rightPotential(x + 0.5 * deltaX, v);
    if (v * v + 2.0 * potentialDiff < -EPS) {
        uMinus = mesh(i, getMeshIndexOfV(-v));
    } else {
        auto vRight = -sqrt(v * v + 2.0 * potentialDiff);
        int k = getMeshIndexOfV(vRight);
        if (k == -1) {
            uMinus = 0.0;
        } else {
            uMinus = (getV(k + 1) - vRight) * mesh(i + 1, k) / deltaV + (vRight - getV(k)) * mesh(i + 1, k + 1) / deltaV;
        }
    }
    return (uMinus - uPlus) / deltaX;
}
template <class T>
T AbsLiouvilleSolver<T>::positiveVfluxDiff(int i, int j) {
    return (mesh(i, j) - mesh(i, j - 1)) / deltaV;
}
template <class T>
T AbsLiouvilleSolver<T>::negativeVfluxDiff(int i, int j) {
    return (mesh(i, j + 1) - mesh(i, j)) / deltaV;
}
#endif
