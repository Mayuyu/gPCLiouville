//
//  CollocationSolver.cpp
//  gPCLiouville
//
//  Created by 马 征 on 14-6-10.
//  Copyright (c) 2014年 mz. All rights reserved.
//

#include "CollocationSolver.h"

Eigen::MatrixXd CollocationSolver::spatialDiff() {
    auto tmp = mesh;
    for (int i = 1; i < mesh.rows() - 1; i++) {
        for (int j = 1; j < mesh.cols() - 1; j++) {
            double x = getX(i), v = getV(j);
            double gradV = (rightPotential(x - 0.5 * deltaX, v) - leftPotential(x + 0.5 * deltaX, v)) / deltaX;
            tmp(i, j) = 0.5 * (-(v + fabs(v)) * positiveXFluxDiff(i, j) - (v - fabs(v)) * negativeXFluxDiff(i, j) +
                        (gradV + fabs(gradV)) * positiveVfluxDiff(i, j) + (gradV - fabs(gradV)) * negativeVfluxDiff(i, j));
        }
    }
    return tmp;
}

void CollocationSolver::solveEquation() {
    for (int t = 0; t < timeSteps; t++) {
//        auto tmp0 = mesh;
//        auto tmp1 = spatialDiff();
//        mesh += deltaT * tmp1;
//        mesh = tmp0 + 0.5 * deltaT * (tmp1 + spatialDiff());
        mesh += deltaT * spatialDiff();
    }
}

void CollocationSolver::initial(TypeFunction initial) {
    for (int i = 0; i < mesh.rows(); i++) {
        for (int j = 0; j < mesh.cols(); j++) {
            mesh(i, j) = initial(getX(i), getV(j));
        }
    }
}

Eigen::MatrixXd CollocationSolver::exactSolution(TypeFunction exactSolution) {
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(mesh.rows(), mesh.cols());
    for (int i = 0; i < mesh.rows(); i++) {
        for (int j = 0; j < mesh.cols(); j++) {
            temp(i, j) = exactSolution(getX(i), getV(j));
        }
    }
    return temp;
}