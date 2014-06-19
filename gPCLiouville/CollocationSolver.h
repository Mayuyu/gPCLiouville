//
//  CollocationSolver.h
//  gPCLiouville
//
//  Created by 马 征 on 14-6-10.
//  Copyright (c) 2014年 mz. All rights reserved.
//

#ifndef __gPCLiouville__CollocationSolver__
#define __gPCLiouville__CollocationSolver__

#include <iostream>
#include "AbsLiouvilleSolver.h"

class CollocationSolver : public AbsLiouvilleSolver<double> {
public:
    CollocationSolver(double xLeft, double xRight,
                      double vBottom, double vTop,
                      double deltaX, double deltaV,
                      double deltaT, double finalTime,
                      TypeFunction leftPotential, TypeFunction rightPotential)
    : AbsLiouvilleSolver(xLeft, vBottom, deltaX, deltaV, deltaT, finalTime,
                              leftPotential, rightPotential) {
        int rows = (xRight - xLeft) / deltaX + 1;
        int cols = (vTop - vBottom) / deltaV + 1;
        mesh = Eigen::MatrixXd::Zero(rows, cols);
    }
    CollocationSolver(double bound, double deltaSpatial,
                      double deltaT, double finalTime,
                      TypeFunction leftPotential, TypeFunction rightPotential)
    : AbsLiouvilleSolver(-bound, -bound, deltaSpatial, deltaSpatial, deltaT,
                              finalTime, leftPotential, rightPotential) {
        int dim = 2 * bound / deltaSpatial + 1;
        mesh = Eigen::MatrixXd::Zero(dim, dim);
    }
    void initial(TypeFunction initial);
    void print() { std::cout << mesh << std::endl; }
    Eigen::MatrixXd getMesh() { return mesh; }
    Eigen::MatrixXd exactSolution(TypeFunction exactSolution);
    Eigen::MatrixXd spatialDiff();
    virtual void solveEquation();
};


#endif /* defined(__gPCLiouville__CollocationSolver__) */
