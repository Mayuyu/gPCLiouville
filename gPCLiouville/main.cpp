//
//  main.cpp
//  gPCLiouville
//
//  Created by 马 征 on 14-6-9.
//  Copyright (c) 2014年 mz. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "CollocationSolver.h"

using namespace std;

double leftPotential(double x, double v) {
    if (x > EPS) {
        return 0.0;
    } else {
        return 0.2;
    }
}

double rightPotential(double x, double v) {
    if (x < -EPS) {
        return 0.2;
    } else {
        return 0.0;
    }
}

double initial(double x, double v) {
    if (x * x + v * v < 1.0 && x >= 0.0 && v < 0.0) {
        return 1.0;
    } else if (x * x + v * v < 1.0 && x <= 0.0 && v > 0.0) {
        return 1.0;
    } else {
        return 0.0;
    }
}

double exactSolution(double x, double v) {
    if (x >= 0 && v < sqrt(0.4) && x < v) {
        return 1.0;
    } else if (x >= 0 && x < 1 && v < 0 && v >0.5 * (x - sqrt(2.0 - x * x)) ) {
        return 1.0;
    } else if (x <= 0 && v < x && v > -sqrt(0.6) && x < (1.0 - sqrt(0.6 - v * v) / sqrt(v * v + 0.4)) * v) {
        return 1.0;
    } else if (x <= 0 && v > 0 && x > -1 && v < 0.5 * (x + sqrt(2.0 - x * x))) {
        return 1.0;
    } else if (x >= 0 && v > sqrt(0.4) && v > x && v < sqrt(1.4) && x > (1.0 - sqrt(1.4 - v * v) / sqrt(v * v - 0.4)) * v) {
        return 1.0;
    } else {
        return 0.0;
    }
}

int main(int argc, const char *argv[]) {
    const double a = 0.03, b = 0.02;
    TypeFunction leftP = leftPotential, rightP = rightPotential, init = initial,
    exactSolu = exactSolution;
    CollocationSolver test(1.5 + 0.5 * a, a, b, 1.0, leftP, rightP);
    test.initial(init);
    test.solveEquation();
    cout << (test.getMesh() - test.exactSolution(exactSolu)).lpNorm<1>() / test.getMesh().size() << endl;
    return 0;
}

