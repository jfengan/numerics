//
// Created by Jiahao FENG on 20/1/2021.
//

#ifndef NUMERICS_LINEAR_INTERP_H
#define NUMERICS_LINEAR_INTERP_H

#include "base_interp.h"
using namespace numerics;

struct Linear_interp: Base_interp{
    Linear_interp(std::vector<double> &xv, std::vector<double> &yv)
    :Base_interp(xv, &yv[0], 2) {}

    double rawinterp(int j, double x) override{
        if(xx[j] == xx[j+1])
            return yy[j];
        else
            return yy[j] + ((x - xx[j])/(xx[j+1] - xx[j])) * (yy[j+1] - yy[j]);
    }
};

#endif //NUMERICS_LINEAR_INTERP_H
