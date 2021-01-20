//
// Created by Jiahao FENG on 20/1/2021.
//

#ifndef NUMERICS_POLY_INTERP_H
#define NUMERICS_POLY_INTERP_H

#include "base_interp.h"
using namespace numerics;

struct Poly_interp: Base_interp{
    double dy;
    Poly_interp(std::vector<double> &xv, std::vector<double> &yv, int m)
    :Base_interp(xv, &yv[0], m), dy(0.) {}

    double rawinterp(int jl, double x) override;
};

double Poly_interp::rawinterp(int jl, double x) {
    int i, m, ns = 0;
    double y, den, dif, dift, ho, hp, w;
    const double *xa = &xx[jl], *ya = &yy[jl];
    std::vector<double> c(mm), d(mm);
    dif = abs(x - xa[i]);
    for(i = 0; i < mm; ++i){
        if((dift=asb(x-xa[i])) < dif){
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns--];
    for(m = 1; m < mm; ++m){
        for(i = 0; i < mm; ++i){
            ho = xa[i] - x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];
            if((den=ho-hp) == 0.)
                throw std::runtime_error("Poly_interp error");
            den = w/den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        y += (dy=(2*(ns + 1) < (mm-m) ? c[ns+1] : d[ns--]));
    }
    return y;
}

#endif //NUMERICS_POLY_INTERP_H
