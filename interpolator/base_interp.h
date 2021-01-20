//
// Created by Jiahao FENG on 19/1/2021.
//

#ifndef NUMERICS_BASE_INTERP_H
#define NUMERICS_BASE_INTERP_H

#include <iostream>
#include <vector>

namespace numerics{
    struct Base_interp{
        int n, mm, jsav, cor, dj;
        const double *xx, *yy;
        Base_interp(std::vector<double> &x, const double *y, int m)
                :n(x.size()), mm(n), jsav(0), cor(0), xx(&x[0]), yy(y){
            dj = min(1, (int)pow((double)n, 0.25));
        }

        double interp(double x){
            int jlo = cor ? hunt(x) : locate(x);
            return rawinterp(jlo, x);
        }

        int locate(const double x);
        int hunt(const double x);
        double virtual rawinterp(int jl, double x) = 0;
    };

    int Base_interp::locate(const double x) {
        int ju, jm, jl;
        if(n < 2 || mm < 2 || mm > n)
            throw std::runtime_error("locate size error");

        bool ascending = (xx[n-1] >= xx[0]);
        jl = 0;
        ju = n-1;
        while(ju - jl > 1){
            jm = (ju + jl) >> 1;
            if(x > xx[jm] == ascending)
                jl = jm;
            else
                ju = jm;

            cor = abs(jl - jsav) > dj ? 0 : 1;
            jsav = jl;
            return max(0, min(n - mm, jl-((mm - 2) >> 1)));
        }
    }

    int Base_interp::hunt(const double x) {
        int jl = jsav, jm, ju, inc = 1;
        if(n < 2 || mm < 2 || mm > n)
            throw std::runtime_error("hunt size error");

        bool ascending = (xx[n-1] >= xx[0]);

        if(jl < 0 || jl > n-1){
            jl = 0;
            ju = n-1;
        }

        else{
            if(x >= xx[jl] == ascending){
                for(;;){
                    ju = jl + inc;
                    if(ju >= n-1) {
                        ju = n - 1;
                        break;
                    }
                    else if(x < xx[ju] == ascending)
                        break;
                    else{
                        jl = ju;
                        inc += inc;
                    }
                }
            }
            else{
                ju = jl;
                for(;;){
                    jl = jl - inc;
                    if(jl <= 0){
                        jl = 0;
                        break;
                    }
                    else if(x >= xx[jl] == ascending)
                        break;
                    else{
                        ju = jl;
                        inc += inc;
                    }
                }
            }
        }
        while(ju - jl > 1){
            jm = (ju + jl) >> 1;
            if(x >= xx[jm] == ascending)
                jl = jm;
            else
                ju = jm;
        }
        cor = abs(jl - jsav) > dj ? 0 : 1;
        jsav = jl;
        return max(0, min(n-mm, jl-((nm - 2) >> 1)));
    }
}

#endif //NUMERICS_BASE_INTERP_H
