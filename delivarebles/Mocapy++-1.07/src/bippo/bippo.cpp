/*
 *  bippo.cpp
 *
 *  Copyright (C) 2008, Thomas Hamelryck, The Bioinformatics Centre,
 *  University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "bippo.h"

using namespace std;

namespace mocapy {

Bippo::Bippo(double lambda, double theta, uint n, uint cache_max) {
    // cache_max: lik/loglik values are cached up to cache_max
    assert(lambda>0);
    this->lambda = lambda;
    assert(theta>0 && theta<1);
    this->theta = theta;
    assert(n>0);
    this->n = n;
    assert(cache_max>0);
    this->cache_max = cache_max;

    // cache loglik & lik
    cache();
}

void Bippo::cache(void) {
    // initialize cache for loglik & lik
    lik_cache.resize(cache_max);
    loglik_cache.resize(cache_max);

    // cache the likelihood and loglik values
    for(uint i=0; i<cache_max; i++)
    {
        double ll;
        ll = calc_loglik(i);
        loglik_cache[i] = ll;
        lik_cache[i] = exp(ll);
    }
}

void Bippo::estimate(double m, double m2, double m3) {
    double th, nn, la;
    double max_nn;
    double max_delta = 1e20;

    // some heuristic boundaries for n
    max_nn = max(5*m, 5.1);

    // loop over n and theta
    // find minimum difference with mean and variance
    for(th=0.01; th<0.91; th+=0.01)
    {
        for(la=0.2; la<10; la+=0.01)
        {
            for(nn=2; nn<max_nn; nn++)
            {
                double delta, e, e2;

                // theoretical moments
                e  = la+nn*th;
                e2 = la+nn*th*(1-th);

                // diff with real moments
                delta = (m-e)*(m-e) + (m2-e2)*(m2-e2);

                if(delta<max_delta)
                {
                    // another optimum - update parameters
                    max_delta = delta;
                    n = (uint) nn+0.5;
                    theta = th;
                    lambda = la;
                }
            }
        }
    }
}

double Bippo::calc_loglik(uint c) {
    // calculate log likelihood
    uint maxval;
    double sum=0;

    if (c>n) maxval=n;
    else maxval=c;

    for(uint i=0; i<=maxval; i++)
    {
        double a,b;
        a = pow(lambda, (double) c-i) * pow(theta, (double) i) * pow(1-theta, (double)n-i);
        b = exp(log_fac.get(c-i)) * exp(log_fac.get(i)) * exp(log_fac.get(n-i));
        sum += a/b;
    }
    
    double ll=-lambda + log_fac.get(n) + log(sum);

    return ll;
}

vector<double> Bippo::get_parameters(){
    vector<double> params(3);

    params[0]=lambda;
    params[1]=theta;
    params[2]=(double) n;

    return params;
}

double Bippo::sample(RandomGen* rg) {
    // First sample from Poisson...
    double sample = poisson_sample(rg);
    // ...then from binomial (n tries)
    for (uint i=0; i<n; i++)
    {
        double r = rg->get_rand();

        if (r<theta) sample++;
    }

    return sample;
}

double Bippo::poisson_sample(RandomGen* rg)
{
    // sample from the poisson distribution
    // use lambda from Bippo
    double l=exp(-lambda);
    double k=0, p=1;

    while(1)
    {
        k++;
        double u = rg->get_rand();
        p *= u;
        if (p<l) break;
    }

    return k-1;
}

}
