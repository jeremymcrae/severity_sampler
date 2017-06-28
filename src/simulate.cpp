// Copyright (c) 2017 Genome Research Ltd.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
// IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <vector>
#include <algorithm>
#include <map>
#include <stdexcept>

#include "simulate.h"

std::map<int, std::map<std::string, double>> prepare_severity(Chooser &choices,
    std::vector<double> &severity) {
    /**
        prepare severity scores, to map from site and alt to a given score
    */
    
    std::map<int, std::map<std::string, double>> sites;
    
    for (int i=0; i < choices.len(); i++) {
        auto x = choices.iter(i);
        sites[x.pos][x.alt] = severity[i];
    }
    
    return sites;
}

double _analyse(Chooser &choices, std::vector<double> severity, double observed,
    int count, int iterations) {
    /**
        simulates the probability of observing n de novos with a combined
        severity score greater than or equal to the obsevered severity total.
    
    */
    int len = choices.len();
    int sev_len = severity.size();
    if (len != sev_len) { throw std::invalid_argument("severity scores do not match rates!"); }
    if (len == 0) { throw std::invalid_argument("no per-base/allele rates supplied!"); }
    if (count == 0) { throw std::invalid_argument("sampling zero de novos!"); }
    
    // figure out how to map sites to severity scores. This requires at a given
    // index position the data within the choices object and the severity vector
    // are for the same site/alt allele.
    auto scores = prepare_severity(choices, severity);
    
    double minimum_p = 1.0/(1.0 + static_cast<double>(iterations));
    double p_value = minimum_p;
    std::vector<double> dist;
    
    while (iterations < 100000000) {
        
        std::vector<double> current;
        int increment = iterations - dist.size();
        for (int n=0; n < increment; n++) {
            double total = 0.0;
            for (int i=0; i < count; i++) {
                auto x = choices.choice();
                total += scores[x.pos][x.alt];
            }
            
            current.push_back(total);
        }
        
        // merge the two sorted lists into a sorted vector
        std::vector<double> v(iterations);
        std::merge(dist.begin(), dist.end(), current.begin(), current.end(),
            v.begin());
        dist = v;
        
        // make sure the scores are sorted, so we can quickly find the position
        std::sort(dist.begin(), dist.end());
        
        // figure out where in the list a random probability would fall
        std::vector<double>::iterator pos;
        pos = std::upper_bound(dist.begin(), dist.end(), observed);
        double position = pos - dist.begin();
        
        // estimate the probability from the position
        p_value = (1.0 + (dist.size() - position))/(1.0 + dist.size());
        
        double z = 2.575829;
        double precision = 0.05;
        if (_halt_permutation(p_value, iterations, z, precision)) { break; }
        
        iterations += 10000000;  // for if we need to run more iterations
    }
    
    return p_value;
}

bool _halt_permutation(double p_val, int iterations, double z, double precision) {
    /**
        halt permutations if the P value is sufficiently precise
        
        given a p-value, and a number of iterations, we can figure out a
        confidence interval around the p-value (assuming a normal approximation
        to the binomal). We calculate the CI delta (the value added or
        subtracted). Once this is less than 5% of the estimated p-value, we
        consider the p-value to be sufficiently precise, and halt the
        permutation testing.
        
        @p_val current simulated P value
        @iterations iterations run in order to obtain the simulated P value
        @z standard normal deviate (eg 1.96 for 95% CI)
        @precision value above which to cease permuting
        @return whether to halt the permuations
    */
    double delta = z * sqrt((p_val * (1 - p_val))/iterations);
    double diff = delta/p_val;
    
    // if the 'variance' of the confidence interval is lower than 5%, we can be
    // confident that the p-value won't change much and can halt.
    return diff < precision;
}
