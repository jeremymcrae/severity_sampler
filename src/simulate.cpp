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

std::map<Site, double> prepare_severity(Chooser &choices, std::vector<double> &severity) {
    /**
        prepare severity scores, to map from site and alt to a given score
    */
    
    std::map<Site, double> sites;
    
    for (int i=0; i < choices.len(); i++) {
        auto x = choices.iter(i);
        sites[Site {x.pos, x.offset, x.alt}] = severity[i];
    }
    
    return sites;
}

double _analyse(Chooser &choices, std::vector<double> severity, double observed,
    int count, int iterations) {
    /**
        simulates the probability of observing n de novos with a combined
        severity score greater than or equal to the obsevered severity total.
    
    */
    
    if (choices.len() != (int)severity.size()) {
        throw std::invalid_argument("severity scores do not match rates!");
    }
    
    if (choices.len() == 0) {
        throw std::invalid_argument("no per-base/allele rates supplied!");
    }
    
    // figure out how to map sites to severity scores. This requires at a given
    // index position the data within the choices object and the severity vector
    // are for the same site/alt allele.
    std::map<Site, double> scores = prepare_severity(choices, severity);
    
    std::vector<double> dist;
    
    for (int n=0; n < iterations; n++) {
        double total = 0.0;
        for (int i=0; i < count; i++) {
            auto x = choices.choice();
            total += scores[Site {x.pos, x.offset, x.alt}];
        }
        dist.push_back(total);
    }
    
    // make sure the scores are sorted, so we can quickly find the position
    std::sort(dist.begin(), dist.end());
    
    // figure out where in the list a random probability would fall
    std::vector<double>::iterator pos;
    pos = std::upper_bound(dist.begin(), dist.end(), observed);
    double position = pos - dist.begin();
    
    // estimate the probability from the position
    return (1.0 + (dist.size() - position))/(1.0 + dist.size());
}
