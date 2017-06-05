
#ifndef SEVERITY_SIMULATE_H_
#define SEVERITY_SIMULATE_H_

#include <vector>
#include <map>

#include "weighted_choice.h"

struct Site {
    int pos, offset;
    std::string alt;
    
    bool operator=(const Site &other) const {
        return pos == other.pos and offset == other.offset and alt == other.alt;
    }
    
    bool operator<(const Site &other) const {
        return pos < other.pos or
            (pos == other.pos and offset < other.offset) or
            (pos == other.pos and offset == other.offset and alt < other.alt);
    }
};

std::map<Site, double> prepare_severity(Chooser &choices, std::vector<double> &severity);
double _analyse(Chooser &choices, std::vector<double> severity, double observed,
    int count, int iterations);

#endif // SEVERITY_SIMULATE_H_
