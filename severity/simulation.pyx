"""
Copyright (c) 2017 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from denovonear.weights cimport WeightedChoice, Chooser

cdef extern from "simulate.h":
    double _analyse(Chooser, vector[double], double, int, int) except +

def analyse(WeightedChoice choices, severity, observed, count,
        iterations=100000000):
    ''' analyse the severity score of de novo mutations in a gene
    
    estimate the chance of observing a total severity euqal to or greater than
    the summed scores for the observed de novos in a gene, given the same number
    of de novo mutations.
    
    Args:
        choices: WeightedChoice object of mutation rates per position and alt
            allele for all possible SNVs within a gene.
        severity: list of severity scores, matching the same position and alt
            allele order as for the choices object.
        observed: summed severity score across the observed de novo mutations.
        count: number of observed de novo mutations.
        iterations: number of iterations to run
    
    Returns:
        probability of getting the observed severity score (or greater) under
        the null distribution.
    '''
    
    return _analyse(deref(choices.thisptr), severity, observed, count, iterations)
