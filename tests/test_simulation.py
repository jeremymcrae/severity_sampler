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

import unittest
import tempfile
from random import randint, uniform, seed

from denovonear.weights import WeightedChoice
from severity.simulation import analyse

class TestSimulationsPy(unittest.TestCase):
    ''' unit test functions for simulations
    '''
    
    def test_analyse(self):
        ''' test that we run the simulations correctly
        '''
        
        rates = WeightedChoice()
        rates.add_choice(200, 1e-5, 'A', 'G')
        rates.add_choice(201, 2e-5, 'C', 'T')
        rates.add_choice(202, 1e-5, 'C', 'G')
        
        severity = [5, 10, 5]
        
        # define a test where the observed score will fall at the midpoint of
        # the simulated null distribution
        p = analyse(rates, severity, 8, 1, iterations=100000)
        self.assertAlmostEqual(p, 0.5, places=2)
        
        # now check when we sample two de novo mutations
        p = analyse(rates, severity, 15, 2, iterations=100000)
        self.assertAlmostEqual(p, 0.25, places=2)
    
    def test_analyse_extreme_p_value(self):
        ''' test when the observed severity score exceeds all possible values
        '''
        
        rates = WeightedChoice()
        rates.add_choice(200, 1e-5, 'A', 'G')
        rates.add_choice(201, 2e-5, 'C', 'T')
        rates.add_choice(202, 1e-5, 'C', 'G')
        
        severity = [5, 10, 5]
        
        # now check when the observed severity score exceeds all possible
        # values from the severity distribution. This test gives an absurd
        # p-value at 1e-6, but that is because the observed value is
        # unachievable given the existsing severity scores. In practice the
        # observed score will always be theoretically achieveable in the null
        # distribution, since the observed score is calculated from the
        # existsing scores.
        p = analyse(rates, severity, 20, 1, iterations=100000)
        self.assertAlmostEqual(p, 1e-6, places=4)
    
    def test_analyse_empty(self):
        ''' check we raise an error if the rates and severity are empty
        '''
        
        with self.assertRaises(ValueError):
            analyse(WeightedChoice(), [], 8, 1, iterations=10000)
    
    def test_analyse_sample_zero(self):
        ''' test we raise an error if the de novo count is zero
        '''
        rates = WeightedChoice()
        rates.add_choice(200, 1e-5, 'A', 'G')
        rates.add_choice(201, 2e-5, 'C', 'T')
        
        severity = [5, 10]
        with self.assertRaises(ValueError):
            analyse(rates, severity, 0, 0, iterations=10000)
    
    def test_analyse_mismatch(self):
        ''' test for error when the rates and severity lengths are different
        '''
        
        rates = WeightedChoice()
        rates.add_choice(200, 1e-5, 'A', 'G')
        rates.add_choice(201, 2e-5, 'C', 'T')
        
        severity = [5, 10, 5]
        
        with self.assertRaises(ValueError):
            analyse(rates, severity, 8, 1, iterations=100000)
    
    def test_analyse_bigger(self):
        ''' test a more realistically sized data set
        '''
        
        seed(0)
        rates = WeightedChoice()
        pos = sorted(set([ randint(1000, 3000) for x in range(2000) ]))
        
        for x in pos:
            rates.add_choice(x, uniform(1e-10, 1e-7), 'A', 'G')
        
        severity = [ randint(0, 40) for x in pos ]
        
        p = analyse(rates, severity, 150, 4, iterations=10000)
        self.assertAlmostEqual(p, 3e-4, places=2)
        
