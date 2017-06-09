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


from intervaltree import Interval, IntervalTree

# ranges and weights for missense sites in regions without regional constraint
unconstrained = [
    ( 0,  5, 0.89378999169559),
    ( 5, 10, 2.45095807132636),
    (10, 15, 1.23783036756664),
    (15, 20, 1.04908176145445),
    (20, 25, 2.13085809157198),
    (25, 30, 4.97619976726511),
    (30, 35, 6.73650544131241),
    (35, 40, 5.68752647734537),
    (40, 1000, 5.68752647734537)]

# ranges and weights for missense sites in regions with regional constraint
constrained = [
    ( 0,  5, 0.0),
    ( 5, 10, 4.05516581596172),
    (10, 15, 2.75708159237827),
    (15, 20, 4.81275329358394),
    (20, 25, 7.41041424690547),
    (25, 30, 16.5358474569603),
    (30, 35, 19.0139355018205),
    (35, 40, 35.7654385873813),
    (40, 1000, 35.7654385873813)]

weights = {
    'altering': {
        'unconstrained': IntervalTree( Interval(*x) for x in unconstrained ),
        'constrained': IntervalTree( Interval(*x) for x in constrained )
        },
    'truncating': 30.4986359738963
    }
