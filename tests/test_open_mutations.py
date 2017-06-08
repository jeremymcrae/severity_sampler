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

from severity.open_mutations import open_mutations, parse_header, parse

class TestLoadGenePy(unittest.TestCase):
    ''' unit test functions to load genes
    '''
    
    def write_temp(self, lines):
        ''' write variant lines to a temporary file and return the file handle
        '''
        temp = tempfile.NamedTemporaryFile()
        
        for x in lines:
            x = '\t'.join(x) + '\n'
            temp.write(x.encode('utf8'))
        
        temp.flush()
        
        return temp
    
    def test_open_mutations(self):
        ''' test that we load mutations correctly
        '''
        
        lines = [['symbol', 'chrom', 'pos', 'ref', 'alt', 'consequence'],
            ['TEST', '1', '200', 'A', 'C', 'missense_variant']]
        
        temp = self.write_temp(lines)
        
        self.assertEqual(open_mutations(temp.name), {'TEST':
            [{'chrom': '1', 'position': 200, 'ref': 'A', 'alt': 'C',
                'consequence': 'missense_variant'}]})
    
    def test_open_mutations_without_header(self):
        ''' test that we fail for files without correct headers
        '''
        
        lines = [['TEST', '1', '200', 'A', 'C', 'missense_variant']]
        
        temp = self.write_temp(lines)
        
        with self.assertRaises(ValueError):
            open_mutations(temp.name)
    
    def test_open_mutations_column_order(self):
        ''' test that we can load files with different column orders
        '''
        
        # scramble the column order
        lines = [['chrom', 'ref', 'alt', 'consequence', 'symbol', 'pos'],
            ['1', 'A', 'C', 'missense_variant', 'TEST', '200']]
        
        temp = self.write_temp(lines)
        
        # check that we still pick up the correct columns based on the header
        self.assertEqual(open_mutations(temp.name), {'TEST':
            [{'chrom': '1', 'position': 200, 'ref': 'A', 'alt': 'C',
                'consequence': 'missense_variant'}]})
        
        # add an extra column, and check that the file is loaded as per usual
        lines = [['chrom', 'ref', 'alt', 'consequence', 'symbol', 'pos', 'extra'],
            ['1', 'A', 'C', 'missense_variant', 'TEST', '200', 'something']]
        
        self.assertEqual(open_mutations(temp.name), {'TEST':
            [{'chrom': '1', 'position': 200, 'ref': 'A', 'alt': 'C',
                'consequence': 'missense_variant'}]})
    
    def test_open_mutations_multiple_variants(self):
        ''' test that we load multiple mutations correctly
        '''
        
        lines = [['symbol', 'chrom', 'pos', 'ref', 'alt', 'consequence'],
            ['TEST', '1', '200', 'A', 'C', 'missense_variant'],
            ['TEST', '1', '300', 'A', 'C', 'missense_variant']]
        
        temp = self.write_temp(lines)
        
        self.assertEqual(open_mutations(temp.name), {'TEST':
            [{'chrom': '1', 'position': 200, 'ref': 'A', 'alt': 'C',
                'consequence': 'missense_variant'},
            {'chrom': '1', 'position': 300, 'ref': 'A', 'alt': 'C',
                'consequence': 'missense_variant'}]})
    
    def test_open_mutations_multiple_genes(self):
        ''' test that we load mutations in multiple genes correctly
        '''
        
        lines = [['symbol', 'chrom', 'pos', 'ref', 'alt', 'consequence'],
            ['TEST', '1', '200', 'A', 'C', 'missense_variant'],
            ['TEST2', '1', '300', 'A', 'C', 'missense_variant']]
        
        temp = self.write_temp(lines)
        
        self.assertEqual(open_mutations(temp.name), {
            'TEST':
                [{'chrom': '1', 'position': 200, 'ref': 'A', 'alt': 'C',
                    'consequence': 'missense_variant'}],
            'TEST2':
                [{'chrom': '1', 'position': 300, 'ref': 'A', 'alt': 'C',
                    'consequence': 'missense_variant'}]
                })
    
    def test_open_mutations_indels(self):
        ''' test that we load indels correctly
        '''
        
        lines = [['symbol', 'chrom', 'pos', 'ref', 'alt', 'consequence'],
            ['TEST', '1', '200', 'A', 'C', 'missense_variant'],
            ['TEST2', '1', '300', 'A', 'CG', 'frameshift_variant']]
        
        temp = self.write_temp(lines)
        
        # check that indels are excluded by default
        self.assertEqual(open_mutations(temp.name), {'TEST':
                [{'chrom': '1', 'position': 200, 'ref': 'A', 'alt': 'C',
                    'consequence': 'missense_variant'}]})
        
        # check that if we pass the indels=True argument, then we pick up indels
        self.assertEqual(open_mutations(temp.name, indels=True), {
            'TEST':
                [{'chrom': '1', 'position': 200, 'ref': 'A', 'alt': 'C',
                    'consequence': 'missense_variant'}],
            'TEST2':
                [{'chrom': '1', 'position': 300, 'ref': 'A', 'alt': 'CG',
                    'consequence': 'frameshift_variant'}]})
    
    def test_open_mutations_consequence(self):
        ''' test that we load indels correctly
        '''
        
        lines = [['symbol', 'chrom', 'pos', 'ref', 'alt', 'consequence'],
            ['TEST', '1', '200', 'A', 'C', 'synonymous_variant']]
        
        temp = self.write_temp(lines)
        
        # check that variants that do not alter protein sequence are excluded
        self.assertEqual(open_mutations(temp.name), {})
    
    def test_open_mutations_consequence_noninteger_positions(self):
        ''' test that we fail if the nucleotide position is not an integer
        '''
        
        lines = [['symbol', 'chrom', 'pos', 'ref', 'alt', 'consequence'],
            ['TEST', '1', '200.0', 'A', 'C', 'missense_variant']]
        
        temp = self.write_temp(lines)
        
        with self.assertRaises(ValueError):
            open_mutations(temp.name)
    
    def test_parse_header(self):
        ''' test that we can parse the header correctly
        '''
        
        header = 'symbol\tchrom\tpos\tref\talt\tconsequence\n'
        self.assertEqual(parse_header(header),
            {'cq': 5, 'pos': 2, 'alt': 4, 'symbol': 0, 'ref': 3, 'chrom': 1})
        
        # scramble the column order
        header = 'chrom\tsymbol\tref\talt\tconsequence\tpos\n'
        self.assertEqual(parse_header(header),
            {'cq': 4, 'pos': 5, 'alt': 3, 'symbol': 1, 'ref': 2, 'chrom': 0})
        
        # add an extra column
        header = 'chrom\tsymbol\tref\talt\tconsequence\tpos\textra\n'
        self.assertEqual(parse_header(header),
            {'cq': 4, 'pos': 5, 'alt': 3, 'symbol': 1, 'ref': 2, 'chrom': 0})
        
        # raise an error if we are missing a column
        with self.assertRaises(ValueError):
            header = 'chrom\tsymbol\tref\talt\tconsequence\n'
            parse_header(header)
    
    def test_parse(self):
        ''' test that we parse lines correctly
        '''
        
        idx = {'cq': 5, 'pos': 2, 'alt': 4, 'symbol': 0, 'ref': 3, 'chrom': 1}
        
        # check a standard line
        line = 'TEST\t1\t200\tA\tC\tmissense_variant\n'
        self.assertEqual(parse(line, idx), ('1', 200, 'A', 'C', 'TEST',
            'missense_variant'))
        
        # check a line with extra columns
        line = 'TEST\t1\t200\tA\tC\tmissense_variant\textra\n'
        self.assertEqual(parse(line, idx), ('1', 200, 'A', 'C', 'TEST',
            'missense_variant'))
        
        # check we raise an error for non-integer positions
        with self.assertRaises(ValueError):
            line = 'TEST\t1\t200.0\tA\tC\tmissense_variant\n'
            parse(line, idx)
        
        # check we raise an error for lines with missing columns
        with self.assertRaises(IndexError):
            line = 'TEST\t1\t200\tA\tC\n'
            parse(line, idx)
        
