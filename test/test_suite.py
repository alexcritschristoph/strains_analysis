#!/usr/bin/env python
'''
Run tests
'''

import os
import glob
import shutil
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from subprocess import call

def load_data_loc():
    return os.path.join(str(os.getcwd()), \
        'test_data/')

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()), \
        'test_backend/testdir/')
    return loc

def get_script_loc(script):
    if script == 'strainRepOneSample.py':
        return os.path.join(str(os.getcwd()), \
            '../strainRepOneSample.py')

class test_strains():
    def setUp(self):
        self.script = get_script_loc('strainRepOneSample.py')

        self.test_dir = load_random_test_dir()

        self.sorted_bam = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa-vs-N5_271_010G1.sorted.bam'
        self.fasta = load_data_loc() + \
            'N5_271_010G1_scaffold_min1000.fa'

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test1()
        self.tearDown()

    def test1(self):
        '''
        Basic test
        '''
        # Set up
        base = self.test_dir + 'test'

        # Run program
        cmd = "{0} {1} {2} -o {3}".format(self.script, self.sorted_bam, \
            self.fasta, base)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produced output
        assert len(glob.glob(base + '*')) == 3

if __name__ == '__main__':
    #test_strainProfiler_breadth().run()
    test_strains().run()
    print('everything is working swimmingly!')
