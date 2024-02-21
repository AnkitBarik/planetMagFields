import unittest
import numpy as np
import os
import sys
import shutil
import time

def readData(file):
    return np.loadtxt(file)


class jupiterTest(unittest.TestCase):

    def __init__(self, testName, dir, precision=1e-2):
        super(jupiterTest, self).__init__(testName)

        self.dir = dir
        self.precision = precision
        self.startDir = os.getcwd()
        self.description = "Jupiter's radial field at 0.85 planetary radius"

    def setUp(self):
        # Cleaning when entering
        print('\nDirectory   :           %s' % self.dir)
        print('Description :           %s' % self.description)
        self.startTime = time.time()
        self.cleanDir(self.dir)
        os.chdir(self.dir)

    def list2reason(self, exc_list):
        if exc_list and exc_list[-1][0] is self:
            return exc_list[-1][1]

    def cleanDir(self,dir):
        if os.path.exists("%s/__pycache__" %dir):
            shutil.rmtree("%s/__pycache__" %dir)

    def tearDown(self):
        # Cleaning when leaving
        os.chdir(self.startDir)
        self.cleanDir(self.dir)

        t = time.time()-self.startTime
        st = time.strftime("%M:%S", time.gmtime(t))
        print('Time used   :                            %s' % st)

        if hasattr(self, '_outcome'): # python 3.4+
            if hasattr(self._outcome, 'errors'):  # python 3.4-3.10
                result = self.defaultTestResult()
                self._feedErrorsToResult(result, self._outcome.errors)
            else:  # python 3.11+
                result = self._outcome.result
        else:  # python 2.7-3.3
            result = getattr(self, '_outcomeForDoCleanups',
                             self._resultForDoCleanups)

        error = self.list2reason(result.errors)
        failure = self.list2reason(result.failures)
        ok = not error and not failure

        if ok:
            print('Validating results..                     OK')
        else:
            if error:
                print('Validating results..                     ERROR!')
                print('\n')
                print(result.errors[-1][-1])
            if failure:
                print('Validating results..                     FAIL!')
                print('\n')
                print(result.failures[-1][-1])

    def outputFileDiff(self):
        sys.path.append(os.path.abspath('../..'))
        from planetmagfields import Planet
        p = Planet(name='jupiter',r=0.85,nphi=256,info=False,model='jrm33')

        br_ref = readData('%s/Br_reference.dat' % self.dir)
        br = p.Br

        percent_err = np.abs( (br_ref - br)/br_ref ) * 100

        np.testing.assert_allclose(percent_err, 0, rtol=self.precision,
                                   atol=0.1)
