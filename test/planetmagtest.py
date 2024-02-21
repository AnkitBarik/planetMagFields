#!/usr/bin/env python3
import argparse
import os
import sys
import unittest
import jupiter.unitTest
import earth.unitTest

__version__="1.0"

def getParser():
    """
    Get script option
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s '+__version__,
                        help="Show program's version number and exit.")

    return parser


def print_logo():
    print("           __                 __  __  ___            _______      __    __    ")
    print("    ____  / /___ _____  ___  / /_/  |/  /___ _____ _/ ____(_)__  / /___/ /____")
    print("   / __ \/ / __ `/ __ \/ _ \/ __/ /|_/ / __ `/ __ `/ /_  / / _ \/ / __  / ___/")
    print("  / /_/ / / /_/ / / / /  __/ /_/ /  / / /_/ / /_/ / __/ / /  __/ / /_/ (__  ) ")
    print(" / .___/_/\__,_/_/ /_/\___/\__/_/  /_/\__,_/\__, /_/   /_/\___/_/\__,_/____/  ")
    print("/_/                                        /____/                             ")



def getSuite(startdir, precision):
    """
    Construct test suite
    """
    suite = unittest.TestSuite()

    suite.addTest(jupiter.unitTest.jupiterTest('outputFileDiff',
                                                '%s/jupiter' %startdir,
                                                precision=precision))
    suite.addTest(earth.unitTest.earthTest('outputFileDiff',
                                                '%s/earth' %startdir,
                                                precision=precision))
    return suite


if __name__ == '__main__':
    precision = 0.1 # relative tolerance between expected and actual result
    startdir = os.getcwd()

    parser = getParser()

    # Initialisation
    print_logo()

    # Run the auto-test suite
    print('                     Running test suite                               ')
    print('----------------------------------------------------------------------')
    suite = getSuite(startdir,precision)
    runner = unittest.TextTestRunner(verbosity=0)
    ret = not runner.run(suite).wasSuccessful()

    sys.exit(ret)
