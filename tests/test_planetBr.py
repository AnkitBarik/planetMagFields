#!/usr/bin/env python3

import numpy as np
import os
import sys

def test_earthBr():
    sys.path.append(os.path.abspath('../'))
    from planetmagfields import Planet

    p = Planet(name='earth',r=1,year=2016,nphi=256,info=False)

    br_ref = np.loadtxt('earth/Br_reference.dat')
    br = p.Br

    br_ref /=1e3

    err = np.abs( (br_ref - br)/br_ref )

    np.testing.assert_allclose(err.mean(), 0, rtol=0.1,
                                atol=0.1)

def test_jupiterBr():
    sys.path.append(os.path.abspath('../'))
    from planetmagfields import Planet
    p = Planet(name='jupiter',r=0.85,nphi=256,info=False,model='jrm33')

    br_ref = np.loadtxt('jupiter/Br_reference.dat')
    br = p.Br

    percent_err = np.abs( (br_ref - br)/br_ref ) * 100

    np.testing.assert_allclose(percent_err, 0, rtol=0.1,
                                atol=0.1)