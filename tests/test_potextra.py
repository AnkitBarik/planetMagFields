import numpy as np
import os
import sys
sys.path.append(os.path.abspath('../'))
from planetmagfields import Planet

def test_jupiterBr():
    p = Planet(name='jupiter',r=0.85,nphi=256,info=False,model='jrm33')

    br_ref = np.loadtxt('jupiter/Br_reference.dat')
    br = p.Br

    percent_err = np.abs( (br_ref - br)/br_ref ) * 100

    np.testing.assert_allclose(percent_err, 0, rtol=0.1,
                                atol=0.1)

def test_potextra():
    p = Planet(name='jupiter',r=10,nphi=256,info=False,model='jrm33')

    br,btheta,bphi = p.extrapolate(np.array([10]))
    br_ref = np.squeeze(br)
    err = np.abs( (br_ref - p.Br) )

    np.testing.assert_allclose(err, 0, rtol=1e-2, atol=1e-2)