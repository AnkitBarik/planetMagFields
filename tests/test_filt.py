import numpy as np
import os
import sys
sys.path.append(os.path.abspath('../'))
from planetmagfields import Planet

def test_filt():
    p = Planet(name='saturn',nphi=256,info=False,model='cassini11+')
    p.plot_filt(r=0.9,lCutMin=3,iplot=False)

    br_ref = np.loadtxt('./saturn/br_filt_ref.dat')
    br = p.Br_filt

    percent_err = np.abs( (br_ref - br)/br_ref ) * 100

    np.testing.assert_allclose(percent_err, 0, rtol=0.1,
                                atol=0.1)