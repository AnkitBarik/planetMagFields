import numpy as np
import os
import sys
sys.path.append(os.path.abspath('../'))
from planetmagfields import Planet

def test_filt_sat():
    p = Planet(name='saturn',nphi=256,info=False,model='cassini11+')
    p.plot_filt(r=0.9,lCutMin=3,iplot=False)

    br_ref = np.loadtxt('./saturn/br_filt_ref_sat.dat')
    br = p.Br_filt

    percent_err = np.abs( (br_ref - br)/br_ref ) * 100

    del p

    np.testing.assert_allclose(percent_err, 0, rtol=0.1,
                                atol=0.1)

def test_filt_earth():
    p = Planet(name='earth',nphi=256,year=2020,info=False,model='igrf14')
    p.plot_filt(larr=[1,2,4],marr=[4],iplot=False)

    br_ref = np.loadtxt('./earth/br_filt_ref_earth.dat')
    br = p.Br_filt

    percent_err = np.abs( (br_ref - br)/br_ref ) * 100

    del p

    np.testing.assert_allclose(percent_err, 0, rtol=0.1,
                                atol=0.1)

def test_filt_consistent():
    from copy import deepcopy

    p = Planet(name="earth",year=2020,nphi=128,info=False,model='igrf14')
    p.plot_filt(lCutMax=p.lmax,iplot=False)
    brfilt1 = deepcopy(p.Br_filt)
    p.plot_filt(lCutMin=0,iplot=False)
    brfilt2 = deepcopy(p.Br_filt)
    err = brfilt1 - brfilt2

    del p

    np.testing.assert_allclose(err, 0, rtol=1e-12, atol=1e-12)
