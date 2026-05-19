import numpy as np
import os
import sys
from copy import deepcopy
sys.path.insert(0, os.path.abspath('../'))
from planetmagfields import Planet


class TestFiltering:
    def test_filt_sat(self):
        p = Planet(name='saturn', nphi=256, info=False, model='cassini11+')
        p.plot_filt(r=0.9, lCutMin=3, iplot=False)

        br_ref = np.loadtxt('./saturn/br_filt_ref_sat.dat')
        percent_err = np.abs((br_ref - p.Br_filt) / br_ref) * 100

        np.testing.assert_allclose(percent_err, 0, rtol=0.1, atol=0.1)

    def test_filt_earth(self):
        p = Planet(name='earth', nphi=256, year=2020, info=False, model='igrf14')
        p.plot_filt(larr=[1, 2, 4], marr=[4], iplot=False)

        br_ref = np.loadtxt('./earth/br_filt_ref_earth.dat')
        percent_err = np.abs((br_ref - p.Br_filt) / br_ref) * 100

        np.testing.assert_allclose(percent_err, 0, rtol=0.1, atol=0.1)

    def test_filt_consistent(self):
        p = Planet(name='earth', year=2020, nphi=128, info=False, model='igrf14')
        p.plot_filt(lCutMax=p.lmax, iplot=False)
        brfilt1 = deepcopy(p.Br_filt)
        p.plot_filt(lCutMin=0, iplot=False)
        brfilt2 = deepcopy(p.Br_filt)

        np.testing.assert_allclose(brfilt1 - brfilt2, 0, rtol=1e-12, atol=1e-12)

