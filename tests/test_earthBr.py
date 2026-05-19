#!/usr/bin/env python3

import numpy as np
import os
import sys
sys.path.insert(0, os.path.abspath('../'))
from planetmagfields import Planet


class TestEarthBr:
    def test_surface_br(self):
        p = Planet(name='earth', r=1, year=2016, nphi=256, info=False)

        br_ref = np.loadtxt('earth/Br_reference.dat')
        br_ref /= 1e3
        err = np.abs((br_ref - p.Br) / br_ref)

        np.testing.assert_allclose(err.mean(), 0, rtol=0.1, atol=0.1)
