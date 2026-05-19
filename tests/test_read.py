import numpy as np
import os
import sys
sys.path.insert(0, os.path.abspath('../'))
from planetmagfields import Planet, utils


class TestRead:
    def test_read(self):
        glm      = []
        hlm      = []
        lmax_arr = []
        lidx     = [1, 6, 17, 8, 1, 1, 1]
        midx     = [0, 3, 11, 6, 0, 0, 0]
        for k, plname in enumerate(utils.planetlist):
            p = Planet(name=plname, info=False)
            l = lidx[k]
            m = midx[k]
            lmax_arr.append(p.lmax)
            m_clamped = min(m, p.mmax)
            glm.append(p.glm[p.idx[l, m_clamped]])
            hlm.append(p.hlm[p.idx[l, m_clamped]])

        res = np.r_[lmax_arr, lidx, midx, glm, hlm]

        ref_val = [
             3.00000000e+00,  1.30000000e+01,  1.80000000e+01,  1.40000000e+01,
             3.00000000e+00,  3.00000000e+00,  2.00000000e+00,  1.00000000e+00,
             6.00000000e+00,  1.70000000e+01,  8.00000000e+00,  1.00000000e+00,
             1.00000000e+00,  1.00000000e+00,  0.00000000e+00,  3.00000000e+00,
             1.10000000e+01,  6.00000000e+00,  0.00000000e+00,  0.00000000e+00,
             0.00000000e+00, -2.15800003e+02, -1.41399994e+02, -1.48899994e+02,
            -1.55000000e+01,  1.18930000e+04,  9.73200000e+03, -7.11000000e+02,
             0.00000000e+00,  6.15400009e+01,  1.76699997e+02,  0.00000000e+00,
             0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
        ]

        np.testing.assert_allclose(res, ref_val, rtol=1e-9, atol=1e-6)
