import numpy as np
import os
import sys
sys.path.insert(0, os.path.abspath('../'))
from planetmagfields import Planet


class TestVtk:
    def test_vtk(self):
        p = Planet(name='earth', r=1, nphi=256, info=False)
        p.writeVtsFile(potExtra=True, ratio_out=2, nrout=32)
        size = os.path.getsize('earth.vts')
        size_ref = 67109652
        np.testing.assert_allclose(abs(size - size_ref), 0, rtol=100, atol=100)
        os.remove('earth.vts')
