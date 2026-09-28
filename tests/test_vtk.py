import numpy as np
import pytest
from planetmagfields import Planet


class TestVtk:
    def test_vtk(self, tmp_path, monkeypatch):
        pytest.importorskip('pyevtk')
        monkeypatch.chdir(tmp_path)
        p = Planet(name='earth', r=1, nphi=256, info=False)
        p.writeVtsFile(potExtra=True, ratio_out=2, nrout=32)
        size = (tmp_path / 'earth.vts').stat().st_size
        size_ref = 67109652
        np.testing.assert_allclose(abs(size - size_ref), 0, rtol=100, atol=100)
