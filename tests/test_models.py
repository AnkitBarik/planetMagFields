import shutil
import numpy as np
import pytest
from planetmagfields import Planet, get_models
from planetmagfields.libgauss import get_grid, getGauss
from planetmagfields.models import PLANETS, MODELS, planetlist, default_model
from planetmagfields.utils import stdDatDir, get_unit


@pytest.mark.parametrize("planet", planetlist)
def test_default_model_available(planet):
    assert default_model(planet) in get_models(planet)


def test_all_models_have_info():
    for planet in PLANETS:
        for model in get_models(planet):
            assert model in MODELS, f"{planet}_{model}.dat has no MODELS entry"


def test_get_models_path_with_underscore(tmp_path):
    datDir = tmp_path / "my_data"
    datDir.mkdir()
    shutil.copy(stdDatDir + "jupiter_jrm09.dat", datDir)
    assert list(get_models("jupiter", str(datDir))) == ["jrm09"]


def test_unknown_inputs():
    with pytest.raises(ValueError):
        Planet(name="pluto", info=False)
    with pytest.raises(ValueError):
        Planet(units="T", info=False)
    with pytest.raises(FileNotFoundError):
        Planet(model="nope", info=False)


def test_units():
    p_nt = Planet(name="jupiter", nphi=32, units="nT", info=False)
    for units in ["muT", "Gauss"]:
        p = Planet(name="jupiter", nphi=32, units=units, info=False)
        np.testing.assert_allclose(p.Br, p_nt.Br * get_unit(units)[0])


def test_filter_out_of_range():
    p = Planet(name="earth", nphi=32, info=False)
    with pytest.raises(ValueError):
        p.plot_filt(lCutMax=p.lmax + 1, iplot=False)


def test_plot_does_not_change_state():
    import matplotlib
    matplotlib.use("Agg")
    p = Planet(name="saturn", nphi=32, info=False)
    br = p.Br.copy()
    p.plot(r=0.9, proj="hammer")
    assert p.r == 1.0
    np.testing.assert_array_equal(p.Br, br)


@pytest.mark.parametrize("planet", ["earth", "jupiter"])
def test_getGauss_inverts_getB(planet):
    p = Planet(name=planet, nphi=64, units="nT", info=False)
    p2D, th2D, phi, theta = get_grid(64, 32)
    lmax = 6
    glm, hlm = getGauss(lmax, p.Br, 1.0, phi, theta, th2D, p2D, planetname=planet)
    n = p.idx[lmax, lmax] + 1
    scale = np.abs(p.glm).max()
    np.testing.assert_allclose(glm, p.glm[:n], atol=1e-10 * scale)
    np.testing.assert_allclose(hlm, p.hlm[:n], atol=1e-10 * scale)
