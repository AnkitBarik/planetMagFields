import numpy as np
import os
import sys
sys.path.append(os.path.abspath('../'))
from planetmagfields import Planet

def test_jupiterBr():
    p = Planet(name='jupiter',r=0.85,nphi=256,info=False,model='jrm33')

    br_ref = np.loadtxt('jupiter/Br_reference085.dat')
    br = p.Br

    percent_err = np.abs( (br_ref - br)/br_ref ) * 100

    np.testing.assert_allclose(percent_err, 0, rtol=0.1,
                                atol=0.1)

def test_potextra_internal():
    p = Planet(name='jupiter',r=10,nphi=256,info=False,model='jrm33')

    p.extrapolate(np.array([10]))
    br_ref = np.squeeze(p.br_ex)
    err = np.abs( (br_ref - p.Br) )

    np.testing.assert_allclose(err, 0, rtol=1e-2, atol=1e-2)

def test_potextra_jupiterMag():

    p = Planet(name='jupiter',nphi=512,info=False,model='jrm33')
    p.extrapolate([2])

    br_ref = np.loadtxt('jupiter/Br_reference.dat')
    bt_ref = np.loadtxt('jupiter/Bt_reference.dat')
    bp_ref = np.loadtxt('jupiter/Bp_reference.dat')

    percent_err = ( np.mean(np.abs(p.br_ex[...,0]*1e3 - br_ref)/br_ref)
                +  np.mean(np.abs(p.btheta_ex[...,0]*1e3 - bt_ref)/bt_ref)
                +  np.mean(np.abs(p.bphi_ex[...,0]*1e3 - bp_ref)/bp_ref)
                ) * 100

    np.testing.assert_allclose(percent_err, 0, rtol=1,
                            atol=1)

def test_potextra_internal_m0():
    p = Planet(name='saturn',r=10,nphi=256,info=False,model='cassini11+')

    p.extrapolate(np.array([10]))
    br_ref = np.squeeze(p.br_ex)
    err = np.abs( (br_ref - p.Br) )

    np.testing.assert_allclose(err, 0, rtol=1e-2, atol=1e-2)

def test_orbit():
    p = Planet(name='earth',r=1,nphi=256,info=False,model='igrf14')
    p.extrapolate([2])
    p.orbit_path([2],[p.theta[10]],[p.phi[10]])

    br_ref,bt_ref,bp_ref = p.br_ex[10,10,0],p.btheta_ex[10,10,0],p.bphi_ex[10,10,0]

    err = np.sqrt((p.br_orb-br_ref)**2
                  +(p.btheta_orb-bt_ref)**2
                  +(p.bphi_orb-bp_ref)**2)

    np.testing.assert_allclose(err, 0, rtol=1e-3, atol=1e-3)


def test_orbit_m0():
    p = Planet(name='saturn',r=1,nphi=256,info=False,model='cassini11+')
    p.extrapolate([2])
    p.orbit_path([2],[p.theta[10]],[p.phi[10]])

    br_ref,bt_ref,bp_ref = p.br_ex[10,10,0],p.btheta_ex[10,10,0],p.bphi_ex[10,10,0]

    err = np.sqrt((p.br_orb-br_ref)**2
                  +(p.btheta_orb-bt_ref)**2
                  +(p.bphi_orb-bp_ref)**2)

    np.testing.assert_allclose(err, 0, rtol=1e-4, atol=1e-4)
