import numpy as np
import pytest
from degali.addons.ambient_vector_exchange import entrained_ambient_rates


def test_exact_mass_momentum_and_energy():
    assert entrained_ambient_rates(2,[3,4,0],specific_enthalpy=-5)==pytest.approx([2,0,6,8,0,15])


def test_rotation_covariance():
    q,_=np.linalg.qr(np.array([[1,2,3],[4,2,1],[2,1,4.]]))
    v=np.array([2.,-3.,4.]);a=entrained_ambient_rates(5,v,specific_enthalpy=8)
    b=entrained_ambient_rates(5,q@v,specific_enthalpy=8)
    assert b[:2]==pytest.approx(a[:2])
    assert b[2:5]==pytest.approx(q@a[2:5])
    assert b[5]==pytest.approx(a[5])


def test_enthalpy_reference_and_zero_mass():
    a=entrained_ambient_rates(.3,[1,2,3],specific_enthalpy=-40)
    b=entrained_ambient_rates(.3,[1,2,3],specific_enthalpy=1e6-40)
    assert b[5]-a[5]==pytest.approx(.3e6)
    assert entrained_ambient_rates(0,[1,2,3],specific_enthalpy=-40)==pytest.approx(np.zeros(6))


def test_broadcast_and_opposed_velocity_not_clipped():
    out=entrained_ambient_rates([1,2],[[-1,0,0],[0,3,0]],specific_enthalpy=0)
    assert out.shape==(2,6)
    assert out[:,2:5]==pytest.approx(np.array([[-1,0,0],[0,6,0]]))
    assert out[:,5]==pytest.approx([.5,9])


@pytest.mark.parametrize('m,v,h',[(-1,[1,2,3],0),(1,[1,2],0),(1,[1,2,np.nan],0),(1,[1,2,3],np.inf)])
def test_bad_inputs(m,v,h):
    with pytest.raises(ValueError): entrained_ambient_rates(m,v,specific_enthalpy=h)
