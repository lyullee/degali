from dataclasses import replace
import numpy as np
import pytest

from degali.addons.binary_mixture_handoff_plane_3d import BinaryMixtureHandoffPlane3DFrame, BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA
from degali.addons.conservative_handoff_hdf5_3d import ConservativeHandoffHdf5Writer3D
from degali.addons.hdf5_dense_layer_interface import Hdf5DenseLayerPacketSource


def history(tmp_path):
    path=tmp_path/'packets.h5'
    ones=np.ones((2,2));zero=ones*0
    base=BinaryMixtureHandoffPlane3DFrame(BINARY_MIXTURE_HANDOFF_PLANE_3D_SCHEMA,
        0.,2,ones,ones*.2,ones*.2,zero,zero,zero,zero,{})
    with ConservativeHandoffHdf5Writer3D(path,source_z_edges_m=np.array((0.,1.,2.)),
            source_y_edges_m=np.array((-1.,0.,1.)),source_id='time-integral') as writer:
        for t,scale in ((0.,1.),(.3,3.),(1.,2.)):
            writer(replace(base,time_s=t,dense_layer_packet={
                'mixture_mass_flux_kg_m2_s':ones*scale,
                'species_mass_flux_kg_m2_s':ones*scale*.2,
                'momentum_x_flux_pa':ones*scale*4.,
                'momentum_y_flux_pa':ones*scale*(-.5),
                'momentum_z_flux_pa':ones*(2*t-1),
            }))
    return path


def test_named_packet_interval_integrals_are_independent_of_receiver_timestep(tmp_path):
    path=history(tmp_path)
    results=[]
    for edges in (np.array((0.,1.)),np.array((0.,.1,.3,.7,1.))):
        with Hdf5DenseLayerPacketSource(path,target_y_edges_m=np.array((-1.,0.,1.))) as source:
            total=np.zeros((4,2))
            for a,b in zip(edges[:-1],edges[1:]):
                total+=source(a,b-a,2,1.)*(b-a)
            results.append(total)
            assert source.absolute_discarded_vertical_impulse_kg_m_s == pytest.approx(2.)
            assert source.discarded_vertical_impulse_kg_m_s == pytest.approx(0.,abs=1e-14)
            assert source.omitted_species_mass_kg == pytest.approx(0.)
    assert results[0] == pytest.approx(results[1],rel=1e-13,abs=1e-13)
    assert results[0][0].sum() == pytest.approx(4*(.3*2+.7*2.5))
    assert results[0][2] == pytest.approx(4.*results[0][0])
    assert results[0][3] == pytest.approx(-.5*results[0][0])


def test_packet_overlap_preserves_inventory_and_zero_fills_exterior(tmp_path):
    path=history(tmp_path)
    target=np.array((-2.8,-1.8,-.8,.2,1.2,2.2))
    with Hdf5DenseLayerPacketSource(path,target_y_edges_m=target) as source:
        packet=source(0.,1.,5,1.)
        assert packet[:,[0,-1]] == pytest.approx(np.zeros((4,2)))
        assert packet[0].sum() == pytest.approx(9.4)
        assert source.omitted_species_mass_kg == pytest.approx(0.)
        with pytest.raises(ValueError,match='outside the complete donor history'):
            source(.9,.2,5,1.)


def test_packet_extension_survives_checkpoint_truncation_and_resume(tmp_path):
    import h5py
    path=history(tmp_path)
    with ConservativeHandoffHdf5Writer3D(path,source_z_edges_m=np.array((0.,1.,2.)),
            source_y_edges_m=np.array((-1.,0.,1.)),source_id='time-integral',resume=True) as writer:
        assert writer.truncate_after(.3)==1
    with h5py.File(path) as f:
        assert f['dense_layer_packet/momentum_x_flux_pa'].shape==(2,2,2)
        assert f['time_s'][-1]==.3


def test_named_donor_packet_enters_layer_with_species_and_xy_momentum_conservation(tmp_path):
    from degali.addons.binary_integral_puff import BinaryPuffThermodynamics
    from degali.addons.binary_mixture_thermodynamics import BinaryIdealGasMixtureClosure
    from degali.addons.dense_layer_patch import solve_dense_layer_patch
    path=history(tmp_path)
    q=np.zeros((4,2,8))
    thermo=BinaryPuffThermodynamics(BinaryIdealGasMixtureClosure('R12',.120913,606.))
    with Hdf5DenseLayerPacketSource(path,target_y_edges_m=np.array((-1.,0.,1.))) as source:
        result=solve_dense_layer_patch(initial_conserved=q,thermo=thermo,
            origin_xy_m=(0.,-1.),cell_size_m=1.,duration_s=1.,wind_xy_m_s=(1.,0.),
            wall_x_height_m=np.zeros((2,9)),wall_y_height_m=np.zeros((3,8)),
            receptor_points_m=np.array([[4.5,0.,.1]]),max_dt_s=.02,
            entrainment_coefficient=0.,drag_coefficient=0.,horizontal_diffusivity_m2_s=0.,
            interface_x_m=2.,interface_source_callback=source)
        assert result.diagnostics['interface_injected_species_kg'] == pytest.approx(source.integrated_species_mass_kg)
        assert result.diagnostics['maximum_species_mass_residual_kg']<1e-11
        assert result.diagnostics['mixture_mass_residual_kg']==pytest.approx(0.,abs=1e-11)
        assert result.diagnostics['momentum_conservation_gate'] is True
        assert source.integrated_horizontal_momentum_kg_m_s[0]==pytest.approx(4.*source.integrated_total_mass_kg)
        assert source.integrated_horizontal_momentum_kg_m_s[1]==pytest.approx(-.5*source.integrated_total_mass_kg)
