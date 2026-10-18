import pytest

from degali.addons.hybrid_interface import (
    HybridInterfaceCoupler,
    InterfaceFluxLedger,
    conservative_interface_transfer,
)


def test_interface_transfer_is_bounded_and_mass_conservative():
    transferred, remainder = conservative_interface_transfer(
        donor_mass_kg=2.5, requested_mass_kg=4.0,
    )
    assert transferred == pytest.approx(2.5)
    assert remainder == pytest.approx(0.0)
    ledger = InterfaceFluxLedger().add(incoming_kg=transferred)
    ledger = ledger.add(outgoing_kg=transferred)
    assert ledger.residual_kg == pytest.approx(0.0)


def test_interface_ledger_rejects_negative_flux():
    with pytest.raises(ValueError):
        InterfaceFluxLedger().add(incoming_kg=-1.0)


def test_coupler_records_both_directions_and_balances():
    coupler = HybridInterfaceCoupler()
    assert coupler.transfer_to_patch(3.0, 1.25) == pytest.approx(1.25)
    assert coupler.transfer_to_coarse(2.0, 1.25) == pytest.approx(1.25)
    coupler.require_balanced()


def test_exchange_step_preserves_total_mass():
    coupler = HybridInterfaceCoupler()
    coarse, patch = coupler.exchange_step(
        coarse_mass_kg=10.0, patch_mass_kg=2.0,
        coarse_to_patch_kg=1.5, patch_to_coarse_kg=0.25,
    )
    assert coarse + patch == pytest.approx(12.0)
