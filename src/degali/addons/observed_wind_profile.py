"""Explicit observed speed/reference-height input without mutating CONTROL."""
from copy import copy
import math

from .energy_crosswind import IndependentEnergyCrosswind
from ..core.atmosphere import friction_velocity


def with_observed_wind(base,*,speed,reference_height):
    if type(base) is not IndependentEnergyCrosswind:
        raise ValueError('explicit density CONTROL required; cannot discard another option')
    if not math.isfinite(speed) or speed < 0:
        raise ValueError('finite nonnegative observed wind required')
    if not math.isfinite(reference_height) or reference_height <= 0:
        raise ValueError('finite positive measurement height required')
    result=copy(base)
    result.jetplume=copy(base.jetplume)
    jp=result.jetplume
    jp.u0=float(speed)
    jp.z0=float(reference_height)
    jp.ustar=friction_velocity(jp.u0,jp.z0,jp.zr,jp.rml)
    return result
