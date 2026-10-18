"""Bounded species-volume reconstruction of a reduced vertical column.

This is an observation operator, not a vertical momentum solver. At fixed
pressure/temperature the integral of mole fraction is S/rho_species, not the
integral of a reshaped mass fraction. Kernels are integrated over vertical
bins before overlapping branches are combined. A capacity-constrained
redistribution preserves that volume instead of clipping it away at 100%.
No observation values, sensor indices or fitted dose multipliers enter here.
"""
from dataclasses import dataclass

import numpy as np
from scipy.special import erfc, ndtr


@dataclass
class ColumnObservation:
    edges_m: np.ndarray
    mole_fraction: np.ndarray
    target_species_volume_m: np.ndarray
    maximum_relative_volume_residual: float
    saturated_column_count: int
    slopes_per_m: np.ndarray
    component_labels: tuple[str, ...] = ()
    component_bin_volume_m: np.ndarray | None = None

    def sample(self, height_m, column_indices=None):
        columns = (np.arange(len(self.edges_m)) if column_indices is None
                   else np.asarray(column_indices, int))
        z = np.broadcast_to(np.asarray(height_m, float), columns.shape)
        edges = self.edges_m[columns]
        indices = np.sum(edges <= z[:, None], axis=1) - 1
        valid = (z >= 0.) & (z < edges[:, -1])
        indices = np.clip(indices, 0, self.mole_fraction.shape[1] - 1)
        centre = .5*(self.edges_m[columns, indices]+self.edges_m[columns, indices+1])
        value = self.mole_fraction[columns, indices]+self.slopes_per_m[columns, indices]*(z-centre)
        return np.where(valid, value, 0.)

    def sample_components(self, height_m, column_indices=None):
        """Attribute the bounded total to pre-mixing bin-volume shares.

        Same-species branches are not physically distinguishable after they
        overlap. These values are therefore a diagnostic attribution, not
        separately conserved prognostic species. Their sum is exactly the
        returned total at every query point.
        """
        if self.component_bin_volume_m is None:
            return np.empty((0, np.size(height_m)))
        columns = (np.arange(len(self.edges_m)) if column_indices is None
                   else np.asarray(column_indices, int))
        z = np.broadcast_to(np.asarray(height_m, float), columns.shape)
        edges = self.edges_m[columns]
        indices = np.sum(edges <= z[:, None], axis=1) - 1
        valid = (z >= 0.) & (z < edges[:, -1])
        indices = np.clip(indices, 0, self.mole_fraction.shape[1] - 1)
        volume = self.component_bin_volume_m[:, columns, indices]
        total_volume = volume.sum(axis=0)
        share = np.divide(volume, total_volume[None, :],
                          out=np.zeros_like(volume), where=total_volume[None, :] > 0.)
        return share * self.sample(z, columns)[None, :] * valid[None, :]


def _limited_slopes(values, widths, edges):
    """Monotone cell-linear reconstruction with unchanged cell integrals."""
    n = widths.shape[1]
    index = np.arange(n)[None, :]
    valid = widths > 0.
    previous = np.maximum.accumulate(np.where(valid, index, -1), axis=1)
    previous = np.concatenate((np.full((len(values), 1), -1), previous[:, :-1]), axis=1)
    following = np.minimum.accumulate(np.where(valid, index, n)[:, ::-1], axis=1)[:, ::-1]
    following = np.concatenate((following[:, 1:], np.full((len(values), 1), n)), axis=1)
    has_neighbours = valid & (previous >= 0) & (following < n)
    previous, following = np.clip(previous, 0, n-1), np.clip(following, 0, n-1)
    centres = .5*(edges[:, :-1]+edges[:, 1:])
    left = np.take_along_axis(values, previous, axis=1)
    right = np.take_along_axis(values, following, axis=1)
    dl = centres-np.take_along_axis(centres, previous, axis=1)
    dr = np.take_along_axis(centres, following, axis=1)-centres
    gl = np.divide(values-left, dl, out=np.zeros_like(values), where=has_neighbours & (dl>0.))
    gr = np.divide(right-values, dr, out=np.zeros_like(values), where=has_neighbours & (dr>0.))
    bound = np.divide(2.*np.minimum(values, 1.-values), widths,
                      out=np.zeros_like(values), where=valid)
    return np.where(gl*gr>0., np.sign(gl)*np.minimum(np.minimum(abs(gl), abs(gr)), bound), 0.)


def _bounded_volume_fractions(volumes, widths):
    """Weighted water filling: sum(c*dz)=sum(volumes), 0<=c<=1."""
    raw = np.divide(volumes, widths, out=np.zeros_like(volumes), where=widths > 0.)
    target = volumes.sum(axis=1)
    affected = np.any(raw > 1., axis=1)
    result = raw.copy()
    if np.any(affected):
        v, w, c = volumes[affected], widths[affected], raw[affected]
        support = np.where(c > 0., w, 0.).sum(axis=1)
        if np.any(target[affected] > support + 1.e-11 * np.maximum(support, 1.)):
            raise ValueError('species volume exceeds vertical profile support capacity')
        order = np.argsort(-c, axis=1)
        ws = np.take_along_axis(w, order, axis=1)
        vs = np.take_along_axis(v, order, axis=1)
        cs = np.take_along_axis(c, order, axis=1)
        capped_width = np.concatenate((np.zeros((len(v), 1)),
                                       np.cumsum(ws, axis=1)[:, :-1]), axis=1)
        remaining_volume = np.cumsum(vs[:, ::-1], axis=1)[:, ::-1]
        multiplier = np.divide(target[affected, None] - capped_width,
                               remaining_volume,
                               out=np.full_like(v, np.inf),
                               where=remaining_volume > 0.)
        product = np.multiply(multiplier, cs, out=np.full_like(cs, np.inf),
                              where=np.isfinite(multiplier))
        admissible = (multiplier >= 0.) & (product <= 1. + 1.e-12)
        filled = np.isclose(target[affected], support, rtol=1.e-12, atol=1.e-14)
        if np.any(~np.any(admissible, axis=1) & ~filled):
            raise RuntimeError('bounded vertical reconstruction has no admissible partition')
        selected = np.argmax(admissible, axis=1)
        scale = multiplier[np.arange(len(v)), selected]
        scale = np.where(filled, 1., scale)
        repaired = np.minimum(c * scale[:, None], 1.)
        result[affected] = np.where(filled[:, None], (c > 0.).astype(float), repaired)
    return result, int(affected.sum())


def reconstruct_column(*, height_m, core_species_volume_m, exponent=0.,
                       tail_kind='none', tail_scale_m=1., tail_extent_m=0.,
                       gaussian_branches=(), core_components=(), bins=64):
    """Return piecewise-constant vol/vol profiles for independent columns.

    Gaussian tuples are (species volume per area, centre, sigma). Their
    support is [0, centre+8 sigma], renormalized explicitly. The mesh depends
    on column/branch geometry only, never on receptor elevation. Each kernel
    has its own nodes; sorting their union resolves narrow branches as well
    as a deep initial cloud. ``bins`` controls numerical observation accuracy.
    """
    height = np.atleast_1d(np.asarray(height_m, float))
    def vector(value):
        return np.broadcast_to(np.asarray(value, float), height.shape)
    core, power = vector(core_species_volume_m), vector(exponent)
    scale, extent = vector(tail_scale_m), vector(tail_extent_m)
    if (height.ndim != 1 or any(np.any(~np.isfinite(a)) for a in
                              (height, core, power, scale, extent))
            or np.any(height < 0.) or np.any(core < 0.) or np.any(power < 0.)
            or np.any(scale <= 0.) or np.any(extent < 0.)
            or np.any((height == 0.) & (core > 0.))):
        raise ValueError('invalid finite column geometry or inventory')
    if isinstance(bins, bool) or not isinstance(bins, int) or bins < 8:
        raise ValueError('bins must be an integer >=8')
    if tail_kind not in {'none', 'exponential', 'diffusive'}:
        raise ValueError('unknown conservative tail kind')
    nodes = np.linspace(0., 1., bins + 1)[None, :]
    edge_sets = [height[:, None] * nodes]
    if tail_kind != 'none':
        edge_sets.append(height[:, None] + extent[:, None] * nodes**2)
    branches = []
    for branch_index, branch in enumerate(gaussian_branches):
        if len(branch)==3:
            volume, centre, sigma = branch
            label=f'gaussian_{branch_index}'
        elif len(branch)==4:
            label, volume, centre, sigma = branch
        else:
            raise ValueError('Gaussian branches require (volume, centre, sigma) or a leading label')
        if not isinstance(label,str) or not label:
            raise ValueError('component labels must be non-empty strings')
        volume, centre, sigma = vector(volume), vector(centre), vector(sigma)
        if (any(np.any(~np.isfinite(a)) for a in (volume, centre, sigma))
                or np.any(volume < 0.) or np.any(centre < 0.) or np.any(sigma <= 0.)):
            raise ValueError('invalid Gaussian branch inventory or geometry')
        branches.append((label, volume, centre, sigma))
        edge_sets.append(np.maximum(centre[:, None] + sigma[:, None]
                                    * (16. * nodes - 8.), 0.))
    edges = np.sort(np.concatenate(edge_sets, axis=1), axis=1)
    # Nearly coincident edges from different kernels must not create a tiny
    # cell whose CDF subtraction is dominated by roundoff at a sensor height.
    tolerance = 64.*np.finfo(float).eps*np.maximum(abs(edges[:, 1:]), 1.e-30)
    distinct = np.concatenate((np.ones((len(height), 1), bool),
                               np.diff(edges, axis=1)>tolerance), axis=1)
    retained = np.maximum.accumulate(np.where(distinct, np.arange(edges.shape[1])[None, :], 0), axis=1)
    edges = np.take_along_axis(edges, retained, axis=1)
    lo, hi = edges[:, :-1], edges[:, 1:]
    widths = hi - lo
    h = np.maximum(height[:, None], 1.e-300)
    # Survival differences avoid cancellation near the top of a steep core.
    weight = (np.maximum(1. - lo/h, 0.)**(power[:, None]+1.)
              - np.maximum(1. - hi/h, 0.)**(power[:, None]+1.))
    if tail_kind != 'none':
        def tail_integral(distance):
            s = scale[:, None]
            if tail_kind == 'exponential':
                return -s * np.expm1(-distance / s)
            u = distance/s
            return s * (u*erfc(u) - np.expm1(-u*u)/np.sqrt(np.pi))
        lower = np.clip(lo-height[:, None], 0., extent[:, None])
        upper = np.clip(hi-height[:, None], 0., extent[:, None])
        integral = tail_integral(extent[:, None])
        total = h + integral
        weight = (weight*h + tail_integral(upper)-tail_integral(lower))/total
    component_labels=['core']
    component_volumes=[core[:, None]*weight]
    target = core.copy()
    for label, volume in core_components:
        if not isinstance(label,str) or not label:
            raise ValueError('component labels must be non-empty strings')
        volume=vector(volume)
        if np.any(~np.isfinite(volume)) or np.any(volume<0.):
            raise ValueError('core component volume must be finite and nonnegative')
        component_labels.append(label);component_volumes.append(volume[:,None]*weight)
        target+=volume
    for label, volume, centre, sigma in branches:
        zlo, zhi = (lo-centre[:, None])/sigma[:, None], (hi-centre[:, None])/sigma[:, None]
        # Evaluate upper tails as survival probabilities, avoiding 1-1.
        probability = np.where(zlo > 0., ndtr(-zlo)-ndtr(-zhi), ndtr(zhi)-ndtr(zlo))
        top = centre+8.*sigma
        probability = np.where(lo >= top[:, None], 0., probability)
        crossing = (lo < top[:, None]) & (hi > top[:, None])
        probability = np.where(crossing, ndtr(-zlo)-ndtr(-8.), probability)
        normalizer = ndtr(8.)-ndtr(-centre/sigma)
        component_labels.append(label)
        component_volumes.append(volume[:, None]*np.maximum(probability, 0.)/normalizer[:, None])
        target += volume
    if len(set(component_labels))!=len(component_labels):
        raise ValueError('component labels must be unique')
    component_volumes=np.stack(component_volumes,axis=0)
    volumes=component_volumes.sum(axis=0)
    fractions, count = _bounded_volume_fractions(volumes, widths)
    error = np.abs(np.sum(fractions*widths, axis=1)-target)
    relative = error/np.maximum(target, 1.e-20)
    maximum = float(np.max(relative, initial=0.))
    if maximum > 1.e-9 or np.any(~np.isfinite(fractions)):
        raise RuntimeError('vertical species-volume conservation gate failed')
    slopes = _limited_slopes(fractions, widths, edges)
    return ColumnObservation(edges, fractions, target, maximum, count, slopes,
                             tuple(component_labels),component_volumes)
