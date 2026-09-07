"""Research boundary-row replacement, with explicit omitted weak residuals.

Not a turbulence model or proof of a convergent tau discretization. Global
conservation is retained; eight selected modal equations are NOT enforced.
"""

import copy
import math
import numpy as np
from numpy.polynomial.legendre import legvander
from .edge_enrichment import SquareEvenBasis
from .enriched_transport import EnrichedModalTransport, panel_cumulative
from .coupled_shape_initialization import FixedMeshCoupledTransport


def lift_same_fields(projection, parameters, degree):
    """Embed the same polynomials, without refitting fields or center values."""
    old = projection.basis
    if degree < old.degree:
        raise ValueError('field lift cannot lower the polynomial degree')
    view = copy.copy(projection)
    view.basis = SquareEvenBasis(degree)
    view.count = 2*view.basis.size
    view._grids, view._edges = {}, {}
    result = []
    for block in np.split(np.asarray(parameters, float), 2):
        raw = old.whitener@block
        lifted = np.zeros(view.basis.size)
        for pair, value in zip(old.pairs, raw):
            lifted[view.basis.pairs.index(pair)] = value
        result.append(np.linalg.solve(view.basis.whitener, lifted))
    return view, np.concatenate(result)


def raw_modal_system(projection, matrix, right):
    size = projection.basis.size
    transform = np.linalg.inv(projection.basis.whitener).T
    raw_matrix = np.vstack([transform@matrix[5:5+size], transform@matrix[5+size:]])
    raw_right = np.r_[transform@right[5:5+size], transform@right[5+size:]]
    indices = sorted(range(size), key=lambda j:(sum(projection.basis.pairs[j]), projection.basis.pairs[j]))
    return raw_matrix, raw_right, indices


class BoundaryMomentTransport:
    """Fixed numerical choice: C/H face P0,P2,P4 and momentum P2,P4."""
    def __init__(self, projection, parameters, *, scalar_mixing, thermal_species_ratio,
                 mechanical_work, order=8, angular_order=16, split_angles=True):
        if projection.basis.size < 5:
            raise ValueError('at least five scalar modes required for eight row replacements')
        self.model = EnrichedModalTransport(projection, parameters, scalar_mixing=scalar_mixing,
            thermal_species_ratio=thermal_species_ratio, mechanical_work=mechanical_work)
        self.mesh = FixedMeshCoupledTransport(projection, parameters, scalar_mixing=scalar_mixing,
            thermal_species_ratio=thermal_species_ratio, mechanical_work=mechanical_work,
            order=order, angular_order=angular_order, face_samples=65, split_angles=split_angles)

    def assemble(self):
        model, mesh = self.model, self.mesh
        p, m, size, count = model.projection, model.projection.mixing, model.size, model.count
        base = mesh.evaluate(model.parameters)
        d = model.local(mesh.a, mesh.b)
        ed = model.local(np.ones_like(mesh.face_t), mesh.face_t)
        q, u, rho = d['q'], d['u'], d['rho']
        integrand = np.zeros((len(q),2,count+1))
        integrand[:,0,1:] = -d['bm']
        integrand[:,1,0], integrand[:,1,1:] = d['force'], -d['bp']
        f, fe = np.empty_like(integrand), np.empty((len(mesh.angles),2,count+1))
        for j,(s,k) in enumerate(zip(mesh.slices,mesh.knots)):
            cumulative, endpoint = panel_cumulative(integrand[s],k,mesh.order)
            f[s], fe[j] = cumulative/(2*q[s,None,None]), endpoint/(2*k[-1])
        fm, fp = f[:,0], f[:,1]
        em, ep, ue = fe[:,0], fe[:,1], ed['u']
        fh = (.5*m.wind**2+.5*ue**2)[:,None]*em-ue[:,None]*ep
        ec = ed['y'][:,None]*em
        eh = ed['specific_h'][:,None]*em-fh
        ec[:,0] -= model.area*ed['rho']*mesh.edge_chi*ed['grad_y'][:,0]/(2*p.q0)
        eh[:,0] -= model.area*ed['rho']*mesh.edge_chi*model.ratio*ed['grad_h'][:,0]/(2*p.q0)
        epx = ep-m.wind*math.cos(m.state[3])*em
        # Physical face residual moments; pointwise normalization is diagnostic
        # only. dt=sec(phi)^2 dphi, same normalization on every mirrored face.
        weight = mesh.angular_weights/np.cos(mesh.angles)**2
        leg = legvander(mesh.face_t,4)[:,[0,2,4]]*np.sqrt([1.,5.,9.])
        moments = [np.einsum('n,nk,nr->kr',weight,leg,e,optimize=True) for e in (ec,eh,epx)]
        boundary = np.vstack([moments[0],moments[1],moments[2][1:]])
        raw_matrix, raw_right, ordered = raw_modal_system(p,base['matrix'],base['right'])
        retained = np.r_[ordered[:-4],size+np.array(ordered[:-4])]
        omitted = np.r_[ordered[-4:],size+np.array(ordered[-4:])]
        matrix = np.vstack([base['matrix'][:5],raw_matrix[retained],boundary[:,1:]])
        right = np.r_[base['right'][:5],raw_right[retained],-boundary[:,0]]
        rates = np.zeros(count)
        rates[5:7] = math.cos(m.state[3]),math.sin(m.state[3])
        scales = np.maximum(np.max(abs(matrix[:,model.active]),axis=1),1e-12)
        scaled = matrix[:,model.active]/scales[:,None]
        left, singular, _ = np.linalg.svd(scaled,full_matrices=False)
        self.linear_diagnostics = dict(matrix=matrix,right=right,scales=scales,
            singular_values=singular,left_smallest_vector=left[:,-1],
            incompatible_smallest_component=float(left[:,-1]@((right-matrix@rates)/scales)),
            momentum_mean_affine=moments[2][0],retained_raw_indices=retained,
            boundary_affine=boundary,active=model.active)
        if singular[-1] <= 1e-12*singular[0]:
            raise ValueError('boundary-row system is rank deficient; no pseudoinverse fallback')
        rates[model.active] = np.linalg.solve(scaled,(right-matrix@rates)/scales)
        at = np.r_[1.,rates]
        edge = np.column_stack([ec@at,eh@at,epx@at])/mesh.edge_scales
        rp = fp-u[:,None]*fm
        chi_p = -(rp@at)/(model.area*rho*d['uq'])
        chi_edge = -((ep-ue[:,None]*em)@at)/(model.area*ed['rho']*ed['uq'])
        production = -(2*q*d['uq'])[:,None]*rp
        ew = 16*p.q0*weight
        heat_axial = mesh.weight@d['bh']@rates
        heat_prod, heat_boundary = mesh.weight@production@at, ew@fh@at
        residual = raw_matrix@rates-raw_right
        # P0 is implied by global projected momentum/mass balances for the
        # retained source model. Verify, don't add a duplicate equation.
        momentum_mean = float(moments[2][0]@at)
        mean_scale = max(float(weight@mesh.edge_scales[:,2]),1.)
        linear_error = float(max(abs(matrix@rates-right)/np.maximum(abs(right),1.)))
        global_error = float(max(abs(base['matrix'][:5]@rates-base['source'])/np.maximum(abs(base['source']),1.)))
        retained_error = float(np.max(abs(residual[retained])/np.maximum(abs(raw_right[retained]),1.)))
        heat_error = float(abs(heat_axial-heat_prod+heat_boundary)/max(1.,abs(heat_axial),abs(heat_prod),abs(heat_boundary)))
        mass_error = float(ew@em@at+base['source'][0])
        return dict(rates=rates,matrix=matrix,right=right,source=base['source'],
            original_rates=base['rates'],original_edge_defects=base['edge_defects'],
            original_matrix=base['matrix'],original_right=base['right'],
            raw_modal_residual=residual,omitted_raw_indices=omitted,
            omitted_pairs=[p.basis.pairs[j] for j in ordered[-4:]],
            omitted_scaled_residual=residual[omitted]/np.maximum(abs(raw_right[omitted]),1.),
            retained_modal_scaled_error=retained_error,
            boundary_moments=boundary@at,momentum_mean_residual=momentum_mean,
            momentum_mean_scaled_error=abs(momentum_mean)/mean_scale,
            matrix_condition=float(singular[0]/singular[-1]),
            linear_scaled_error=linear_error,global_scaled_error=global_error,
            weak_heat_scaled_error=heat_error,
            heat_terms=dict(axial=heat_axial,production=heat_prod,boundary=heat_boundary),
            mass_boundary_error=mass_error,
            edge_defects=dict(zip(('hydrogen','heat','momentum'),np.max(abs(edge),axis=0))),
            edge_residual=edge,minimum_chi_momentum=float(min(np.min(chi_p),np.min(chi_edge))),
            maximum_outward_mass=float(np.max(em@at)),
            curvature_half_width=float(abs(rates[3])*math.sqrt(2*p.q0)*m.sn),
            local_numerics_passed=bool(max(linear_error,global_error,retained_error,heat_error,
                abs(momentum_mean)/mean_scale,abs(mass_error)/max(abs(base['source'][0]),1.)) <= 1e-8),
            all_original_weak_equations_enforced=False,adopted=False)
