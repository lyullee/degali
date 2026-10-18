"""Named donor packets, conservative overlaps and interval integration.

The marked-branch closure is exported by the donor and is distinct from its
signed ambient anomaly. A horizontal receiver accepts px/py; vertical impulse
is audited separately. This is a one-way, vertically collapsed diagnostic.
"""
from pathlib import Path
import math
import numpy as np


class Hdf5DenseLayerPacketSource:
    """Return interval-mean rates [mixture, species, px, py] per y bin."""

    def __init__(self, path: str | Path, *, target_y_edges_m: np.ndarray) -> None:
        import h5py
        self.path = Path(path).resolve()
        self.calls = 0
        self.integrated_species_mass_kg = 0.0
        self.integrated_total_mass_kg = 0.0
        self.integrated_horizontal_momentum_kg_m_s = np.zeros(2)
        self.discarded_vertical_impulse_kg_m_s = 0.0
        self.absolute_discarded_vertical_impulse_kg_m_s = 0.0
        self.omitted_species_mass_kg = 0.0
        with h5py.File(self.path, 'r') as handle:
            if handle.attrs.get('dense_layer_packet_schema') != 'degali.marked-layer-donor-packet.v1':
                raise ValueError('handoff lacks named donor x/y momentum; regenerate donor packets')
            self.time_s = np.asarray(handle['time_s'][:],float)
            if (self.time_s.ndim != 1 or self.time_s.size<2
                    or np.any(~np.isfinite(self.time_s)) or np.any(np.diff(self.time_s)<=0)):
                raise ValueError('handoff time history must be finite and increasing')
            self.source_z_edges_m = np.asarray(handle['source_z_edges_m'][:],float)
            self.source_y_edges_m = np.asarray(handle['source_y_edges_m'][:],float)
            self.target_y_edges_m = np.asarray(target_y_edges_m,float)
            for edges in (self.source_z_edges_m,self.source_y_edges_m,self.target_y_edges_m):
                if (edges.ndim!=1 or len(edges)<2 or np.any(~np.isfinite(edges))
                        or np.any(np.diff(edges)<=0)):
                    raise ValueError('interface edges must be finite and strictly increasing')
            overlap = np.maximum(0.,np.minimum(self.target_y_edges_m[1:,None],
                self.source_y_edges_m[None,1:])-np.maximum(self.target_y_edges_m[:-1,None],
                self.source_y_edges_m[None,:-1]))
            dz, source_dy = np.diff(self.source_z_edges_m), np.diff(self.source_y_edges_m)
            fields = []
            for name in ('mixture_mass_flux_kg_m2_s','species_mass_flux_kg_m2_s',
                         'momentum_x_flux_pa','momentum_y_flux_pa','momentum_z_flux_pa'):
                values = np.asarray(handle['dense_layer_packet/'+name][:],float)
                if values.shape != (len(self.time_s),len(dz),len(source_dy)) or np.any(~np.isfinite(values)):
                    raise ValueError('donor packet fields must match the finite time/z/y grid')
                collapsed = np.sum(values*dz[None,:,None],axis=1)
                fields.append(collapsed @ overlap.T)
                if name == 'species_mass_flux_kg_m2_s':
                    omitted = collapsed @ (source_dy-overlap.sum(axis=0))
        self._rates = np.asarray(fields)
        if np.min(self._rates[:2]) < -1e-10 or np.max(self._rates[1]-self._rates[0])>1e-10:
            raise ValueError('one-way layer donor packets must have 0 <= species <= mixture')
        # Sixth component is a scalar audit ledger, stored in the first bin.
        audit = np.zeros_like(self._rates[:1]); audit[0,:,0] = omitted
        self._rates = np.concatenate((self._rates,audit),axis=0)
        dt = np.diff(self.time_s)
        self._integrals = np.concatenate((np.zeros_like(self._rates[:,:1]),
            np.cumsum(.5*(self._rates[:,1:]+self._rates[:,:-1])*dt[None,:,None],axis=1)),axis=1)
        abs_step = self._absolute_linear_integral(self._rates[4,:-1],self._rates[4,1:],dt[:,None])
        self._absolute_vertical_integrals = np.concatenate((np.zeros_like(self._rates[4,:1]),
            np.cumsum(abs_step,axis=0)),axis=0)

    @staticmethod
    def _absolute_linear_integral(a,b,delta):
        ordinary = .5*delta*(np.abs(a)+np.abs(b))
        split = np.divide(delta*(a*a+b*b),2*np.abs(b-a),out=ordinary.copy(),where=np.abs(b-a)>0)
        return np.where(a*b<0,split,ordinary)

    def _integral_at(self, query: float):
        value = float(query)
        tol = 1e-10*max(1.,abs(value))
        if value<self.time_s[0]-tol or value>self.time_s[-1]+tol:
            raise ValueError('requested interval is outside the complete donor history')
        value = np.clip(value,self.time_s[0],self.time_s[-1])
        k = min(max(int(np.searchsorted(self.time_s,value,side='right'))-1,0),len(self.time_s)-2)
        delta = value-self.time_s[k]
        slope = (self._rates[:,k+1]-self._rates[:,k])/(self.time_s[k+1]-self.time_s[k])
        integral = self._integrals[:,k]+delta*self._rates[:,k]+.5*delta*delta*slope
        a = self._rates[4,k]
        absolute = self._absolute_vertical_integrals[k]+self._absolute_linear_integral(a,a+slope[4]*delta,delta)
        return integral,absolute

    def __call__(self,time_s:float,dt_s:float,ny:int,cell_size_m:float)->np.ndarray:
        if not math.isfinite(float(time_s)) or not math.isfinite(float(dt_s)) or dt_s<=0:
            raise ValueError('interface time must be finite and dt positive')
        if not math.isfinite(float(cell_size_m)) or cell_size_m<=0:
            raise ValueError('interface cell size must be positive and finite')
        if ny!=len(self.target_y_edges_m)-1:
            raise ValueError('interface callback ny does not match target edges')
        if not np.allclose(np.diff(self.target_y_edges_m),cell_size_m,rtol=1e-12,atol=1e-12):
            raise ValueError('receiving layer cell size must match target y-bin widths')
        left,left_abs = self._integral_at(time_s)
        right,right_abs = self._integral_at(time_s+dt_s)
        packet = right-left
        self.calls += 1
        self.integrated_total_mass_kg += float(packet[0].sum())
        self.integrated_species_mass_kg += float(packet[1].sum())
        self.integrated_horizontal_momentum_kg_m_s += packet[2:4].sum(axis=1)
        self.discarded_vertical_impulse_kg_m_s += float(packet[4].sum())
        self.absolute_discarded_vertical_impulse_kg_m_s += float((right_abs-left_abs).sum())
        self.omitted_species_mass_kg += float(packet[5].sum())
        return packet[:4]/float(dt_s)

    def close(self):
        """No handle remains after loading the collapsed donor history."""

    def __enter__(self):
        return self

    def __exit__(self,*args):
        self.close()


__all__ = ['Hdf5DenseLayerPacketSource']
