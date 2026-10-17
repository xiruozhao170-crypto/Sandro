"""Explicit numerical below-ground wind extensions, not recovered observations.

Harmonic extension solves a spherical finite-volume graph in Cartesian vector
components. Nearest-neighbour extension provides a different sensitivity case.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

FILL_VERSION = "spherical-cartesian-v1"


class WindFiller:
    def __init__(self, lat, lon):
        self.lat, self.lon = np.asarray(lat), np.asarray(lon)
        if not (np.allclose(np.diff(lat), np.diff(lat)[0]) and lat[0] == 90 and lat[-1] == -90):
            raise ValueError("Extension requires a regular north-to-south grid including both poles")
        if not np.allclose(np.diff(lon), 360 / len(lon)):
            raise ValueError("Extension requires complete regular global longitudes")
        ny, nx = len(lat), len(lon)
        # Each geographic pole is ONE node, avoiding a coordinate singularity.
        ids = np.empty((ny, nx), dtype=int)
        ids[0] = 0
        ids[1:-1] = np.arange(1, 1 + (ny-2)*nx).reshape(ny-2, nx)
        ids[-1] = 1 + (ny-2)*nx
        self.ids = ids
        self.n = int(ids.max()) + 1
        phi, lam = np.deg2rad(lat)[:, None], np.deg2rad(lon)[None, :]
        self.east = np.stack(np.broadcast_arrays(-np.sin(lam), np.cos(lam), np.zeros_like(phi)), axis=-1)
        self.north = np.stack(np.broadcast_arrays(-np.sin(phi)*np.cos(lam), -np.sin(phi)*np.sin(lam), np.cos(phi)), axis=-1)
        xyz = np.stack(np.broadcast_arrays(np.cos(phi)*np.cos(lam), np.cos(phi)*np.sin(lam), np.sin(phi)), axis=-1)
        self.xyz = self._nodes(xyz)
        dphi, dlam = abs(float(phi[1, 0]-phi[0, 0])), 2*np.pi/nx
        a, b, weights = [], [], []
        for j in range(ny-1):
            a.extend(ids[j]); b.extend(ids[j+1])
            weights.extend(np.full(nx, np.cos((phi[j, 0]+phi[j+1, 0])/2)*dlam/dphi))
        for j in range(1, ny-1):
            a.extend(ids[j]); b.extend(np.roll(ids[j], -1))
            weights.extend(np.full(nx, dphi/(np.cos(phi[j, 0])*dlam)))
        a, b, weights = np.asarray(a), np.asarray(b), np.asarray(weights)
        self.laplacian = coo_matrix((np.r_[weights, weights, -weights, -weights],
                                   (np.r_[a, b, a, b], np.r_[a, b, b, a])),
                                  shape=(self.n, self.n)).tocsr()

    def _nodes(self, values):
        return np.concatenate([values[:1].mean(axis=1), values[1:-1].reshape(-1, 3),
                               values[-1:].mean(axis=1)])

    def fill(self, u, v, method="harmonic"):
        if method not in ("harmonic", "nearest"):
            raise ValueError(f"Unknown wind extension: {method}")
        u, v = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
        if u.shape != self.ids.shape or v.shape != u.shape:
            raise ValueError("Wind shape does not match extension grid")
        if not np.array_equal(np.isnan(u), np.isnan(v)):
            raise ValueError("u/v missing masks differ; inspect source before extending")
        if np.isinf(u).any() or np.isinf(v).any():
            raise ValueError("Infinite wind is not a below-ground missing value")
        missing = np.isnan(u)
        if not missing.any():
            return u.copy(), v.copy()
        if missing.all():
            raise ValueError("Cannot extend an entirely missing global field")
        cartesian = u[..., None]*self.east + v[..., None]*self.north
        vectors = self._nodes(np.nan_to_num(cartesian))
        known = np.ones(self.n, dtype=bool)
        known[1:-1] = ~missing[1:-1].ravel()
        for row, node in ((0, 0), (-1, self.n-1)):
            valid = ~missing[row]
            known[node] = valid.any()
            if valid.any():
                vectors[node] = cartesian[row, valid].mean(axis=0)
        absent = ~known
        if not absent.any():
            # Only some duplicate longitude representations of a pole were
            # missing; the geographic node is already fixed by observations.
            pass
        elif method == "harmonic":
            matrix = self.laplacian[absent][:, absent]
            rhs = -(self.laplacian[absent][:, known] @ vectors[known])
            vectors[absent] = spsolve(matrix, rhs)
            residual = matrix @ vectors[absent] - rhs
            if np.max(abs(residual), initial=0) > 1e-8 * max(1, np.max(abs(rhs), initial=0)):
                raise ArithmeticError("Harmonic extension failed its linear-system residual check")
        elif method == "nearest":
            tree = cKDTree(self.xyz[known])
            _, nearest = tree.query(self.xyz[absent])
            vectors[absent] = vectors[known][nearest]
        completed = vectors[self.ids]
        uf, vf = u.copy(), v.copy()
        uf[missing] = np.sum(completed*self.east, axis=-1)[missing]
        vf[missing] = np.sum(completed*self.north, axis=-1)[missing]
        if not np.isfinite(uf).all() or not np.isfinite(vf).all():
            raise ArithmeticError("Nonfinite extended wind")
        return uf, vf
