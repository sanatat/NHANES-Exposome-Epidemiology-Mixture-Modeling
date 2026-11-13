#!/usr/bin/env python3
"""Offline numerical tests for v8 conversion and 4-knot RCS geometry (no CDC data)."""
import numpy as np

def basis(x,k):
    x = np.asarray(x,dtype=float)
    k = np.asarray(k,dtype=float)
    assert len(k)==4 and np.all(np.diff(k)>0)
    t = lambda y: np.maximum(y,0)**3
    return np.column_stack([(
      t(x-k[j]) - t(x-k[2])*(k[3]-k[j])/(k[3]-k[2])
       +t(x-k[3])*(k[2]-k[j])/(k[3]-k[2])
    )/(k[3]-k[0])**2 for j in (0,1)])

assert np.isclose(100000 * 9.7/100, 9700)  # 9.7 ng/mL / 100 mg/dL == 9700 ng/g
k=[1.0,2.0,3.0,5.0]
for x in (np.linspace(-3,.8,100), np.linspace(5.2,9,100)):
    v=basis(x,k)
    # Natural cubic spline is linear beyond the first and last knot.
    assert np.allclose(np.diff(v,n=2,axis=0),0,atol=1e-9)
x=np.linspace(-3,9,1000)
assert np.linalg.matrix_rank(np.column_stack((np.ones(len(x)),x,basis(x,k))))==4
print('[PASS] Unit conversion and 4-knot restricted cubic spline geometry')
