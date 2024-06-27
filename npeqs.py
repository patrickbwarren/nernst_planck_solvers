#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The species / function mapping is as follows (default)
#              z    conc    x=0   x=1    diffusion coef
#     polymer  -1    cp      0     1     Dp (eg 0.0 for PSS)
#      co-ion  -1    cs      α     α     Ds (eg 2.03 for Cl)
#  counterion  +1  cp + cs               Dc (eg 1.33 for Na)
#   potential       φ      0    ---
#h ere are five boundary conditions in total.

# The NP equations being solved are:
#  d(φ)/dx = [ (1-Dp/Dc) Jp + (1-Ds/Dc) Js ] / (2 cp + 2 cs)
#  d(cp)/dx = - Jp + cp d(φ)/dx
#  d(cs)/dx = - Js + cs d(φ)/dx
# where Jp and Js are the constant scaled fluxes, which are treated as
# unknown parameters in the problem.  There are three functions (φ,
# cp, cs) and two unknown parameters, which match the above.

# Can run calculations at cs = 0 but need to make sure that cp is not
# also zero as a boundary condition, and vice versa.  Under such
# conditions the Henderson approximate solution becomes exact.

import argparse
import numpy as np
import pandas as pd
from scipy.integrate import solve_bvp

parser = argparse.ArgumentParser("Nernst-Planck steady-state")
parser.add_argument('-D', '--Darr', default='0.0,2.03,1.33', help='diffusion coeffs')
parser.add_argument('-p', '--cpoly', default='0.0:1.0', help='polymer concentrations')
parser.add_argument('-c', '--csalt', default='1.0', help='background salt, default 1.0')
parser.add_argument('-v', '--verbosity', action='count', default=0, help='increasing verbosity')
parser.add_argument('-s', '--show', action='store_true', help='plot the density profile')
parser.add_argument('-o', '--output', help='output data for xmgrace, etc')
args = parser.parse_args()

def double_up(s):
    return s.replace(':', ',') if ':' in s else f'{s},{s}'

Dp, Ds, Dc = eval(args.Darr)
cpI, cpII = eval(double_up(args.cpoly))
csI, csII = eval(double_up(args.csalt))
Δcp = cpII - cpI
Δcs = csII - csI

print('cpI, cpII =', cpI, cpII)
print('csI, csII =', csI, csII)

x0 = np.linspace(0, 1, 5) # initial grid

φ0 = np.zeros_like(x0) # zero electric field, d(φ)/dx = 0
cp0 = cpI + Δcp*x0 # linear gradient , d(cp)/dx = (cp1-cp0)
cs0 = csI + Δcs*x0 # uniform concentration alpha, d(cs)/dx = 0

Jp0, Js0 = Δcp, Δcs # to correspond to the above solution

y0, p0 = np.vstack((φ0, cp0, cs0)), np.array([Jp0, Js0])

def func(x, y, p):
    φ, cp, cs, Jp, Js = y[0], y[1], y[2], p[0], p[1]
    gradφ = ((1-Dp/Dc)*Jp + (1-Ds/Dc)*Js) / (2*cp + 2*cs)
    return np.vstack((gradφ, -Jp+cp*gradφ, -Js+cs*gradφ))

def bc(ya, yb, p):
    φ0, cp0, cs0 = ya[0], ya[1], ya[2]
    φ1, cp1, cs1 = yb[0], yb[1], yb[2]
    return np.array([φ0, cp0-cpI, cp1-cpII, cs0-csI, cs1-csII])

res = solve_bvp(func, bc, x0, y0, p=p0, verbose=min(2, args.verbosity))

print('found parameters Jp, Js =', res.p[0], res.p[1])
print('res.sol(0) =', res.sol(0))
print('res.sol(1) =', res.sol(1))
print('Δφ = ', res.sol(1)[0]-res.sol(0)[0])
    
x = np.linspace(0, 1, 101) ; y = res.sol(x)
φ, cp, cs = y[0], y[1], y[2]

# Approximate (Henderson) solution to the NP equations

σ  = (Dc+Dp)*(cpI+Δcp*x) + (Dc+Ds)*(csI+Δcs*x)
σI = (Dc+Dp)*cpI + (Dc+Ds)*csI
Δσ = (Dc+Dp)*Δcp + (Dc+Ds)*Δcs
Δg = (Dc-Dp)*Δcp + (Dc-Ds)*Δcs 

φH = -Δg/Δσ * np.log(σ/σI)

if args.show:
    import matplotlib.pyplot as plt
    plt.plot(x, φ-np.min(φ), 'g-', label='φ')
    plt.plot(x, φH-np.min(φH), 'g--', label='φ (Henderson)')
    plt.plot(x, cp, 'k-', label='cp')
    plt.plot(x0, cp0, 'k--', label='cp (linear)')
    plt.plot(x, cs, 'b-', label='cs')
    plt.plot(x0, cs0, 'b--', label='cs (linear)')
    # plt.plot(x, cp+cs, 'r-', label='cp + cs')
    plt.legend() ; plt.xlabel("x / L") ; plt.show()
