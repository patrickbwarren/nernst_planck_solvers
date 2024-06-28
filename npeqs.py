#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# The species mapping is as follows
#              z    conc     x=0     x=1       diffusion coef
#     polymer  -1    cp      cpI     cpII     Dp (eg 0.0 for PSS)
#      co-ion  -1    cs      csI     csII     Ds (eg 2.03 for Cl)
#  counterion  +1  cp + cs                    Dc (eg 1.33 for Na)
#   potential        φ        0      ---
# There are five boundary conditions in total.

# The NP equations being solved are:
#  d(φ)/dx = [ (1-Dp/Dc) Jp' + (1-Ds/Dc) Js' ] / (2 cp + 2 cs)
#  d(cp)/dx = - Jp' + cp d(φ)/dx
#  d(cs)/dx = - Js' + cs d(φ)/dx
# where Jp' = Jp/Dp and Js'=Js/Ds are the constant scaled fluxes,
# which are treated as unknown parameters in the problem.  There are
# three functions (φ, cp, cs) and two unknown parameters.

# Can run calculations at cs = 0 but need to make sure that cp is not
# also zero as a boundary condition, and vice versa.  Under such
# conditions the Henderson approximate solution becomes exact.

# The --cpoly and --csalt arguments can be specified as either a pair
# of concentrations as cI:cII, or as a single value cI = cII.

import argparse
import numpy as np
from scipy.integrate import solve_bvp

parser = argparse.ArgumentParser("Nernst-Planck steady-state")
parser.add_argument('-D', '--Darr', default='0.0,2.03,1.33', help='diffusion coeffs, default 0.0,2.03,1.33')
parser.add_argument('-p', '--cpoly', default='1.0:0.0', help='polymer concentrations, default 1.0:0.0')
parser.add_argument('-c', '--csalt', default='1.0', help='background salt, default 1.0')
parser.add_argument('-v', '--verbosity', action='count', default=0, help='increasing verbosity')
parser.add_argument('-H', '--henderson', action='store_true', help='show the Henderson profile')
parser.add_argument('-t', '--total', action='store_true', help='plot the co-ion profile')
parser.add_argument('-l', '--legend', action='store_true', help='plot the legend')
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

x0 = np.linspace(0, 1, 5) # initial grid

φ0 = np.zeros_like(x0) # zero electric field, d(φ)/dx = 0
cp0 = cpI + Δcp*x0 # linear gradient , d(cp)/dx = (cp1-cp0)
cs0 = csI + Δcs*x0 # also a linear gradient

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
    
x = np.linspace(0, 1, 101) ; y = res.sol(x)
φ, cp, cs = y[0], y[1], y[2]
Δφ = res.sol(1)[0] - res.sol(0)[0]

# Approximate (Henderson) solution to the NP equations

σ  = (Dc+Dp)*(cpI+Δcp*x) + (Dc+Ds)*(csI+Δcs*x)
σII = (Dc+Dp)*cpII + (Dc+Ds)*csII
σI = (Dc+Dp)*cpI + (Dc+Ds)*csI
Δσ = (Dc+Dp)*Δcp + (Dc+Ds)*Δcs
Δg = (Dc-Dp)*Δcp + (Dc-Ds)*Δcs 
φH = -Δg/Δσ * np.log(σ/σI)
ΔφH = -Δg/Δσ * np.log(σII/σI)

if args.show:

    import matplotlib.pyplot as plt
    plt.plot(x, φ-np.min(φ), 'g--', label='φ')
    plt.plot(x, cp, 'k-', label='cp')
    plt.plot(x, cs, 'b-', label='cs')
    if args.henderson:
        plt.plot(x, φH-np.min(φH), 'g-.', label='φ (Henderson)')
        plt.plot(x0, cp0, 'k--', label='cp (linear)')
        plt.plot(x0, cs0, 'b--', label='cs (linear)')
    if args.total:
        plt.plot(x, cp+cs, 'r-', label='cp + cs')
    if args.legend:
        plt.legend()
    plt.xlabel("x / L") ; plt.show()

elif args.output:

    import pandas as pd
    df = pd.DataFrame({'x':x, 'phi':φ-φ[0], 'cp':cp, 'cs':cs, 'cp+cs':(cp+cs), 'phi(Hend)':φH})
    headers = [f'{col}[{i+1}]' for i, col in enumerate(df.columns)]
    header_row = '#  ' + '  '.join(headers)
    data_rows = df.to_string(index=False).split('\n')[1:]
    with open(args.output, 'w') as f:
        print('\n'.join([header_row] + data_rows), file=f)
    print(f'Data written to {args.output}')

else:

    print('cpI, cpII =', cpI, cpII)
    print('csI, csII =', csI, csII)
    print('Dp, Ds, Dc =', Dp, Ds, Dc)

    print('Jp/Dp, Js/Ds =', res.p[0], res.p[1])
    print('Δφ       =\t', Δφ)
    print('Δφ(Hend) =\t', ΔφH)
