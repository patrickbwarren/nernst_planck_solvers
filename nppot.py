#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# As npeqs.py but run for a series of background salt concentrations,
# reporting the liquid junction potential.

# The --cpoly and --csalt arguments can be specified as either a pair
# of concentrations as cI:cII, or as a single value cI = cII.

import argparse
import numpy as np
import pandas as pd
from numpy import log10
from scipy.integrate import solve_bvp

parser = argparse.ArgumentParser("Nernst-Planck steady-state")
parser.add_argument('-D', '--Darr', default='0.0,2.03,1.33', help='diffusion coeffs, default 0.0,2.03,1.33')
parser.add_argument('-p', '--cpoly', default='0.0:1.0', help='polymer concentrations, default 0.0:1.0')
parser.add_argument('-c', '--csalt', default='0.01,10.0,41', help='background salt range, default 0.01,10.0,41')
parser.add_argument('-v', '--verbosity', action='count', default=0, help='increasing verbosity')
parser.add_argument('-s', '--show', action='store_true', help='plot the density profile')
parser.add_argument('-o', '--output', help='output data for xmgrace, etc')
args = parser.parse_args()

Dp, Ds, Dc = eval(args.Darr)
cpI, cpII = eval(args.cpoly.replace(':', ','))
Δcp = cpII - cpI

# convert --csalt to a logspace

vals = args.csalt.split(',')
start, end, npt = float(vals[0]), float(vals[1]), int(vals[2])

def func(x, y, p):
    φ, cp, cs, Jp, Js = y[0], y[1], y[2], p[0], p[1]
    gradφ = ((1-Dp/Dc)*Jp + (1-Ds/Dc)*Js) / (2*cp + 2*cs)
    return np.vstack((gradφ, -Jp+cp*gradφ, -Js+cs*gradφ))

def bc(ya, yb, p):
    φ0, cp0, cs0 = ya[0], ya[1], ya[2]
    φ1, cp1, cs1 = yb[0], yb[1], yb[2]
    return np.array([φ0, cp0-cpI, cp1-cpII, cs0-cs, cs1-cs])

results = []

for cs in np.logspace(log10(start), log10(end), npt):

    x0 = np.linspace(0, 1, 5) # initial grid

    φ0 = np.zeros_like(x0) # zero electric field, d(φ)/dx = 0
    cp0 = cpI + Δcp*x0 # linear gradient , d(cp)/dx = (cp1-cp0)
    cs0 = cs * np.ones_like(cp0) # uniform concentration, d(cs)/dx = 0

    Jp0, Js0 = Δcp, 0.0 # to correspond to the above solution

    y0, p0 = np.vstack((φ0, cp0, cs0)), np.array([Jp0, Js0])

    res = solve_bvp(func, bc, x0, y0, p=p0, verbose=min(2, args.verbosity))

    Jp, Js = res.p[0], res.p[1]
    Δφ = res.sol(1)[0] - res.sol(0)[0]
    
    # Approximate (Henderson) solution to the NP equations

    σII  = (Dc+Dp)*cpII + (Dc+Ds)*cs
    σI = (Dc+Dp)*cpI + (Dc+Ds)*cs
    Δσ = (Dc+Dp)*Δcp
    Δg = (Dc-Dp)*Δcp
    ΔφH = -Δg/Δσ * np.log(σII/σI)

    results.append((cs, Jp, Js, Δφ, ΔφH))

cols =  ['cs', 'Jp/Dp', 'Js/Ds', 'Δφ', 'Δφ(Hend)']
df = pd.DataFrame(results, columns=cols)

if args.show:
    import matplotlib.pyplot as plt
    df.set_index('cs', inplace=True)
    df[['Δφ', 'Δφ(Hend)']].plot(logx=True)
    plt.show()

elif args.output:

    headers = [f'{col}[{i+1}]' for i, col in enumerate(df.columns)]
    header_row = '#  ' + '  '.join(headers)
    data_rows = df.to_string(index=False).split('\n')[1:]
    with open(args.output, 'w') as f:
        print('\n'.join([header_row] + data_rows), file=f)
    print(f'Data written to {args.output}')
    
else:
    
    print(df.set_index('cs'))
