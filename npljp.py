#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Code to solve one-dimensional steady-state Nernst-Planck equations
# for a restricted-diffusion junction connecting polymer + salt
# reservoirs.  See J. Newman and K. E. Thomas-Alyea, Electrochemical
# Systems, 3rd ed. (Wiley, Hoboken, NJ, 2004) for more details.

# This code is copyright (c) 2024 Patrick B Warren (STFC).
# Email: patrick.warren{at}stfc.ac.uk.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# As npeqs.py but run for a series of background salt concentrations,
# reporting the liquid junction potential.

import argparse
import numpy as np
import pandas as pd
from numpy import log10
from scipy.optimize import root
from numpy import log as ln

parser = argparse.ArgumentParser("Nernst-Planck steady-state")
parser.add_argument('-D', '--Darr', default='0.0,2.03,1.33', help='diffusion coeffs, default 0.0,2.03,1.33')
parser.add_argument('-p', '--cpoly', default='1.0:0.1', help='polymer concentrations, default 1.0:0.01')
parser.add_argument('-c', '--csalt', default='1e-3,10.0,41', help='background salt range, default 1e-3,10.0,41')
parser.add_argument('-v', '--verbosity', action='count', default=0, help='increasing verbosity')
parser.add_argument('-H', '--henderson', action='store_true', help='show the Henderson profile')
parser.add_argument('-l', '--legend', action='store_true', help='plot the legend')
parser.add_argument('-s', '--show', action='store_true', help='plot the density profile')
parser.add_argument('-o', '--output', help='output data for xmgrace, etc')
args = parser.parse_args()

Dp, Ds, Dc = eval(args.Darr)
cpI, cpII = eval(args.cpoly.replace(':', ','))
Δcp = cpII - cpI

vals = args.csalt.split(',') 
start, end, npt = float(vals[0]), float(vals[1]), int(vals[2])

def func(x, a, b):
    A, Jp, Js = x # extract variables
    K = ((1-Dp/Dc)*Jp + (1-Ds/Dc)*Js)/2 # a derived quantity
    cp0 = Jp/(Jp+Js)*a + A * a**(K/b) # cp at x = 0
    cp1 = Jp/(Jp+Js)*(a+b) + A * (a+b)**(K/b) # cp at x = L
    return([cp0-cpI, cp1-cpII, b+Jp+Js-K]) # bcs ; constraint on K

results = []

for cs in np.logspace(log10(start), log10(end), npt):

    a, b = cpI + cs, Δcp  # as in cp + cs = a + b x

    x0 = [0, -Δcp, 0] # initial guess
    sol = root(func, x0, args=(a, b)) # solve the non-linear equation set
    A, Jp, Js = sol.x # extract the solution
    K = ((1-Dp/Dc)*Jp + (1-Ds/Dc)*Js)/2 # as above, used in next
    Δφ = (K/b)*ln(1 + b/a) # liquid junction potential

    # Approximate solution to the NP equations

    σII  = (Dc+Dp)*cpII + (Dc+Ds)*cs
    σI = (Dc+Dp)*cpI + (Dc+Ds)*cs
    Δσ = (Dc+Dp)*Δcp
    Δg = (Dc-Dp)*Δcp
    Δφ_approx = -Δg/Δσ * ln(σII/σI)

    results.append((cs, Jp, Js, Δφ, Δφ_approx))

cols =  ['cs', 'Jp/Dp', 'Js/Ds', 'Δφ', 'Δφ (approx)']
df = pd.DataFrame(results, columns=cols)

if args.show:

    import matplotlib.pyplot as plt
    df.set_index('cs', inplace=True)
    df[['Δφ', 'Δφ (approx)']].plot(logx=True, style=['-', '-.'], color=['k', 'k'])
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
