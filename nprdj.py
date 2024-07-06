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

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

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

# Since cp + cs is a linear in x, the above ODEs can be solved
# exactly.  This leaves only a non-linear equation set for the
# constants of integration.

# Can run calculations at cs = 0 but need to make sure that cp is not
# also zero as a boundary condition, and vice versa.  Under such
# conditions the approximate solution becomes exact.

# The --cpoly and --csalt arguments can be specified as either a pair
# of concentrations as cI:cII, or as a single value cI = cII.

import argparse
import numpy as np
from scipy.optimize import root
from numpy import log as ln

parser = argparse.ArgumentParser("Nernst-Planck steady-state")
parser.add_argument('-D', '--Darr', default='0.0,2.03,1.33', help='diffusion coeffs, default 0.0,2.03,1.33')
parser.add_argument('-p', '--cpoly', default='1.0:0.0', help='polymer concentrations, default 1.0:0.0')
parser.add_argument('-c', '--csalt', default='1.0', help='background salt, default 1.0')
parser.add_argument('-v', '--verbosity', action='count', default=0, help='increasing verbosity')
parser.add_argument('-a', '--approx', action='store_true', help='show the approximate soln too')
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

a = cpI + csI # the constant in cp + cs = a + b x
b = Δcp + Δcs # the gradient in the same

def func(x):
    A, Jp, Js = x # extract variables
    K = ((1-Dp/Dc)*Jp + (1-Ds/Dc)*Js)/2 # a derived quantity
    cp0 = Jp/(Jp+Js)*a + A * a**(K/b) # cp at x = 0
    cp1 = Jp/(Jp+Js)*(a+b) + A * (a+b)**(K/b) # cp at x = L
    return([cp0-cpI, cp1-cpII, b+Jp+Js-K]) # bcs ; constraint on K

x0 = [0, -Δcp, -Δcs] # initial guess
sol = root(func, x0) # solve the non-linear equation set
A, Jp, Js = sol.x # extract the solution
K = ((1-Dp/Dc)*Jp + (1-Ds/Dc)*Js)/2 # as above, used below
Δφ = (K/b)*ln(1 + b/a) # liquid junction potential

x = np.linspace(0, 1, 41)
cp = (a+b*x)*Jp/(Jp+Js) + A * (a+b*x)**(K/b) # as above
cs = (a+b*x)*Js/(Jp+Js) - A * (a+b*x)**(K/b) # as above
φ = (K/b)*ln(1 + b*x/a) # the solution for the potential

# Approximate solution to the NP equations

σ  = (Dc+Dp)*(cpI+Δcp*x) + (Dc+Ds)*(csI+Δcs*x)
σII = (Dc+Dp)*cpII + (Dc+Ds)*csII
σI = (Dc+Dp)*cpI + (Dc+Ds)*csI
Δσ = (Dc+Dp)*Δcp + (Dc+Ds)*Δcs
Δg = (Dc-Dp)*Δcp + (Dc-Ds)*Δcs 
φ_approx = -Δg/Δσ * ln(σ/σI)
Δφ_approx = -Δg/Δσ * ln(σII/σI)

if args.show:

    import matplotlib.pyplot as plt
    plt.plot(x, φ-φ[0], 'k-', label='φ')
    plt.plot(x, cp, 'g-', label='cp')
    plt.plot(x, cs, 'b-', label='cs')
    if args.total:
        plt.plot(x, cp+cs, 'r-', label='cp + cs')
    if args.approx:
        plt.plot([0, 1], [cpI, cpII], 'g--', label='cp (linear)')
        plt.plot([0, 1], [csI, csII], 'b--', label='cs (linear)')
        plt.plot(x, φ_approx-φ_approx[0], 'k-.', label='φ (approx)')
    if args.legend:
        plt.legend()
    plt.xlabel("x / L") ; plt.show()

elif args.output:

    import pandas as pd
    df = pd.DataFrame({'x':x, 'phi':φ-φ[0], 'cp':cp, 'cs':cs, 'cp+cs':(cp+cs),
                       'phi(approx)':φ_approx-φ_approx[0]})
    headers = [f'{col}[{i+1}]' for i, col in enumerate(df.columns)]
    header_row = '#  ' + '  '.join(headers)
    data_rows = df.to_string(index=False).split('\n')[1:]
    with open(args.output, 'w') as f:
        print('\n'.join([header_row] + data_rows), file=f)
    print('Data', ', '.join(headers), 'written to', args.output)

else:

    print('cpI, cpII =', cpI, cpII)
    print('csI, csII =', csI, csII)
    print('Dp, Ds, Dc =', Dp, Ds, Dc)

    print('Jp/Dp, Js/Ds, A =', Jp, '\t', Js, '\t', A)
    print('Δφ           =\t', Δφ)
    print('Δφ (approx)  =\t', Δφ_approx)
