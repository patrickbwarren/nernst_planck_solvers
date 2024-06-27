#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import argparse
import numpy as np
from numpy import exp
from scipy.special import erf
from crossgrad import crossgrad as cg

parser = argparse.ArgumentParser()
parser.add_argument("name", action='store', help="run name")
parser.add_argument("--ionx", action='store', default="Na1Cl1", help="list of ions")
parser.add_argument("--iony", action='store', default="K1Br1", help="list of ions")
parser.add_argument("--ionb", action='store', default="Li1OH1", help="list of ions")
parser.add_argument("--cx", action='store', default="0.1,1", help="ion gradient in x")
parser.add_argument("--cy", action='store', default="0.1,1", help="ion gradient in y")
parser.add_argument("--back", action='store', default=None, help="background ion concentration")
parser.add_argument("--shape", action='store', default="10,25,5", help="shape of gradient in x")
parser.add_argument("--finger", action='store', default="7,11,14", help="fingerprint controls")
parser.add_argument('--nx', action='store', default=20, type=int, help='length of side')
parser.add_argument('--ny', action='store', default=20, type=int, help='length of side')
parser.add_argument('--njacobi', action='store', default=5000, type=int, help='number of jacobi sweeps')
parser.add_argument('--niter', action='store', default=20, type=int, help='number of iterations')
parser.add_argument('--nmon', action='store', default=1, type=int, help='number of iterations')
parser.add_argument('--nstep', action='store', default='200', help='number of iterations')
parser.add_argument('--ndiag', action='store', default=10, type=int, help='number of iterations')
parser.add_argument('--nsave', action='store', default=20, type=int, help='number of iterations')
parser.add_argument('--dt', action='store', default=0.1, type=float, help='time step')
parser.add_argument('--gaussians', action='store_true', help='use Gaussians')
parser.add_argument('--solve', action='store_true', help='run the solver')
parser.add_argument('--trace', action='store_true', help='extract trace')
args = parser.parse_args()

property = {'d':0, 'z':1, 'name':2}

ion_props = { # base units here are μm^2/ms
    'H'   : (9.31, 1, 'H+'),   # proton (H+)
    'K'   : (1.96, 1, 'K+'),   # cation (K+)
    'NHf' : (1.96, 1, 'NH4+'), # cation (NH4+)
    'Na'  : (1.33, 1, 'Na+'),  # cation (Na+)
    'Li'  : (1.03, 1, 'Li+'),  # cation (Li+)
    'Pb'  : (0.94, 2, 'Pb2+'), # cation lead (2+)
    'Me'  : (0.26, 1, 'Me+'),  # methylene blue
    'OH'  : (5.28, -1, 'OH-'),   # anion (OH-)
    'Cl'  : (2.03, -1, 'Cl-'),   # anion (Cl-)
    'Br'  : (2.01, -1, 'Br-'),   # anion (Br-)
    'COt' : (0.95, -2, 'CO32-'), # anion carbonate (2-)
    'Fl'  : (0.63, -2, 'Fl2-'),  # anion fluorescein (2-)
    'SDS' : (0.28, -1, 'SDS-'),  # anion SDS (-)
}

def ion_split(s):
    '''split a string into ion / number pairs and return a dict of the same'''
    l = filter(None, re.split(r'(\d+)', s)) ; it = iter(l)
    return dict([(k, int(v)) for k, v in zip(it, it)])

ionxd = ion_split(args.ionx)
ionyd = ion_split(args.iony)
ionbd = ion_split(args.ionb) if args.back else {}

# Extract a list of ions from the keys of a merged dictionary

iond = {**ionxd, **ionyd, **ionbd}
ions = list(iond.keys())

ionx = [int(ionxd[ion]) if ion in ionxd else 0 for ion in ions]
iony = [int(ionyd[ion]) if ion in ionyd else 0 for ion in ions]

if args.back:
    ionb = [int(ionbd[ion]) if ion in ionbd else 0 for ion in ions]
    back = eval(args.back)
else:
    ionb = [0] * len(ions)

d = np.array([ion_props[ion][property['d']] for ion in ions])
z = np.array([ion_props[ion][property['z']] for ion in ions])
ion_names = [ion_props[ion][property['name']] for ion in ions]

print('ions = ', ion_names)
print('z = ', z)
print('d = ', d)

ion_list = zip(['x', 'y', 'back'], [ionx, iony, ionb]) if args.back else zip(['x', 'y'], [ionx, iony])

for what, ion in ion_list:
    beta = np.sum(ion*d*np.sign(z)) / np.sum(ion*d*np.abs(z))
    salt = ':'.join([ion_names[i] for i in range(len(ion)) if ion[i]])
    dd = [d[i] for i in range(len(ion)) if ion[i]]
    print('%5s %8s : beta = %0.3f' % (what, salt, beta), '; D =', dd)

bc = [0] * 8 # default is all zero flux

rc_file = args.name + '.rc'

with open(rc_file, 'w') as fp:
    fp.write(args.name + '\n')
    fp.write(' '.join(ion_names) + '\n')
    fp.write('%i %i %i' % (1+args.nx, 1+args.ny, len(ions)) + ' :: nx ny ncmp\n')
    fp.write(' '.join(['%i' % val for val in bc]) + ' :: boundary_conditions\n')
    fp.write('\t'.join(['%i' % val for val in z]) + ' :: z\n')
    fp.write('\t'.join(['%g' % val for val in d]) + ' :: d\n')
    fp.write('%g :: dt\n' % args.dt)

if args.solve:

    cg.read_rc(rc_file) # first call to FORTRAN sector

    dx, dy = args.nx / (cg.nx-1), args.ny / (cg.ny-1)

    print('dx, dy =', dx, dy)

    x, y = np.arange(cg.nx)*dx, np.arange(cg.ny)*dy
    xm, ym = np.mean(x), np.mean(y)
    xx, yy = np.meshgrid(x, y)
axes
    d, R, w = [v for v in eval(f'[{args.shape}]')]

    if args.gaussians:

        c0, c1 = eval(f'[{args.cx}]')
        rr2 = (xx-xm+d)**2+(yy-ym)**2
        cx = c0 + (c1-c0)*exp(-rr2/(2*R**2))
        print(' first Gaussian d, R =', d, R, ' c0, c1 =', c0, c1)

        c0, c1 = eval(f'[{args.cy}]')
        rr2 = (xx-xm-d)**2+(yy-ym)**2
        cy = c0 + (c1-c0)*exp(-rr2/(2*R**2))
        print('second Gaussian d, R =', d, R, ' c0, c1 =', c0, c1)

    else:

        c0, c1 = eval(f'[{args.cx}]')
        rr = np.sqrt((xx-xm+d)**2+(yy-ym)**2)
        cx = 0.5*(c0 + c1) - 0.5*(c1 - c0)*erf((rr-R)/w)
        print(' first disc d, R, w =', d, R, w, ' c0, c1 =', c0, c1)

        c0, c1 = eval(f'[{args.cy}]')
        rr = np.sqrt((xx-xm-d)**2+(yy-ym)**2)
        cy = 0.5*(c0 + c1) - 0.5*(c1 - c0)*erf((rr-R)/w)
        print('second disc d, R, w =', d, R, w, ' c0, c1 =', c0, c1)

    if args.back:
        back = eval(args.back) * np.ones_like(cx)
    else:
        back = np.zeros_like(cx)

    for k in range(cg.ncmp):
        cg.rho[:, :, k] = 0.0
        if ionx[k]: cg.rho[:, :, k] += cx[:, :] * ionx[k]
        if iony[k]: cg.rho[:, :, k] += cy[:, :] * iony[k]
        if ionb[k]: cg.rho[:, :, k] += back[:, :] * ionb[k]

    for k in range(cg.ncmp):
        rhomin = np.min(cg.rho[:, :, k])
        rhomax = np.max(cg.rho[:, :, k])
        print(ion_names[k], ' min/max = ', rhomin, rhomax)

    if args.trace:
        for i in range(min(cg.nx, cg.ny)):
            print(i, ':', ' '.join([str(cg.rho[i, i, k]) for k in range(len(ions))]))
    else:
        cg.converge(args.njacobi, args.niter, args.nmon)

    if not cg.converged:
        raise Exception('not converged')

    cg.diagnostics()

    data_file = args.name + '_init.dat'
    fingerprint = data_file + '.fingerprint'

    cg.write_array(data_file)

    with open(rc_file, 'a') as fp:
        fp.write('%i %i %i' % (args.njacobi, args.niter, args.nmon) + ' :: njacobi niter nmon\n')
        fp.write('%i %i %i' % (eval(args.nstep), args.ndiag, args.nsave) + ' :: nstep ndiag nsave\n')
        fp.write(data_file + '\n')
        fp.write(fingerprint + '\n')

    cg.write_fingerprint(fingerprint, *eval(f'[{args.finger}]'))

print('finalised run control', rc_file)
