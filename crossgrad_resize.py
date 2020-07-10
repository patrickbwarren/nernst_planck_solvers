#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os.path
import argparse
import numpy as np
from scipy.io import FortranFile as fortran_file
from scipy.interpolate import RectBivariateSpline
from crossgrad import crossgrad as cg

parser = argparse.ArgumentParser()
parser.add_argument("name", action='store', help="run name")
parser.add_argument("--frame", help="fortran data file")
parser.add_argument('--rc', action='store', help='run control file')
parser.add_argument("--finger", action='store', default="[7, 11, 14]", help="fingerprint controls")
parser.add_argument('--nside', action='store', default=100, type=int, help='length of square side')
parser.add_argument('--njacobi', action='store', default=5000, type=int, help='number of jacobi sweeps')
parser.add_argument('--niter', action='store', default=20, type=int, help='number of iterations')
parser.add_argument('--nstep', action='store', default=200, type=int, help='number of iterations')
parser.add_argument('--ndiag', action='store', default=10, type=int, help='number of iterations')
parser.add_argument('--nsave', action='store', default=20, type=int, help='number of iterations')
parser.add_argument('--nmon', action='store', default=0, type=int, help='number of iterations')
parser.add_argument('--nmon2', action='store', default=1, type=int, help='number of iterations')
parser.add_argument('--dt', action='store', default=0.1, type=float, help='time step')
parser.add_argument('--run', action='store_true', help='run the code')
args = parser.parse_args()

if not args.rc:
    raise Exception('--rc run control file required')
    
with open(args.rc) as fp:
    run_name = fp.readline().split()[0]
    ions = fp.readline().split()
    nx, ny, ncmp = [int(v) for v in fp.readline().split()[0:3]]
    bc = [int(v) for v in fp.readline().split()[0:8]]
    z = [int(v) for v in fp.readline().split()[0:ncmp]]
    d = [float(v) for v in fp.readline().split()[0:ncmp]]
    for skip in range(3): fp.readline()
    init_file = fp.readline().strip()
    fingerprint = fp.readline().strip()

print('run_control: nx, ny, ncmp = ', nx, ny, ncmp)
print('run_control: boundary_conditions = ', bc)
print('run_control: z = ', z)
print('run_control: d = ', d)

data_file = args.frame if args.frame else init_file

print('reading data from', data_file)

with fortran_file(data_file) as f:
    r = f.read_reals()

rr = r.reshape((ncmp+5, nx*ny))

rho = np.zeros((ncmp, nx, ny))

for k in range(ncmp):
    rho[k, :, :] = rr[k, :].reshape((nx, ny)).transpose().copy()

phi = rr[ncmp, :].reshape((nx, ny)).transpose().copy()
fix = rr[ncmp+1, :].reshape((nx, ny)).transpose().copy()
fiy = rr[ncmp+2, :].reshape((nx, ny)).transpose().copy()
div = rr[ncmp+3, :].reshape((nx, ny)).transpose().copy()
curl = rr[ncmp+4, :].reshape((nx, ny)).transpose().copy()

g = np.zeros((nx, ny))
for k in range(ncmp): g += z[k] * d[k] * rho[k, :, :]

sigma = np.zeros((nx, ny))
for k in range(ncmp): sigma += z[k]**2 * d[k] * rho[k, :, :]

if fingerprint:
    with open(fingerprint) as fp:
        i0, j0, j1 = [int(v) for v in fp.readline().split()]
    print('fingerprint:', fingerprint)
    print('fingerprint:', i0, j0, j1)
    print('fingerprint: sigma = ', sigma[i0, j0:j1+1])
    print('fingerprint:     g = ', g[i0, j0:j1+1])
    print('fingerprint:   phi = ', phi[i0, j0:j1+1])

gy, gx = np.gradient(g)
ey, ex = np.gradient(-phi)
ix, iy = [np.multiply(sigma, e) - gradg for e, gradg in zip([ex, ey], [gx, gy])]
ixy, ixx = np.gradient(ix)
iyy, iyx = np.gradient(iy)

lnsigmay, lnsigmax = np.gradient(np.log(sigma))
approx_curl = np.multiply(lnsigmax, gy) - np.multiply(lnsigmay, gx)

sumrhoz = 0.0
for k in range(ncmp): sumrhoz += z[k] * np.sum(rho[k, :, :])

print('diagnostic: (fortran)     |I| =', np.sqrt(np.sum(fix**2 + fiy**2)))
print('diagnostic: (python)      |I| =', np.sqrt(np.sum(ix**2 + iy**2)))
print('diagnostic: (fortran)   |∇·I| =', np.sqrt(np.sum(div**2)))
print('diagnostic: (python)    |∇·I| =', np.sqrt(np.sum((ixx + iyy)**2)))
print('diagnostic: (fortran) |∇ ⨯ I| =', np.sqrt(np.sum(curl**2)))
print('diagnostic: (python)  |∇ ⨯ I| =', np.sqrt(np.sum((iyx - ixy)**2)))
print('diagnostic: (estd)    |∇ ⨯ I| =', np.sqrt(np.sum(approx_curl**2)))
print('diagnostic: (python) Σi zi ρi =', sumrhoz)

# resize the density data

nx2, ny2 = args.nside+1, args.nside+1

Lx, Ly = nx2 - 1.0, ny2 - 1.0
x = np.linspace(0.0, Lx, nx)
y = np.linspace(0.0, Ly, ny)
x2 = np.linspace(0.0, Lx, nx2)
y2 = np.linspace(0.0, Ly, ny2)
 
rho2 = np.zeros((ncmp, nx2, ny2))

for k in range(ncmp):
    spline = RectBivariateSpline(y, x, rho[k, :, :])
    rho2[k, :, :] = spline(y2, x2)

# create and save the revised data and run controls

rc_file = args.name + '.rc'

with open(rc_file, 'w') as fp:
    fp.write(args.name + '\n')
    fp.write(' '.join(ions) + '\n')
    fp.write('%i %i %i' % (nx2, ny2, ncmp) + ' :: nx ny ncmp\n')
    fp.write(' '.join(['%i' % val for val in bc]) + ' :: boundary_conditions\n')
    fp.write('\t'.join(['%i' % val for val in z]) + ' :: z\n')
    fp.write('\t'.join(['%g' % val for val in d]) + ' :: d\n')
    fp.write('%g :: dt\n' % args.dt)

exit() if not args.run else None

cg.read_rc(rc_file) # first call to FORTRAN sector

for k in range(ncmp):
    cg.rho[:, :, k] = rho2[k, :, :]

cg.converge(args.njacobi, args.niter, args.nmon2)

if not cg.converged:
    raise Exception('not converged')

cg.diagnostics()

data_file = args.name + '_init.dat'
fingerprint = data_file + '.fingerprint'

cg.write_array(data_file)

with open(rc_file, 'a') as fp:
    fp.write('%i %i %i' % (args.njacobi, args.niter, args.nmon) + ' :: njacobi niter nmon\n')
    fp.write('%i %i %i' % (args.nstep, args.ndiag, args.nsave) + ' :: nstep ndiag nsave\n')
    fp.write(data_file + '\n')
    fp.write(fingerprint + '\n')

cg.write_fingerprint(fingerprint, *eval(args.finger))

print('finalised run control', rc_file)
