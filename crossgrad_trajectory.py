#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from scipy.io import FortranFile as fortran_file
from scipy.interpolate import RectBivariateSpline as spline2d

parser = argparse.ArgumentParser()
parser.add_argument("rc", nargs='?', help='run control file')
parser.add_argument('--Lx', action='store', default=100.0, type=float, help='span in microns, default 100.0')
parser.add_argument('--Ly', action='store', default=100.0, type=float, help='span in microns, default 100.0')
parser.add_argument('--xarr', action='store', default='[5, 100, 10]', help='initial position array, default=[5, 100, 10]')
parser.add_argument('--yarr', action='store', default='[5, 100, 10]', help='initial position array, default=[5, 100, 10]')
parser.add_argument('--xy', action='store', default='[50, 50]', help='initial position, default [50, 50]')
parser.add_argument('--zeta', action='store', default=-50.0, type=float, help='zeta potential, default -50.0')
parser.add_argument('--single', action='store_true', help='do a single trajectory rather than a mesh')
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
    dt = float(fp.readline().split()[0])
    fp.readline()
    nstep, _, nsave = [int(v) for v in fp.readline().split()[0:3]]
    fp.readline()
    fingerprint = fp.readline().strip()

print('run_control: run_name = ', run_name)
print('run_control: ions = ', ions)
print('run_control: nx, ny, ncmp = ', nx, ny, ncmp)
print('run_control: boundary_conditions = ', bc)
print('run_control: z = ', z)
print('run_control: d = ', d)
print('run_control: dt = ', dt)
print('run_control: nstep, nsave = ', nstep, nsave)
print('run_control: fingerprint = ', fingerprint)

kB, T, e = 1.38e-23, 298.0, 1.602e-19
epsr, eps0, eta = 78.0, 8.854e-12, 0.89e-3 
kTbye = kB*T/e
prefac = epsr*eps0/eta * kTbye**2 * 1e9 # in um^2 / ms
ezetabykT = args.zeta*1e-3 / kTbye # zeta in mV
A = 4.0*prefac*np.log(np.cosh(ezetabykT/4.0))
B = prefac*ezetabykT

print('kTbye [mV] = %0.2f' % (kTbye*1e3))
print('prefac [um^2 / ms] = %0.3f' % prefac)
print('drift coefficients A, B [um^2/ms] = ', A, B)

deltax = args.Lx / (nx - 1)
deltat = dt * deltax**2
print('deltax, deltat = ', deltax, deltat)

nframes = nstep // nsave

x = np.linspace(0.0, args.Lx, nx)
y = np.linspace(0.0, args.Ly, ny)

if args.single:
    px, py = [np.array([v], dtype=np.float) for v in eval(args.xy)]
else:
    xarr, yarr = map(eval, [args.xarr, args.yarr])
    x0 = np.arange(*xarr, dtype=np.float)
    y0 = np.arange(*yarr, dtype=np.float)
    xx, yy = np.meshgrid(x0, y0)
    px, py = np.array([[x, y] for x, y in zip(xx.flatten(), yy.flatten())]).transpose()

for frame in range(nframes+1):

    data_file = '%s_f%05i.dat' % (run_name, frame)

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

    rhotot = np.zeros((nx, ny))
    for k in range(ncmp): rhotot += rho[k, :, :]

    if frame == 0:
        with open(fingerprint) as fp:
            i0, j0, j1 = [int(v) for v in fp.readline().split()]
        print('fingerprint:', fingerprint)
        print('fingerprint:', i0, j0, j1)
        print('fingerprint: sigma = ', sigma[i0, j0:j1+1])
        print('fingerprint:     g = ', g[i0, j0:j1+1])
        print('fingerprint:   phi = ', phi[i0, j0:j1+1])

    gy, gx = np.gradient(g)
    ey, ex = [-v for v in np.gradient(phi)]
    lnry, lnrx = np.gradient(np.log(rhotot))

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

    lnrx_spline = spline2d(x, y, lnrx)
    lnry_spline = spline2d(x, y, lnry)

    lnrpx = lnrx_spline(px, py, grid=False)
    lnrpy = lnry_spline(px, py, grid=False)

    ex_spline = spline2d(x, y, ex)
    ey_spline = spline2d(x, y, ey)
    
    epx = ex_spline(px, py, grid=False)
    epy = ey_spline(px, py, grid=False)

    if args.single:
        i, j = [np.rint(p/deltax).astype(np.int)[0] for p in [px, py]]
        print('spline2d: grad ln s(%0.2f, %0.2f) = %g %g' % (px[0], py[0], lnrpx[0], lnrpy[0]))
        print('spline2d: grad ln s(%2i, %2i)       = %g %g' % (i, j, lnrx[i, j], lnry[i, j]))
        print('spline2d: e(%0.2f, %0.2f) = %g %g' % (px[0], py[0], epx[0], epy[0]))
        print('spline2d: e(%2i, %2i)       = %g %g' % (i, j, ex[i, j], ey[i, j]))

    # correct for the gradients assuming unit lattice spacing
    
    ux = (A*lnrpx + B*epx) / deltax
    uy = (A*lnrpy + B*epy) / deltax

    t = frame * nsave * deltat

    dx = nsave * deltat * ux
    dy = nsave * deltat * uy

    for k in range(len(px)):
        print('%03i %10.3f %3i : px, py = %10.3f %10.3f (position)' % (frame, t, k, px[k], py[k]))
        print('%03i %10.3f %3i : ux, uy = %10.5f %10.5f (chemiphoresis)' % (frame, t, k, A*lnrpx[k], A*lnrpy[k]))
        print('%03i %10.3f %3i : ux, uy = %10.5f %10.5f (electrophoresis)' % (frame, t, k, B*epx[k], B*epy[k]))
        print('%03i %10.3f %3i : ux, uy = %10.5f %10.5f (total)' % (frame, t, k, ux[k], uy[k]))
        print('%03i %10.3f %3i : dx, dy = %10.5f %10.5f (displacement)' % (frame, t, k, dx[k], dy[k]))

    px += dx
    py += dy
