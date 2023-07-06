#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from scipy.io import FortranFile as fortran_file
from scipy.ndimage import laplace as del2
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("frame", nargs='?', help="fortran data file")
parser.add_argument('--rc', action='store', help='run control file')
parser.add_argument('--Lx', action='store', default=100.0, type=float, help='span in microns, default 100.0')
parser.add_argument('--Ly', action='store', default=100.0, type=float, help='span in microns, default 100.0')
parser.add_argument('--phi-max', action='store', default=70.0, type=float, help='max expected phi,default 70.0')
parser.add_argument('--zeta', action='store', default=-50.0, type=float, help='zeta potential, default -50.0')
parser.add_argument('--zoom', action='store', default=1.0, type=float, help='zoom factor, default 1.0')
parser.add_argument('--nn', action='store', default=5, type=int, help='spacing for quivers, default 5')
parser.add_argument('--drift', action='store_true', help='plot phoresis drift data')
parser.add_argument('--detail', action='store_true', help='plot conductivity and g fields')
parser.add_argument('--show', action='store_true', help='show the plot')
parser.add_argument('--save', action='store_true', help='save the plot')
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
    fp.readline() ; fp.readline() ; fp.readline()
    init_file = fp.readline().strip()
    fingerprint = fp.readline().strip()

print('run_control: run_name = ', run_name)
print('run_control: ions = ', ions)
print('run_control: nx, ny, ncmp = ', nx, ny, ncmp)
print('run_control: boundary_conditions = ', bc)
print('run_control: z = ', z)
print('run_control: d = ', d)
print('run_control: init_file = ', init_file)
print('run_control: fingerprint = ', fingerprint)

deltax = args.Lx / (nx - 1)
print('deltax = ', deltax)

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

rhotot = np.zeros((nx, ny))
for k in range(ncmp): rhotot += 0.25*rho[k, :, :]

rhox = 0.5*(rho[0, :, :] + rho[1, :, :])
rhoy = 0.5*(rho[2, :, :] + rho[3, :, :])

if fingerprint:
    with open(fingerprint) as fp:
        i0, j0, j1 = [int(v) for v in fp.readline().split()]
    print('fingerprint:', fingerprint)
    print('fingerprint:', i0, j0, j1)
    print('fingerprint: sigma = ', sigma[i0, j0:j1+1])
    print('fingerprint:     g = ', g[i0, j0:j1+1])
    print('fingerprint:   phi = ', phi[i0, j0:j1+1])

gy, gx = np.gradient(g/deltax)
ey, ex = np.gradient(-phi/deltax)

ix, iy = [np.multiply(sigma, e) - gradg for e, gradg in zip([ex, ey], [gx, gy])]
ixy, ixx = np.gradient(ix/deltax)
iyy, iyx = np.gradient(iy/deltax)

lnsigmay, lnsigmax = np.gradient(np.log(sigma)/deltax)
approx_curl = np.multiply(lnsigmax, gy) - np.multiply(lnsigmay, gx)

sumrhoz = 0.0
for k in range(ncmp): sumrhoz += z[k] * np.sum(rho[k, :, :])

print('diagnostic: (fortran)     |I| =', np.sqrt(np.sum(fix**2 + fiy**2))/deltax)
print('diagnostic: (python)      |I| =', np.sqrt(np.sum(ix**2 + iy**2)))
print('diagnostic: (fortran)   |∇·I| =', np.sqrt(np.sum(div**2))/deltax**2)
print('diagnostic: (python)    |∇·I| =', np.sqrt(np.sum((ixx + iyy)**2)))
print('diagnostic: (fortran) |∇ ⨯ I| =', np.sqrt(np.sum(curl**2))/deltax**2)
print('diagnostic: (python)  |∇ ⨯ I| =', np.sqrt(np.sum((iyx - ixy)**2)))
print('diagnostic: (estd)    |∇ ⨯ I| =', np.sqrt(np.sum(approx_curl**2)))
print('diagnostic: (python) Σi zi ρi =', sumrhoz)

Lx, Ly = nx - 1.0, ny - 1.0

x = np.linspace(0.0, Lx, nx)
y = np.linspace(0.0, Ly, ny)

phi -= np.min(phi)

kB, T, e = 1.38e-23, 298.0, 1.602e-19
epsr, eps0, eta = 78.0, 8.854e-12, 0.89e-3
kTbye = kB*T/e
prefac = epsr*eps0/eta * kTbye**2 * 1e9 # in um^2 / ms
ezetabykT = args.zeta*1e-3 / kTbye # recall that zeta in mV

A = 4.0*prefac*np.log(np.cosh(0.25*ezetabykT))
B = prefac*ezetabykT

print('kTbye [mV] = %0.2f' % (kTbye*1e3))
print('prefac [um^2 / ms] = %0.3f' % prefac)
print('drift coefficients A, B [um^2/ms] = ', A, B)

# np.gradient returns a list, which is stacked in reverse order

grad_g = np.stack(np.gradient(g)[::-1])
E = np.stack(np.gradient(-phi)[::-1])
grad_ln_rho = np.stack(np.gradient(np.log(rhotot))[::-1])
double_sigma = np.stack([sigma, sigma]) # double up for vector operations

I = double_sigma * E - grad_g

a = A * grad_ln_rho
b = B * grad_g / double_sigma
c = B * I / double_sigma
u = a + b + c

c_mod = np.sqrt(c[0, :, :]**2 + c[1, :, :]**2)
u_mod = np.sqrt(u[0, :, :]**2 + u[1, :, :]**2)

cbyu = c_mod / u_mod
cbyumin, cbyumax = np.min(cbyu), np.max(cbyu)
print('minmax: |c|/|u| = ', cbyumin, cbyumax)

phi *= kTbye * 1e3 # mV

spc = del2(phi) # space charge (a.u)

phimin, phimax = np.min(phi), np.max(phi)
print('minmax: phi = ', phimin, phimax)
phi_levels = list(np.linspace(0.0, args.phi_max, 10))
print('phi levels', phi_levels)

curl *= 1e3 # rescale

curlmin, curlmax = np.min(curl), np.max(curl)
print('minmax: curl = ', curlmin, curlmax)
curl_levels = list(np.linspace(-5.0, 5.0, 7))

rmin, rmax = np.min(rhotot), np.max(rhotot)
print('minmax: rhotot = ', rmin, rmax)
r_levels = list(np.linspace(0.0, 1.0, 11))

sigmin, sigmax = np.min(sigma), np.max(sigma)
print('minmax: conductivity = ', sigmin, sigmax)
sig_levels = list(np.linspace(0.0, 16.0, 16))
print('sigma levels', sig_levels)

gmin, gmax = np.min(g), np.max(g)
print('minmax: g = ', gmin, gmax)
g_levels = list(np.linspace(-2.0, 0.0, 16))
print('g levels', g_levels)

spcmin, spcmax = np.min(spc), np.max(spc)
print('minmax: space charge = ', spcmin, spcmax)
spclim = max([-spcmin, spcmax])
#sig_levels = list(np.linspace(0.0, 16.0, 16))
#print('sigma levels', sig_levels)

for k in range(ncmp):
    rhomin, rhomax = np.min(rho[k, :, :]), np.max(rho[k, :, :])
    print('minmax: rho', ions[k], '=', rhomin, rhomax)

rho_levels = list(np.linspace(0.0, 1.0, 11))
print('rho levels', rho_levels)

if not args.show and not args.save: exit()

if args.drift or args.detail:

    plt.subplot(2, 3, 1)
    plt.contour(x, y, phi, phi_levels, colors='black')
    plt.imshow(spc, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='RdBu',
               vmin=-spclim, vmax=spclim, alpha=0.5)
    if args.zoom > 1.0:
        zf = 1.0 / args.zoom
        plt.axis([0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx, 0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx])
    plt.title('φ and ∇²φ')

    plt.subplot(2, 3, 2)
    streamlines = plt.streamplot(x, y, ix, iy, color='black')
    plt.imshow(curl, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='cool',
               vmin=curl_levels[0], vmax=curl_levels[-1], alpha=0.5)
    # plt.colorbar()
    if args.zoom > 1.0:
        zf = 1.0 / args.zoom
        plt.axis([0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx, 0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx])
    plt.title('I and ∇×I')

    if args.detail:
        plt.subplot(2, 3, 3)
        # plt.contour(x, y, rhox, r_levels, colors='black')
        print('r_levels =', r_levels)
        plt.imshow(rhox, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='viridis',
                   vmin=r_levels[0], vmax=r_levels[-1], alpha=0.5)
        if args.zoom > 1.0:
            zf = 1.0 / args.zoom
            plt.axis([0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx, 0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx])
        plt.plot(x, 20+20*rhox[100, :], 'k') 
        plt.title(' : '.join(ions[0:2]))
    else:
        plt.subplot(2, 3, 3)
        #plt.contour(x, y, rhotot, r_levels, colors='black')
        plt.imshow(rhotot, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='viridis',
                   vmin=r_levels[0], vmax=r_levels[-1], alpha=0.5)
        plt.title('ρ_tot')

    plt.subplot(2, 3, 4)
    plt.contour(x, y, sigma, sig_levels, colors='black')
    plt.imshow(sigma, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='viridis',
               vmin=sig_levels[0], vmax=sig_levels[-1], alpha=0.5)
    if args.zoom > 1.0:
        zf = 1.0 / args.zoom
        plt.axis([0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx, 0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx])
    plt.title('σ (conductivity)')

    if args.detail:
        plt.subplot(2, 3, 5)
        plt.contour(x, y, g, g_levels, colors='black')
        plt.imshow(g, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='viridis',
                   vmin=g_levels[0], vmax=g_levels[-1], alpha=0.5)
        if args.zoom > 1.0:
            zf = 1.0 / args.zoom
            plt.axis([0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx, 0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx])
        plt.title('g field')
    else:
        plt.subplot(2, 3, 5)
        quivers1 = plt.quiver(x[::4], y[::4], a[0, ::4, ::4], a[1, ::4, ::4], angles='xy', scale_units='xy', scale=0.001, color='red')
        quivers2 = plt.quiver(x[::4], y[::4], b[0, ::4, ::4], b[1, ::4, ::4], angles='xy', scale_units='xy', scale=0.001, color='blue')
        quivers3 = plt.quiver(x[::4], y[::4], c[0, ::4, ::4], c[1, ::4, ::4], angles='xy', scale_units='xy', scale=0.001, color='yellow')
        quivers4 = plt.quiver(x[::4], y[::4], u[0, ::4, ::4], u[1, ::4, ::4], angles='xy', scale_units='xy', scale=0.001, color='black')
        plt.imshow(100*c_mod/u_mod, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='cool',
                   vmin=0.0, vmax=30.0, alpha=0.5)
        plt.title('non local fraction')

    if args.detail:
        plt.subplot(2, 3, 6)
        #plt.contour(x, y, rhoy, r_levels, colors='black')
        plt.imshow(rhoy, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='viridis',
                   vmin=r_levels[0], vmax=r_levels[-1], alpha=0.5)
        if args.zoom > 1.0:
            zf = 1.0 / args.zoom
            plt.axis([0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx, 0.5*(1-zf)*Lx, 0.5*(1+zf)*Lx])
        plt.plot(x, 20+20*rhoy[100, :], 'k') 
        plt.title(' : '.join(ions[2:4]))
    else:
        plt.subplot(2, 3, 6)
        quivers = plt.quiver(x[::2], y[::2], ex[::2, ::2], ey[::2, ::2])
        plt.contour(x, y, phi, phi_levels, colors='black')
        plt.imshow(phi, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='plasma',
                   vmin=phi_levels[0], vmax=phi_levels[-1], alpha=0.5)
        plt.title('φ and −∇φ')

else: # not args.drift or args.detail

    plt.subplot(2, 3, 1)
    quivers = plt.quiver(x[::2], y[::2], ex[::2, ::2], ey[::2, ::2])
    plt.contour(x, y, phi, phi_levels, colors='black')
    plt.imshow(phi, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='plasma',
               vmin=phi_levels[0], vmax=phi_levels[-1], alpha=0.5)
    # plt.colorbar()
    plt.title('phi and E')

    plt.subplot(2, 3, 4)
    streamlines = plt.streamplot(x, y, ix, iy, color='black')
    plt.imshow(curl, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='cool',
               vmin=curl_levels[0], vmax=curl_levels[-1], alpha=0.5)
    # plt.colorbar()
    plt.title('I and curl I')

    for k in range(ncmp):
        trho = rho[k, :, :].copy()
        plt.subplot(2, 3, 2+k if k < 2 else 3+k)
        plt.contour(x, y, trho, rho_levels, colors='black')
        plt.imshow(trho, extent=[0.0, Lx, 0.0, Ly], origin='lower',
                   cmap='Blues', vmin=rho_levels[0], vmax=rho_levels[-1], alpha=0.5)
        plt.title(ions[k])

    if ncmp == 2:

        beta = np.sum(d*np.sign(z)) / np.sum(d*np.abs(z))
        print('binary: beta = %0.3f' % beta)

        rhos = - rho[0, :, :] / z[1]
        phiplus = phi + beta * np.log(rhos) * kTbye * 1e3 # mV

        plt.subplot(2, 3, 5)
        plt.contour(x, y, rhos, rho_levels, colors='black')
        plt.imshow(rhos, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='Blues',
                   vmin=rho_levels[0], vmax=rho_levels[-1], alpha=0.5)
        plt.title('rho_s')

        plt.subplot(2, 3, 6)
        plt.contour(x, y, phiplus, phi_levels, colors='black')
        plt.imshow(phiplus, extent=[0.0, Lx, 0.0, Ly], origin='lower', cmap='plasma',
                   vmin=phi_levels[0], vmax=phi_levels[-1], alpha=0.5)
        plt.title('phi + beta ln rho_s')

if args.show:
    plt.show()

if args.save:
    png_file = data_file + '.png'
    plt.savefig(png_file)
    print('written', png_file)
