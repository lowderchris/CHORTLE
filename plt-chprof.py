# plt-chprof.py
# Script to generate selected time-latitude profiles of computed coronal hole quantities

import sunpy.visualization.colormaps
import sunpy.coordinates
import scipy.signal
import configparser

import sunpy.io

import numpy as np

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

# Read configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

# Specify directory structures from configuration file
outdir = config['paths']['outdir']

# Load processed coronal hole data files
chprof0 = np.load(outdir + 'dat/chprof0.npy')
chprof = np.load(outdir + 'dat/chprof.npy')
improf = np.load(outdir + 'dat/improf.npy')
sfxprof = np.load(outdir + 'dat/sfxprof.npy')
ufxprof = np.load(outdir + 'dat/ufxprof.npy')
crs = np.load(outdir + 'dat/crs.npy')

# Generate array of time objects
dts = sunpy.coordinates.sun.carrington_rotation_time(crs)
dts.format = 'datetime'
dts = dts.value

ndts = mdates.date2num(dts)

# Begin plotting

# Coronal hole profile plot (fractional coverage)
f, (ax) = plt.subplots(1, figsize=[6,2.5])
im = ax.imshow(chprof0, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', vmin=0, cmap='Greys_r')
cb = plt.colorbar(im, ax=ax, label='Fractional CH Coverage', extend='max')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
plt.tight_layout()
plt.savefig(outdir + 'plt/chprof0.pdf')
plt.savefig(outdir + 'plt/chprof0.png', dpi=300)

# Coronal hole profile plot (coverage)
f, (ax) = plt.subplots(1, figsize=[6,2.5])
im = ax.imshow(chprof, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap='Greys_r')
cb = plt.colorbar(im, ax=ax, label='CH Coverage')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
plt.tight_layout()
plt.savefig(outdir + 'plt/chprof.pdf')
plt.savefig(outdir + 'plt/chprof.png', dpi=300)

# Input EUV data profile plot
f, (ax) = plt.subplots(1, figsize=[6,2.5])
im = ax.imshow(improf, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap=sunpy.visualization.colormaps.cm.sdoaia193)
cb = plt.colorbar(im, ax=ax, label='Input EUV Data')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
plt.tight_layout()
plt.savefig(outdir + 'plt/improf.pdf')
plt.savefig(outdir + 'plt/improf.png', dpi=300)

# Coronal hole enclosed signed magnetic flux profile plot
f, (ax) = plt.subplots(1, figsize=[6,2.5])
im = ax.imshow(sfxprof, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap='RdBu_r', vmin=-1e20, vmax=1e20)
cb = plt.colorbar(im, ax=ax, label='CH Signed Flux [Mx]', extend='both')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
plt.tight_layout()
plt.savefig(outdir + 'plt/sfxprof.pdf')
plt.savefig(outdir + 'plt/sfxprof.png', dpi=300)

# Coronal hole enclosed unsigned magnetic flux profile plot
f, (ax) = plt.subplots(1, figsize=[6,2.5])
im = ax.imshow(ufxprof, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap='magma')
cb = plt.colorbar(im, ax=ax, label='CH Unsigned Flux [Mx]')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
plt.tight_layout()
plt.savefig(outdir + 'plt/ufxprof.pdf')
plt.savefig(outdir + 'plt/ufxprof.png', dpi=300)

# Coronal hole enclosed magnetic flux plot
f, (ax) = plt.subplots(1, figsize=[6,2.5])
ax.plot(dts, ufxprof[:,:].sum(0) / 1e23, 'k', label='Total')
ax.plot(dts, ufxprof[580:719,:].sum(0) / 1e23, 'C0', label='NP')
ax.plot(dts, ufxprof[141:579,:].sum(0) / 1e23, 'C4', label='EQ')
ax.plot(dts, ufxprof[0:140,:].sum(0) / 1e23, 'C3', label='SP')
axb = ax.twiny()
ax.set_xlim(dts[0], dts[-1])
axb.set_xlim(crs[0], crs[-1])
ax.set_ylim(0,1.25)
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Unsigned Magnetic Flux \n [1e23 Mx]')
plt.tight_layout()
ax.legend(bbox_to_anchor=(1.03,0),loc=3, ncol=1)
ax.set_position([0.12265046296296299, 0.23311111111111105, 0.685, 0.545])
plt.savefig(outdir + 'plt/flx.pdf')
plt.savefig(outdir + 'plt/flx.png', dpi=300)