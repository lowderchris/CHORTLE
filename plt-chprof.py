import sunpy.visualization.colormaps
import sunpy.coordinates
import scipy.signal
import configparser

import sunpy.io

import numpy as np

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

chprof0 = np.load('dat/chprof0.npy')
chprof = np.load('dat/chprof.npy')
improf = np.load('dat/improf.npy')
sfxprof = np.load('dat/sfxprof.npy')
ufxprof = np.load('dat/ufxprof.npy')
crs = np.load('dat/crs.npy')

# Generate array of time objects
dts = sunpy.coordinates.sun.carrington_rotation_time(crs)
dts.format = 'datetime'
dts = dts.value

ndts = mdates.date2num(dts)

# 1920x1080 recommended size

f, (ax) = plt.subplots(1, figsize=[6,2.5])
im = ax.imshow(chprof0, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', vmin=0, cmap='viridis')
cb = plt.colorbar(im, ax=ax, label='Fractional CH Coverage', extend='max')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
plt.tight_layout()
plt.savefig('plt/chprof0.pdf')
plt.savefig('plt/chprof0.png', dpi=300)

f, (ax) = plt.subplots(1, figsize=[6,2.5])
im = ax.imshow(chprof, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap='viridis')
cb = plt.colorbar(im, ax=ax, label='CH Coverage')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
plt.tight_layout()
plt.savefig('plt/chprof.pdf')
plt.savefig('plt/chprof.png', dpi=300)

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
plt.savefig('plt/improf.pdf')
plt.savefig('plt/improf.png', dpi=300)

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
plt.savefig('plt/sfxprof.pdf')
plt.savefig('plt/sfxprof.png', dpi=300)

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
plt.savefig('plt/ufxprof.pdf')
plt.savefig('plt/ufxprof.png', dpi=300)

f, (ax) = plt.subplots(1, figsize=[6,2.5])
ax.plot(dts, ufxprof[:,:].sum(0), 'k', label='Total')
ax.plot(dts, ufxprof[0:140,:].sum(0), label='SP')
ax.plot(dts, ufxprof[580:719,:].sum(0), label='NP')
ax.plot(dts, ufxprof[141:579,:].sum(0), label='EQ')
#gw = scipy.signal.gaussian(10,std=3)
#ax.plot(crs, scipy.signal.convolve(ufxprof[:,:].sum(0), gw/gw.sum(), mode='same'), 'k')
#ax.plot(crs, scipy.signal.convolve(ufxprof[0:140,:].sum(0), gw/gw.sum(), mode='same'))
#ax.plot(crs, scipy.signal.convolve(ufxprof[580:719,:].sum(0), gw/gw.sum(), mode='same'))
#ax.plot(crs, scipy.signal.convolve(ufxprof[141:579,:].sum(0), gw/gw.sum(), mode='same'))
ax.legend(ncol=2, loc=2)
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Unsigned CH Magnetic Flux [Mx]')
plt.tight_layout()
plt.savefig('plt/flx.pdf')
plt.savefig('plt/flx.png', dpi=300)