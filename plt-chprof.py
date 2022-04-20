import sunpy.visualization.colormaps
import sunpy.coordinates
import scipy.signal
import configparser

import sunpy.io

import numpy as np

import matplotlib.dates as mdates

light_figs()
#keynote_figs()

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

f, (ax) = subplots(1, figsize=[6,2.5])
im = ax.imshow(chprof0, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', vmin=0, vmax=0.02)
cb = colorbar(im, ax=ax, label='Fractional CH Coverage', extend='max')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
tight_layout()
savefig('plt/chprof0.pdf')
savefig('plt/chprof0.png', dpi=300)

f, (ax) = subplots(1, figsize=[6,2.5])
im = ax.imshow(chprof, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto')
cb = colorbar(im, ax=ax, label='CH Coverage')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
tight_layout()
savefig('plt/chprof.pdf')
savefig('plt/chprof.png', dpi=300)

f, (ax) = subplots(1, figsize=[6,2.5])
im = ax.imshow(improf, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap=sunpy.visualization.colormaps.cm.sdoaia193)
cb = colorbar(im, ax=ax, label='Input EUV Data')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
tight_layout()
savefig('plt/improf.pdf')
savefig('plt/improf.png', dpi=300)

f, (ax) = subplots(1, figsize=[6,2.5])
im = ax.imshow(sfxprof, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap='RdBu_r', vmin=-1e20, vmax=1e20)
cb = colorbar(im, ax=ax, label='CH Signed Flux [Mx]', extend='both')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
tight_layout()
savefig('plt/sfxprof.pdf')
savefig('plt/sfxprof.png', dpi=300)

f, (ax) = subplots(1, figsize=[6,2.5])
im = ax.imshow(ufxprof, extent=[ndts[0],ndts[-1], -90, 90], aspect='auto', cmap='magma')
cb = colorbar(im, ax=ax, label='CH Unsigned Flux [Mx]')
axb = ax.twiny()
axb.set_xlim(crs[0], crs[-1])
ax.set_xlabel('Date')
axb.set_xlabel('Carrington Rotation')
ax.set_ylabel('Latitude')
ax.xaxis_date()
tight_layout()
savefig('plt/ufxprof.pdf')
savefig('plt/ufxprof.png', dpi=300)

f, (ax) = subplots(1, figsize=[6,2.5])
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
tight_layout()
savefig('plt/flx.pdf')
savefig('plt/flx.png', dpi=300)

# Load and plot one CR
#config = configparser.ConfigParser()
#config.read('config.cfg')

# Specify directory structures from configuration file
#datdir = config['paths']['datdir']
#magdir = config['paths']['magdir']
#outdir = config['paths']['outdir']

#cr = 2190

#fname = outdir + 'chmap/chmap-' + str(cr) + '.fits'
#[imdat, imhdr] = (sunpy.io.read_file(fname))[0]

# Persistence map plotting
