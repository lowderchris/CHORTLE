import matplotlib.pyplot as plt
import numpy as np

import scipy.stats
import scipy.signal
import scipy.ndimage

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import datetime

from reproject import reproject_interp

import sunpy.map
import sunpy.io
import sunpy.visualization.colormaps

import glob
import os

from sunpy.net import Fido, attrs as a
import drms
import sunpy.visualization.colormaps

import configparser

# def genprof(cr0, cr1):

# Read configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

# Specify directory structures from configuration file
datdir = config['paths']['datdir']
magdir = config['paths']['magdir']
outdir = config['paths']['outdir']

cr0 = 2097
cr1 = 2232
crs = np.arange(cr0,cr1+1)

oshape = [720,1440]
ncrs = len(crs)

chprof0 = np.zeros([oshape[0], ncrs], dtype=np.double)
chprof0[:,:] = np.nan
chprof = np.zeros([oshape[0], ncrs], dtype=np.double)
chprof[:,:] = np.nan
improf = np.zeros([oshape[0], ncrs], dtype=np.double)
improf[:,:] = np.nan

for cr in crs:

    fname = outdir+'chmap/chmap-'+str(cr)+'.fits'
    if not os.path.exists(fname): continue

    [chdat, chhdr] = (sunpy.io.read_file(fname))[0]

    fname = outdir+'chmap/chmap-'+str(cr)+'-chim.fits'
    [imdat, imhdr] = (sunpy.io.read_file(fname))[0]

    chprof0[:,cr-cr0] = chdat.sum(1)/oshape[1]
    chprof[:,cr-cr0] = (chdat!=0).sum(1)/oshape[1]
    improf[:,cr-cr0] = imdat.sum(1)/oshape[1]

np.save('chprof0.npy',chprof0)
np.save('chprof.npy',chprof)
np.save('improf.npy',improf)
np.save('crs.npy', crs)

chprof0 = np.load('chprof0.npy')
chprof = np.load('chprof.npy')
improf = np.load('improf.npy')
crs = np.load('crs.npy')
