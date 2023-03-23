"""
==================
CHORTLE
==================

Coronal Hole Observer and Regional Tracker for Long-term Examination

"""

# Import libraries

import configparser
import datetime
import os

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

import scipy.ndimage
import scipy.signal
import scipy.stats

import sunpy.coordinates
import sunpy.image
import sunpy.io.fits
import sunpy.map
import sunpy.visualization.colormaps

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from reproject import reproject_interp
from sunpy.net import Fido, attrs as a

import sunpy
import sunpy.io
import sunpy.coordinates
import datetime
from sunpy.net import Fido, attrs as a
import drms
import os
import glob


def chortle(cr, plot=False):

    # Read configuration file
    config = configparser.ConfigParser()
    config.read('config.cfg')

    # Specify directory structures from configuration file
    datdir = config['paths']['datdir']
    magdir = config['paths']['magdir']
    outdir = config['paths']['outdir']

    # Specify timing parameters
    nsday = 1

    # Grab the start and end dates for this rotation
    t0 = sunpy.coordinates.sun.carrington_rotation_time(cr)
    t0.format = 'datetime'
    t0 = t0.value
    t1 = sunpy.coordinates.sun.carrington_rotation_time(cr + 1)
    t1.format = 'datetime'
    t1 = t1.value

    # For now, round these to the nearest day
    if t0.hour > 12:
        t0 = t0 + datetime.timedelta(days=1)
    t0 = datetime.datetime(t0.year, t0.month, t0.day)
    if t1.hour > 12:
        t1 = t1 + datetime.timedelta(days=1)
    t1 = datetime.datetime(t1.year, t1.month, t1.day)

    # Download appropriate data
    search_aia = (a.Instrument('AIA') & a.Sample(nsday * u.day) & a.Time(t0, t1))
    res_aia = Fido.search(a.Wavelength(19.3 * u.nm), search_aia)

    search_sta = (a.Source('STEREO_A') & a.Instrument('SECCHI') & a.Detector('EUVI') & a.Sample(
        nsday * u.day) & a.Time(t0, t1))
    res_sta = Fido.search(a.Wavelength(19.5 * u.nm), search_sta)

    search_stb = (a.Source('STEREO_B') & a.Instrument('SECCHI') & a.Detector('EUVI') & a.Sample(
        nsday * u.day) & a.Time(t0, t1))
    res_stb = Fido.search(a.Wavelength(19.5 * u.nm), search_stb)

    # Running into some errors with automatically fetching HMI data
    # Manually download and load for now...
    # c = drms.Client()
    # c.pkeys('hmi_synoptic_mr_polfil_720s')
    # res_hmi = Fido.search(a.jsoc.Series('hmi_synoptic_mr_polfil_720s') & a.jsoc.PrimeKey('crnum', '2193'))

    files_aia = Fido.fetch(res_aia, path=datdir + 'aia/')
    files_aia.sort()

    files_sta = Fido.fetch(res_sta, path=datdir + 'sta/')
    files_sta.sort()

    files_stb = Fido.fetch(res_stb, path=datdir + 'stb/')
    files_stb.sort()

    skip_aia = skip_sta = skip_stb = False

    if len(files_aia) == 0:
        skip_aia = True
    if len(files_sta) == 0:
        skip_sta = True
    if len(files_stb) == 0:
        skip_stb = True

    # Generate some blank storage arrays
    oshape = [720, 1440]
    chmap = np.zeros(oshape, dtype=np.double)
    chmap[:, :] = np.nan

    # Grab a synoptic magnetogram
    br0 = sunpy.io.fits.read(magdir + 'hmi.synoptic_mr_polfil_720s.' + str(cr) + '.Mr_polfil.fits')[1].data
    br = sunpy.image.resample.resample(br0, oshape, method='linear')

    # Iterate through this data, reprojecting as you go
    for file in files_aia:

        # Read in data
        map_aia = sunpy.map.Map(file)

        # Check for bad files
        if map_aia.exposure_time == 0:
            continue

        # Construct an output header
        header = sunpy.map.make_fitswcs_header(np.empty(oshape),
                                               SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.AU,
                                                        frame="heliographic_carrington",
                                                        obstime=map_aia.date,
                                                        observer='self'),
                                               scale=[180 / oshape[0],
                                                      360 / oshape[1]] * u.deg / u.pix,
                                               projection_code="CAR")

        # header['crval2'] = 0
        out_wcs = WCS(header)

        # Reproject
        array, footprint = reproject_interp(
            (map_aia.data, map_aia.wcs), out_wcs, shape_out=oshape)
        omap_aia = sunpy.map.Map((array, header))
        omap_aia.plot_settings = map_aia.plot_settings

        # Normalize data to exposure time
        omap_aia_data = (omap_aia.data / (map_aia.exposure_time / u.second)).value

        # Condense the reprojected map minimums into chmap
        chmap = np.fmin(chmap, omap_aia_data)

    chmap_aia = np.copy(chmap)

    # Generate some blank storage arrays
    oshape = [720, 1440]
    chmap = np.zeros(oshape, dtype=np.double)
    chmap[:, :] = np.nan

    # Iterate through this data, reprojecting as you go
    for file in files_sta:

        # Read in data
        map_sta = sunpy.map.Map(file)

        # Check for bad files
        if map_sta.exposure_time == 0:
            continue

        # Construct an output header
        header = sunpy.map.make_fitswcs_header(np.empty(oshape),
                                               SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.AU,
                                                        frame="heliographic_carrington",
                                                        # frame="heliographic_stonyhurst",
                                                        obstime=map_sta.date,
                                                        observer='self'),
                                               scale=[180 / oshape[0],
                                                      360 / oshape[1]] * u.deg / u.pix,
                                               projection_code="CAR")

        # header['crval2'] = 0
        out_wcs = WCS(header)

        # Reproject
        array, footprint = reproject_interp(
            (map_sta.data, map_sta.wcs), out_wcs, shape_out=oshape)
        omap_sta = sunpy.map.Map((array, header))
        omap_sta.plot_settings = map_sta.plot_settings

        # Normalize data to exposure time
        omap_sta_data = (omap_sta.data / (map_sta.exposure_time / u.second)).value

        # Condense the reprojected map minimums into chmap
        chmap = np.fmin(chmap, omap_sta_data)

    chmap_sta = np.copy(chmap)

    # Generate some blank storage arrays
    oshape = [720, 1440]
    chmap = np.zeros(oshape, dtype=np.double)
    chmap[:, :] = np.nan

    # Iterate through this data, reprojecting as you go
    for file in files_stb:

        # Read in data
        map_stb = sunpy.map.Map(file)

        # Check for bad files
        if map_stb.exposure_time == 0:
            continue

        # Construct an output header
        header = sunpy.map.make_fitswcs_header(np.empty(oshape),
                                               SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.AU,
                                                        frame="heliographic_carrington",
                                                        # frame="heliographic_stonyhurst",
                                                        obstime=map_stb.date,
                                                        observer='self'),
                                               scale=[180 / oshape[0],
                                                      360 / oshape[1]] * u.deg / u.pix,
                                               projection_code="CAR")

        # header['crval2'] = 0
        out_wcs = WCS(header)

        # Reproject
        array, footprint = reproject_interp(
            (map_stb.data, map_stb.wcs), out_wcs, shape_out=oshape)
        omap_stb = sunpy.map.Map((array, header))
        omap_stb.plot_settings = map_stb.plot_settings

        # Normalize data to exposure time
        omap_stb_data = (omap_stb.data / (map_stb.exposure_time / u.second)).value

        # Condense the reprojected map minimums into chmap
        chmap = np.fmin(chmap, omap_stb_data)

    chmap_stb = np.copy(chmap)

    # Shift things and generate coordinates
    chmap_aia[:, 0] = (chmap_aia[:, -1] + chmap_aia[:, 1]) / 2
    chmap_aia[0, :] = np.nan
    chmap_aia = np.roll(chmap_aia, [0, int(oshape[1] / 2)])

    chmap_sta[:, 0] = (chmap_sta[:, -1] + chmap_sta[:, 1]) / 2
    chmap_sta[0, :] = np.nan
    chmap_sta = np.roll(chmap_sta, [0, int(oshape[1] / 2)])

    chmap_stb[:, 0] = (chmap_stb[:, -1] + chmap_stb[:, 1]) / 2
    chmap_stb[0, :] = np.nan
    chmap_stb = np.roll(chmap_stb, [0, int(oshape[1] / 2)])

    pscale = [oshape[0] / 180, oshape[1] / 360]

    # Generate threshold values for AIA
    qs = np.nanmedian(chmap_aia)
    if np.isfinite(qs):
        thrsh = np.array([])
        [dlat, dlon] = [60, 60]
        for ilat in np.arange(0, 180 / dlat):
            for ilon in np.arange(0, 360 / dlon):
                plat0 = int(ilat * dlat * pscale[0])
                plat1 = int((ilat + 1) * dlat * pscale[0])
                plon0 = int(ilon * dlon * pscale[1])
                plon1 = int((ilon + 1) * dlon * pscale[1])
                sarr = chmap_aia[plat0:plat1, plon0:plon1]

                sarr_hist = np.histogram(sarr[np.where(np.isfinite(sarr))].flatten(), bins=100, range=[0, qs])
                sh_x = sarr_hist[1][0:-1]
                sh_y = sarr_hist[0]
                sh_y2 = scipy.signal.convolve(sh_y, scipy.signal.hann(20), mode='same') / sum(scipy.signal.hann(20))
                pks = scipy.signal.find_peaks_cwt(sh_y2, np.arange(1, 20))
                if len(pks) >= 2:
                    sh_x2 = sh_x[pks[0]:pks[-1] - 1]
                    sh_y3 = sh_y2[pks[0]:pks[-1] - 1]
                else:
                    minval = int(len(sh_x) / 4)
                    maxval = int(len(sh_x) * 0.9)
                    sh_x2 = sh_x[minval:maxval]
                    sh_y3 = sh_y2[minval:maxval]
                pks2 = scipy.signal.find_peaks_cwt(-1 * sh_y3, np.arange(1, 20))
                if len(pks2) != 0:
                    thrsh = np.append(thrsh, sh_x2[pks2[0]])
                else:
                    thrsh = np.append(thrsh, np.nan)

        chval_aia = np.nanmean(thrsh)

    # Generate threshold values for STA
    qs = np.nanmedian(chmap_sta)
    if np.isfinite(qs):
        thrsh = np.array([])
        [dlat, dlon] = [60, 60]
        for ilat in np.arange(0, 180 / dlat):
            for ilon in np.arange(0, 360 / dlon):
                plat0 = int(ilat * dlat * pscale[0])
                plat1 = int((ilat + 1) * dlat * pscale[0])
                plon0 = int(ilon * dlon * pscale[1])
                plon1 = int((ilon + 1) * dlon * pscale[1])
                sarr = chmap_sta[plat0:plat1, plon0:plon1]

                sarr_hist = np.histogram(sarr[np.where(np.isfinite(sarr))].flatten(), bins=100, range=[0, qs])
                sh_x = sarr_hist[1][0:-1]
                sh_y = sarr_hist[0]
                sh_y2 = scipy.signal.convolve(sh_y, scipy.signal.hann(20), mode='same') / sum(scipy.signal.hann(20))
                pks = scipy.signal.find_peaks_cwt(sh_y2, np.arange(1, 20))
                if len(pks) >= 2:
                    sh_x2 = sh_x[pks[0]:pks[-1] - 1]
                    sh_y3 = sh_y2[pks[0]:pks[-1] - 1]
                else:
                    minval = int(len(sh_x) / 4)
                    maxval = int(len(sh_x) * 0.9)
                    sh_x2 = sh_x[minval:maxval]
                    sh_y3 = sh_y2[minval:maxval]
                pks2 = scipy.signal.find_peaks_cwt(-1 * sh_y3, np.arange(1, 20))
                if len(pks2) != 0:
                    thrsh = np.append(thrsh, sh_x2[pks2[0]])
                else:
                    thrsh = np.append(thrsh, np.nan)

        chval_sta = np.nanmean(thrsh)

    # Generate threshold values for STB
    qs = np.nanmedian(chmap_stb)
    if np.isfinite(qs):
        thrsh = np.array([])
        [dlat, dlon] = [60, 60]
        for ilat in np.arange(0, 180 / dlat):
            for ilon in np.arange(0, 360 / dlon):
                plat0 = int(ilat * dlat * pscale[0])
                plat1 = int((ilat + 1) * dlat * pscale[0])
                plon0 = int(ilon * dlon * pscale[1])
                plon1 = int((ilon + 1) * dlon * pscale[1])
                sarr = chmap_stb[plat0:plat1, plon0:plon1]

                sarr_hist = np.histogram(sarr[np.where(np.isfinite(sarr))].flatten(), bins=100, range=[0, qs])
                sh_x = sarr_hist[1][0:-1]
                sh_y = sarr_hist[0]
                sh_y2 = scipy.signal.convolve(sh_y, scipy.signal.hann(20), mode='same') / sum(scipy.signal.hann(20))
                pks = scipy.signal.find_peaks_cwt(sh_y2, np.arange(1, 20))
                if len(pks) >= 2:
                    sh_x2 = sh_x[pks[0]:pks[-1] - 1]
                    sh_y3 = sh_y2[pks[0]:pks[-1] - 1]
                else:
                    minval = int(len(sh_x) / 4)
                    maxval = int(len(sh_x) * 0.9)
                    sh_x2 = sh_x[minval:maxval]
                    sh_y3 = sh_y2[minval:maxval]
                pks2 = scipy.signal.find_peaks_cwt(-1 * sh_y3, np.arange(1, 20))
                if len(pks2) != 0:
                    thrsh = np.append(thrsh, sh_x2[pks2[0]])
                else:
                    thrsh = np.append(thrsh, np.nan)

        chval_stb = np.nanmean(thrsh)

    chmap0_aia = np.copy(chmap_aia)
    if not skip_aia:
        chmap0_aia[np.where(chmap_aia > chval_aia)] = 0

    chmap0_sta = np.copy(chmap_sta)
    if not skip_sta:
        chmap0_sta[np.where(chmap_sta > chval_sta)] = 0

    chmap0_stb = np.copy(chmap_stb)
    if not skip_stb:
        chmap0_stb[np.where(chmap_stb > chval_stb)] = 0

    # Generate a merged chmap
    chmap0 = (np.nansum(np.dstack((chmap0_aia, chmap0_sta, chmap0_stb)), 2) != 0)

    ocstruct = np.ones([3, 3])
    chmap1 = scipy.ndimage.binary_opening(scipy.ndimage.binary_closing(chmap0, structure=ocstruct, iterations=1),
                                          structure=ocstruct, iterations=1) * chmap0

    # Label each dark region for comparison
    labstruct = np.ones([3, 3])
    lchmap = (scipy.ndimage.label(chmap1, structure=labstruct))[0]

    for i in np.arange(2, lchmap.max() + 1):
        chmag = br[np.where(lchmap == i)]
        if len(chmag) < 10:
            chmap1[np.where(lchmap == i)] = 0
            continue
        if abs(scipy.stats.skew(chmag)) < 0.5:
            chmap1[np.where(lchmap == i)] = 0
            continue

    chmap = chmap1

    # For visualization, scale the images to create a merged view
    im_aia = np.copy(chmap_aia)
    im_sta = np.copy(chmap_sta)
    im_stb = np.copy(chmap_stb)

    im_aia = im_aia - np.nanmin(im_aia)
    im_aia = im_aia / np.nanmax(im_aia)

    im_sta = im_sta - np.nanmin(im_sta)
    im_sta = im_sta / np.nanmax(im_sta)

    im_stb = im_stb - np.nanmin(im_stb)
    im_stb = im_stb / np.nanmax(im_stb)

    chim = np.fmin(im_sta, im_aia)
    chim = np.fmin(chim, im_stb)

    # Generate latitude and longitude coordinates
    lats = np.linspace(-90,90,oshape[0])
    lons = np.linspace(0,360,oshape[1])

    # Save everything out to file
    fname = outdir + 'chmap/chmap-' + str(cr) + '.fits'
    sunpy.io.write_file(fname, chmap * chim, header, overwrite=True)

    fname = outdir + 'chmap/chmap-' + str(cr) + '-chim.fits'
    sunpy.io.write_file(fname, chim, header, overwrite=True)

    if plot:

        # Some plotting
        f, (ax) = plt.subplots(1, figsize=[6,3])
        ax.imshow(chim, extent=[0,360,-90,90], cmap=sunpy.visualization.colormaps.cm.sdoaia193, vmin=0, vmax=0.25)
        ax.contour(lons, lats, chmap, colors='teal',linewidths=0.5)
        ax.set_xlabel('Longitude (degrees)')
        ax.set_ylabel('Latitude (degrees)')
        ax.set_title('CH - AIA/EUVI - CR '+str(cr))
        plt.tight_layout()
        plt.savefig(outdir+'plt/plt-chmap-'+str(cr)+'.pdf')
        plt.savefig(outdir+'plt/plt-chmap-'+str(cr)+'.png', dpi=150)
        plt.close('all')


def grab_data():

    # Specify any directories
    hmidat = os.path.expanduser('~/data/hmi.Synoptic_Mr.polfil/')

    # Sort out the last downloaded rotation
    crfiles = glob.glob(hmidat+'*.fits')
    crfiles.sort()
    crlist = [int(i[-19:-15]) for i in crfiles]

    # Specify requested rotations
    if len(crlist) == 0:
        cr0 = 2096
    else:
        cr0 = max(crlist) + 1
    cr1 = int(sunpy.coordinates.sun.carrington_rotation_number(t='now'))

    if (cr0 - 1) == cr1:
        print('what?')

    # Start the client
    c = drms.Client()

    # Generate a search
    crots = a.jsoc.PrimeKey('CAR_ROT', str(cr0) + '-' + str(cr1))
    res = Fido.search(a.jsoc.Series('hmi.Synoptic_Mr_polfil_720s'), crots, a.jsoc.Notify(os.environ["JSOC_EMAIL"]))

    # Once the query is made and trimmed down...
    download = Fido.fetch(res, path=hmidat+'{file}.fits')


def genprof(cr0, cr1):
    # Read configuration file
    config = configparser.ConfigParser()
    config.read('config.cfg')

    # Specify directory structures from configuration file
    datdir = config['paths']['datdir']
    magdir = config['paths']['magdir']
    outdir = config['paths']['outdir']

    crs = np.arange(cr0, cr1 + 1)

    oshape = [720, 1440]
    ncrs = len(crs)

    chprof0 = np.zeros([oshape[0], ncrs], dtype=np.double)
    chprof0[:, :] = np.nan

    chprof = np.zeros([oshape[0], ncrs], dtype=np.double)
    chprof[:, :] = np.nan

    improf = np.zeros([oshape[0], ncrs], dtype=np.double)
    improf[:, :] = np.nan

    sfxprof = np.zeros([oshape[0], ncrs], dtype=np.double)
    sfxprof[:, :] = np.nan

    ufxprof = np.zeros([oshape[0], ncrs], dtype=np.double)
    ufxprof[:, :] = np.nan

    sfxprof90 = np.zeros([oshape[0], ncrs], dtype=np.double)
    sfxprof90[:, :] = np.nan

    ufxprof90 = np.zeros([oshape[0], ncrs], dtype=np.double)
    ufxprof90[:, :] = np.nan

    sfxprof80 = np.zeros([oshape[0], ncrs], dtype=np.double)
    sfxprof80[:, :] = np.nan

    ufxprof80 = np.zeros([oshape[0], ncrs], dtype=np.double)
    ufxprof80[:, :] = np.nan

    sfxprof70 = np.zeros([oshape[0], ncrs], dtype=np.double)
    sfxprof70[:, :] = np.nan

    ufxprof70 = np.zeros([oshape[0], ncrs], dtype=np.double)
    ufxprof70[:, :] = np.nan

    for cr in crs:

        fname = outdir + 'chmap/chmap-' + str(cr) + '.fits'
        if not os.path.exists(fname):
            continue

        [chdat, chhdr] = (sunpy.io.read_file(fname))[0]

        fname = outdir + 'chmap/chmap-' + str(cr) + '-chim.fits'
        [imdat, imhdr] = (sunpy.io.read_file(fname))[0]

        br0 = sunpy.io.fits.read(magdir + 'hmi.synoptic_mr_polfil_720s.' + str(cr) + '.Mr_polfil.fits')[1].data
        br = sunpy.image.resample.resample(br0, oshape, method='linear')
        pscale = 4 * np.pi * 6.957e10 ** 2 / (oshape[0] * oshape[1])

        chtop = np.nanmax(chdat)
        chbot = np.nanmin(chdat[np.where(chdat != 0)])
        chdep = chtop - chbot

        ch90 = np.logical_and(chdat >= chbot, chdat <= (0.9 * chtop))
        ch80 = np.logical_and(chdat >= chbot, chdat <= (0.8 * chtop))
        ch70 = np.logical_and(chdat >= chbot, chdat <= (0.7 * chtop))

        chprof0[:, cr - cr0] = chdat.sum(1) / oshape[1]
        chprof[:, cr - cr0] = (chdat != 0).sum(1, where=np.isfinite(chdat)) / oshape[1]
        improf[:, cr - cr0] = imdat.sum(1) / oshape[1]
        sfxprof[:, cr - cr0] = (br * pscale * np.logical_and(chdat != 0, np.isfinite(chdat))).sum(1)
        ufxprof[:, cr - cr0] = (np.abs(br * pscale * np.logical_and(chdat != 0, np.isfinite(chdat)))).sum(1)
        sfxprof90[:, cr - cr0] = (br * pscale * ch90).sum(1)
        ufxprof90[:, cr - cr0] = (np.abs(br * pscale * ch90)).sum(1)
        sfxprof80[:, cr - cr0] = (br * pscale * ch80).sum(1)
        ufxprof80[:, cr - cr0] = (np.abs(br * pscale * ch80)).sum(1)
        sfxprof70[:, cr - cr0] = (br * pscale * ch70).sum(1)
        ufxprof70[:, cr - cr0] = (np.abs(br * pscale * ch70)).sum(1)

    np.save(outdir + 'dat/chprof0.npy', chprof0)
    np.save(outdir + 'dat/chprof.npy', chprof)
    np.save(outdir + 'dat/improf.npy', improf)
    np.save(outdir + 'dat/sfxprof.npy', sfxprof)
    np.save(outdir + 'dat/ufxprof.npy', ufxprof)
    np.save(outdir + 'dat/sfxprof90.npy', sfxprof90)
    np.save(outdir + 'dat/ufxprof90.npy', ufxprof90)
    np.save(outdir + 'dat/sfxprof80.npy', sfxprof80)
    np.save(outdir + 'dat/ufxprof80.npy', ufxprof80)
    np.save(outdir + 'dat/sfxprof70.npy', sfxprof70)
    np.save(outdir + 'dat/ufxprof70.npy', ufxprof70)
    np.save(outdir + 'dat/crs.npy', crs)

def chortle_eit(cr, plot=False):

    # Read configuration file
    config = configparser.ConfigParser()
    config.read('config.cfg')

    # Specify directory structures from configuration file
    datdir = config['paths']['datdir']
    magdir = config['paths']['magdir']
    outdir = config['paths']['outdir']

    # Specify timing parameters
    nshour = 6

    # Grab the start and end dates for this rotation
    t0 = sunpy.coordinates.sun.carrington_rotation_time(cr)
    t0.format = 'datetime'
    t0 = t0.value
    t1 = sunpy.coordinates.sun.carrington_rotation_time(cr+1)
    t1.format = 'datetime'
    t1 = t1.value

    # For now, round these to the nearest day
    if t0.hour > 12:
        t0 = t0 + datetime.timedelta(days=1)
    t0 = datetime.datetime(t0.year, t0.month, t0.day)
    if t1.hour > 12:
        t1 = t1 + datetime.timedelta(days=1)
    t1 = datetime.datetime(t1.year, t1.month, t1.day)

    # Download appropriate data
    search_eit = (a.Instrument('EIT') & a.vso.Sample(nshour*u.hour) & a.Time(t0,t1))
    res_eit = Fido.search(a.Wavelength(19.5 * u.nm), search_eit)

    files_eit = Fido.fetch(res_eit, path=datdir+'eit/')
    files_eit.sort()

    ## Grab a synoptic magnetogram
    br = sunpy.io.fits.read(magdir + 'synop_Mr_0.polfil.'+str(cr)+'.fits')[0].data

    ## Generate some blank storage arrays
    oshape = [720,1440]
    chmap = np.zeros(oshape,dtype=np.double)
    chmap[:,:] = np.nan

    ## Iterate through this data, reprojecting as you go
    for file in files_eit:

        # Read in data
        map_eit = sunpy.map.Map(file)

        # Check for bad files
        if map_eit.exposure_time == 0: continue

        # Remove any missing data
        temp_header = map_eit.meta
        map_eit = sunpy.map.Map((map_eit.data.astype(np.float), temp_header))
        map_eit.data[np.where(map_eit.data == 0)] = np.nan
        
        # Construct an output header
        header = sunpy.map.make_fitswcs_header(np.empty(oshape),
                       SkyCoord(0, 0, unit=u.deg,
                                frame="heliographic_carrington",
                                obstime=map_eit.date),
                       scale=[180 / oshape[0],
                              360 / oshape[1]] * u.deg / u.pix,
                       projection_code="CAR")

        out_wcs = WCS(header)

        # Reproject
        array, footprint = reproject_interp(
            (map_eit.data, map_eit.wcs), out_wcs, shape_out=oshape)
        omap_eit = sunpy.map.Map((array, header))
        omap_eit.plot_settings = map_eit.plot_settings

        # Normalize data to exposure time
        omap_eit_data = (omap_eit.data / (map_eit.exposure_time / u.second)).value

        # Condense the reprojected map minimums into chmap
        chmap = np.fmin(chmap, omap_eit_data)

    chmap_eit = np.copy(chmap)

    # Shift things and generate coordinates
    chmap_eit[:,0] = (chmap_eit[:,-1]+chmap_eit[:,1])/2
    chmap_eit[0,:] = np.nan
    chmap_eit = np.roll(chmap_eit, [0,int(oshape[1]/2)])

    lats = np.linspace(-90,90,oshape[0])
    lons = np.linspace(0,360,oshape[1])
    pscale = [oshape[0]/180, oshape[1]/360]

    # Generate threshold values for AIA
    qs = np.median(chmap_eit[np.isfinite(chmap_eit)])
    thrsh = np.array([])
    [dlat, dlon] = [60,60]
    for ilat in np.arange(0,180/dlat):
        for ilon in np.arange(0,360/dlon):
            plat0 = int(ilat*dlat*pscale[0])
            plat1 = int((ilat+1)*dlat*pscale[0])
            plon0 = int(ilon*dlon*pscale[1])
            plon1 = int((ilon+1)*dlon*pscale[1])
            sarr = chmap_eit[plat0:plat1, plon0:plon1]

            sarr_hist = np.histogram(sarr[np.where(np.isfinite(sarr))].flatten(), bins=100, range=[np.nanmin(sarr),np.nanmax(sarr)])
            sh_x = sarr_hist[1][0:-1]
            sh_y = sarr_hist[0]
            sh_y2 = scipy.signal.convolve(sh_y, scipy.signal.hann(20), mode='same')/sum(scipy.signal.hann(20))
            pks = scipy.signal.find_peaks_cwt(sh_y2,np.arange(1,20))
            if len(pks) >= 2:
                sh_x2 = sh_x[pks[0]:pks[-1]-1]
                sh_y3 = sh_y2[pks[0]:pks[-1]-1]
            else:
                minval = int(len(sh_x)/4)
                maxval = int(len(sh_x)*0.9)
                sh_x2 = sh_x[minval:maxval]
                sh_y3 = sh_y2[minval:maxval]
            pks2 = scipy.signal.find_peaks_cwt(-1*sh_y3,np.arange(1,20))
            if len(pks2) != 0:
                thrsh = np.append(thrsh, sh_x2[pks2[0]])
            else:
                thrsh = np.append(thrsh, np.nan)

    chval_eit = np.nanmean(thrsh)

    chmap0_eit = np.copy(chmap_eit)
    chmap0_eit[np.where(chmap_eit>chval_eit)] = 0

    # Generate a merged chmap
    chmap0 = chmap_eit

    ocstruct = np.ones([3,3])
    chmap1 = scipy.ndimage.binary_opening(scipy.ndimage.binary_closing(chmap0, structure=ocstruct, iterations=1), structure=ocstruct, iterations=1) * chmap0

    # Label each dark region for comparison
    labstruct = np.ones([3,3])
    lchmap = (scipy.ndimage.label(chmap1, structure=labstruct))[0]

    for i in np.arange(2, lchmap.max()+1):
        chmag = br[np.where(lchmap==i)]
        if len(chmag) < 10:
            chmap1[np.where(lchmap==i)] = 0
            continue
        if abs(scipy.stats.skew(chmag)) < 0.5:
            chmap1[np.where(lchmap==i)] = 0
            continue

    chmap = chmap1

    # For visualization, scale the images to create a merged view
    im_eit = np.copy(chmap_eit)

    im_eit = im_eit - np.nanmin(im_eit)
    im_eit = im_eit / np.nanmax(im_eit)

    chim = im_eit

    # Save everything out to file
    fname = outdir+'chmap/chmap-'+str(cr)+'-eit.fits'
    sunpy.io.write_file(fname, chmap*chim, header, overwrite=True)

    if plot:

        # Some plotting
        f, (ax) = plt.subplots(1, figsize=[6,3])
        ax.imshow(chim, extent=[0,360,-90,90], cmap=sunpy.visualization.colormaps.cm.sohoeit195, vmin=0, vmax=0.5)
        ax.contour(lons, lats, chmap, colors='teal',linewidths=0.5)
        ax.set_xlabel('Longitude (degrees)')
        ax.set_ylabel('Latitude (degrees)')
        ax.set_title('CH - EIT - CR '+str(cr))
        plt.tight_layout()
        plt.savefig(outdir+'chmap/plt-chmap-'+str(cr)+'-eit.pdf')
        plt.savefig(outdir+'chmap/plt-chmap-'+str(cr)+'-eit.png', dpi=150)