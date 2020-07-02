"""
==================
CHORTLE
==================

Coronal Hole Observer and Regional Tracker for Long-term Examination
Quick bit of code to work with EIT data - eventually merged into one grand CHORTLE code that will automatically choose an appropriate dataset

"""

def chortle(cr):

    # Import libraries

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

    from sunpy.net import Fido, attrs as a
    import drms
    import sunpy.visualization.colormaps

    # Specify directory structures
    datdir = '/Users/clowder/data/chortle/'
    outdir = '/Users/clowder/data/chortle/'

    # Specify timing parameters
    #cr = 2193
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
    br = sunpy.io.fits.read('/Users/clowder/data/mdi.Synoptic_Mr.polfil/synop_Mr_0.polfil.'+str(cr)+'.fits')[0].data

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
        # This is a weird workaround since I can't add nans to int16 arrs
        temp_header = map_eit.meta
        map_eit = sunpy.map.Map((map_eit.data.astype(np.float), temp_header))
        map_eit.data[where(map_eit.data == 0)] = np.nan
        
        # Construct an output header
        header = sunpy.map.make_fitswcs_header(np.empty(oshape),
                       SkyCoord(0, 0, unit=u.deg,
                                frame="heliographic_carrington",
                                obstime=map_eit.date),
                       scale=[180 / oshape[0],
                              360 / oshape[1]] * u.deg / u.pix,
                       projection_code="CAR")

        #header['crval2'] = 0
        out_wcs = WCS(header)

        # Reproject
        array, footprint = reproject_interp(
            (map_eit.data, map_eit.wcs), out_wcs, shape_out=oshape)
        omap_eit = sunpy.map.Map((array, header))
        omap_eit.plot_settings = map_eit.plot_settings

        # Normalize data to exposure time
        omap_eit_data = (omap_eit.data / (map_eit.exposure_time / u.second)).value

        # Condense the reprojected map minimums into chmap
        chmap = numpy.fmin(chmap, omap_eit_data)

    chmap_eit = np.copy(chmap)

    # Shift things and generate coordinates
    # CL - Temporary fix for weird start NaN column and row
    # CL - Might also need to shift coordinates to sine latitude...
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

            sarr_hist = histogram(sarr[where(np.isfinite(sarr))].flatten(), bins=100, range=[np.nanmin(sarr),np.nanmax(sarr)])
            #sarr_dist = scipy.stats.rv_histogram(sarr_hist)
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
    chmap0_eit[where(chmap_eit>chval_eit)] = 0

    # Generate a merged chmap
    # CL - make changes here to restore to original behavior of measuring CH depth. Might need to create a normalized merged of AIA / STA data to fill the gaps
    chmap0 = chmap_eit

    ocstruct = np.ones([3,3])
    chmap1 = scipy.ndimage.binary_opening(scipy.ndimage.binary_closing(chmap0, structure=ocstruct, iterations=1), structure=ocstruct, iterations=1) * chmap0

    # Label each dark region for comparison
    labstruct = np.ones([3,3])
    lchmap = (scipy.ndimage.label(chmap1, structure=labstruct))[0]

    for i in np.arange(2, lchmap.max()+1):
        chmag = br[where(lchmap==i)]
        if len(chmag) < 10:
            chmap1[where(lchmap==i)] = 0
            continue
        if abs(scipy.stats.skew(chmag)) < 0.5:
            chmap1[where(lchmap==i)] = 0
            continue

    chmap = chmap1

    # For visualization, scale the images to created a merged view
    im_eit = np.copy(chmap_eit)

    im_eit = im_eit - np.nanmin(im_eit)
    im_eit = im_eit / np.nanmax(im_eit)

    chim = im_eit

    # Save everything out to file
    fname = outdir+'chmap/chmap-'+str(cr)+'-eit.fits'
    sunpy.io.write_file(fname, chmap*chim, header, overwrite=True)

    # Some plotting
    f, (ax) = subplots(1, figsize=[6,3])
    ax.imshow(chim, extent=[0,360,-90,90], cmap=sunpy.visualization.colormaps.cm.sohoeit195, vmin=0, vmax=0.5)
    ax.contour(lons, lats, chmap, colors='teal',linewidths=0.5)
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title('CH - EIT - CR '+str(cr))
    tight_layout()
    savefig(outdir+'chmap/plt-chmap-'+str(cr)+'-eit.pdf')
    savefig(outdir+'chmap/plt-chmap-'+str(cr)+'-eit.png', dpi=150)
