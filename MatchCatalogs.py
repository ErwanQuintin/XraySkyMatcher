import subprocess
import shlex
import matplotlib
matplotlib.use('TkAgg')
from LoadMasterSources import *
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
from matplotlib import rc
plt.rcParams.update({'font.size': 15})
rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})


#TODO: change the name_of_pointing for each different pointing you do, and its position (in degrees)
name_of_pointing = 'J0437-4715'
pointing = (69.251, -47.2501000) #Must be in degrees
pointing_coords = SkyCoord(pointing[0] * u.degree, pointing[1] * u.degree, frame="icrs")


#Workflow: Once you input ra, dec, and radius search, it will do a cone search in the multi-instrument catalog (to get
#source names) in the error region. It will also cut the catalogs to only keep relevant sources.
#Then, it loads master_source objects (which are custom objects for multi-instrument sources). These have all the
#information you need.

#NOTE: There are also upper limits. For the NICER pointing application, I removed them. If they are needed, you can
#ask me

path_to_file = os.getcwd()
catalog_datapath=os.path.join(path_to_file, 'FullCatalogs')
outputname= os.path.join(path_to_file, 'Pointings', name_of_pointing)
os.mkdir(outputname)

catalogs = ["XMM","Chandra","Swift","eRosita","Slew","RASS","WGACAT","Stacked"]

def select_master_sources_around_region(ra, dec, radius, outputname):
    """Radius is in arcminutes"""
    print(f"Extracting sources around region: RA {ra} and Dec {dec}")
    cmd= (f"stilts tpipe {os.path.join(catalog_datapath,'Master_source.fits')} cmd='"+
          f'select skyDistanceDegrees({ra},{dec},MS_RA,MS_DEC)*60<{radius} '+
          f"' out={os.path.join(outputname, 'Master_source_cone.fits')}")
    cmd=shlex.split(cmd)
    subprocess.run(cmd)

def select_catalogsources_around_region(outputname):
    print('Selecting catalog sources')
    for cat in catalogs:
        path_to_cat = os.path.join(catalog_datapath, cat)
        cmd = (f"stilts tmatch2 matcher=exact in1='{os.path.join(outputname, 'Master_source_cone.fits')}'\
         in2='{os.path.join(catalog_datapath,cat)}.fits' out='{os.path.join(outputname,cat)}.fits'\
         values1='{cat}' values2='{cat}_IAUNAME' find=all progress=none")
        cmd = shlex.split(cmd)
        subprocess.run(cmd)

select_master_sources_around_region(pointing[0], pointing[1], 60, outputname)
select_catalogsources_around_region(outputname)
master_sources = load_master_sources(outputname)

for multi_instrument_source in list(master_sources.values())[:5]:
    #Each multi_instrument_source is an object with the underlying catalog sources associated with it

    #Here we compute the off-axis angle between the source and the pointing
    source_coords = SkyCoord(multi_instrument_source.ra*u.degree, multi_instrument_source.dec*u.degree, frame="icrs")
    off_axis = pointing_coords.separation(source_coords)

    plt.figure()
    for catalog in multi_instrument_source.sources.keys():
        #If a given catalog is contained in this source, it will be in the "sources" dictionary, catalog as key,
        #source object as value
        catalog_source = multi_instrument_source.sources[catalog]
        tab_width = 2 * np.array(band_half_width[catalog])
        for band_det in range(len(catalog_source.band_flux)):
            #The band fluxes are stored in catalog_source.band_flux. They're in erg/s/cm2, so divide by tab_width to
            #be in erg/s/cm2/keV. Here I plot them, but you can do whatever you want with those
            plt.step(band_edges[catalog],
                     [catalog_source.band_flux[band_det][0] / tab_width[0]]
                     + list(catalog_source.band_flux[band_det] / tab_width),
                     c=colors[catalog], where='pre')
            plt.errorbar(band_center[catalog], catalog_source.band_flux[band_det] / tab_width,
                         yerr=[catalog_source.band_fluxerr[0][band_det] / tab_width,
                               catalog_source.band_fluxerr[1][band_det] / tab_width],
                         fmt="o", markeredgecolor='gray', c=colors[catalog], alpha=0.4)
        plt.step([],[],c=colors[catalog],label=catalog_source.name)
    plt.xlabel("Energy (keV)")
    plt.ylabel(r"$F_{\nu}$ ($\mathrm{erg.s}^{-1}.\mathrm{cm}^{-2}.\mathrm{keV}^{-1}$)")
    plt.legend()
    plt.loglog()
    plt.show()