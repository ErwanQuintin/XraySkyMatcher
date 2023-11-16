import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm
import os
import webbrowser
import seaborn as sns
from itertools import combinations
from matplotlib import rc
plt.rcParams.update({'font.size': 15})
rc('text', usetex=True)


style="bmh"
cmap_to_use="turbo"

def click_action(ra, dec, xmm_name, swift_name, has_efeds=False, has_sdss=False):
    url_esaskyDSS = "http://sky.esa.int/?target="+str(np.round(ra,4))+" "+str(np.round(dec,4))+"&hips=DSS2+color&fov=0.1&cooframe=J2000&sci=true&lang=en"
    url_esaskyXMM = "http://sky.esa.int/?target=" + str(np.round(ra, 4)) + " " + str(
        np.round(dec, 4)) + "&hips=XMM-Newton+EPIC+color&fov=0.1&cooframe=J2000&sci=true&lang=en"
    url_esaskyChandra = "http://sky.esa.int/?target=" + str(np.round(ra, 4)) + " " + str(
        np.round(dec, 4)) + "&hips=Chandra+RGB&fov=0.1&cooframe=J2000&sci=true&lang=en"
    webbrowser.get('firefox').open(url_esaskyDSS)
    webbrowser.get('firefox').open(url_esaskyXMM, new=0)
    webbrowser.get('firefox').open(url_esaskyChandra, new=0)
    if has_efeds:
        url_esaskyEFEDS = "http://sky.esa.int/?target=" + str(np.round(ra, 4)) + " " + str(
            np.round(dec, 4)) + "&hips=eFEDS+RGB&fov=0.1&cooframe=J2000&sci=true&lang=en"
        webbrowser.get('firefox').open(url_esaskyEFEDS)
    if xmm_name != []:
        xmm_name = xmm_name[5:].replace(' ', '%2B')
        url_xmmssc = f"http://xmm-catalog.irap.omp.eu/sources?f={xmm_name}"
        webbrowser.get('firefox').open(url_xmmssc, new=0)
    if swift_name != []:
        url_swift = "https://www.swift.ac.uk/2SXPS/" + swift_name
        webbrowser.get('firefox').open(url_swift, new=0)
    if has_sdss:
        url_sdss = "https://skyserver.sdss.org/dr16/en/tools/explore/summary.aspx?ra="+str(ra)+"&dec="+str(dec)
        webbrowser.get('firefox').open(url_sdss, new=0)
    url_simbad = "http://simbad.u-strasbg.fr/simbad/sim-coo?Coord="+str(ra)+"+"+str(dec)+"&Radius=1&Radius.unit=arcmin&submit=submit+query"
    webbrowser.get('firefox').open(url_simbad)
    url_rosat= f"http://xmm-ssc.irap.omp.eu/claxson/xray_analyzer2.php?srcquery={ra}%20{dec}"
    webbrowser.get('firefox').open(url_rosat)

catalogs = ["XMM","Chandra","Swift","eRosita","Slew","RASS","WGACAT","Stacked"]#,"NewXMM"]

posErr_Names = {}
src_names={}
colors = {}
for ind,cat in enumerate(catalogs):
    posErr_Names[cat]=f"{cat}_PosErr"
    src_names[cat] = f"{cat}_IAUNAME"
    colors[cat] = matplotlib.cm.get_cmap(cmap_to_use)(ind / len(catalogs))

#Defining all the catalog-related column names
flux_names={"XMM": "EP_8_FLUX",
            "Chandra":"flux_powlaw_aper_b",
            "Swift":"Flux",
            "eRosita":"ML_FLUX",
            "Slew":"Flux",
            "Stacked":"EP_FLUX",
            "RASS":"Flux",
            "WGACAT":"Flux"}
flux_error_names={"XMM": ["EP_8_FLUX_ERR","EP_8_FLUX_ERR"],
                  "Chandra":["flux_powlaw_aper_b_negerr","flux_powlaw_aper_b_poserr"],
                  "Swift":["FluxErr_neg","FluxErr_pos"],
                  "eRosita":["ML_FLUX_ERR","ML_FLUX_ERR"],
                  "Slew":["FluxErr","FluxErr"],
                  "Stacked":["EP_FLUX_ERR","EP_FLUX_ERR"],
                  "RASS":["FluxErr","FluxErr"],
                  "WGACAT":["FluxErr","FluxErr"]}
conv_factors = {"XMM": 1/0.999,
                "NewXMM": 1/0.999,
                "Chandra":1/0.69,
                "Swift":1/0.9,
                "eRosita":1/0.39,
                "Slew":1/0.999,
                "Stacked":1/0.999,
                "RASS":1/0.35,
                "WGACAT":1/0.35}
time_names={"XMM": "MJD_START",
            "Chandra":"gti_mjd_obs",
            "Swift":"MidTime_MJD",
            "eRosita":"MJD_OBS",
            "Slew":"DATE_OBS",
            "Stacked":"MJD_FIRST",
            "RASS":"OBS_DATE_1",
            "WGACAT":"StartDate",
            "OM":"MJD_START",
            "UVOT":"DATE_MIN"}

obsid_names={"XMM":"OBS_ID","Swift":"ObsID", "Stacked":"OBS_ID", "OM":"OBSID","UVOT":"OBSID"}

band_flux_names = {"XMM":["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"],
                  "Chandra":["flux_powlaw_aper_s","flux_powlaw_aper_m","flux_powlaw_aper_h"],
                  "Swift":["Flux1","Flux2","Flux3"],
                  "Slew":["Flux6","Flux7"],
                  "eRosita":['ML_FLUX_b1','ML_FLUX_b2','ML_FLUX_b3','ML_FLUX_b4'],
                  "RASS":["Flux1","Flux3","Flux4"],
                  "WGACAT":["Flux1","Flux2","Flux3"],
                  "Stacked":["EP_1_FLUX","EP_2_FLUX","EP_3_FLUX","EP_4_FLUX","EP_5_FLUX"],
                  "OM":["UVW2_AB_FLUX","UVM2_AB_FLUX","UVW1_AB_FLUX","U_AB_FLUX","B_AB_FLUX","V_AB_FLUX"],
                  "UVOT":["UVW2_FLUX","UVM2_FLUX","UVW1_FLUX","U_FLUX","B_FLUX","V_FLUX"]}
band_fluxerr_names = {"XMM":[["EP_1_FLUX_ERR","EP_2_FLUX_ERR","EP_3_FLUX_ERR","EP_4_FLUX_ERR","EP_5_FLUX_ERR"],
                            ["EP_1_FLUX_ERR","EP_2_FLUX_ERR","EP_3_FLUX_ERR","EP_4_FLUX_ERR","EP_5_FLUX_ERR"]],
                     "Chandra":[["flux_powlaw_aper_s_negerr","flux_powlaw_aper_m_negerr","flux_powlaw_aper_h_negerr"],
                                ["flux_powlaw_aper_s_poserr","flux_powlaw_aper_m_poserr","flux_powlaw_aper_h_poserr"]],
                     "Swift":[["FluxErr1_neg","FluxErr2_neg","FluxErr3_neg"],["FluxErr1_pos","FluxErr2_pos","FluxErr3_pos"]],
                     "Slew":[["Flux6Err","Flux7Err"],["Flux6Err","Flux7Err"]],
                     "eRosita":[['ML_FLUX_ERR_b1','ML_FLUX_ERR_b2','ML_FLUX_ERR_b3','ML_FLUX_ERR_b4'],
                                      ['ML_FLUX_ERR_b1','ML_FLUX_ERR_b2','ML_FLUX_ERR_b3','ML_FLUX_ERR_b4']],
                     "RASS":[["FluxErr1","FluxErr3","FluxErr4"],["FluxErr1","FluxErr3","FluxErr4"]],
                     "WGACAT":[["FluxErr1","FluxErr2","FluxErr3"],["FluxErr1","FluxErr2","FluxErr3"]],
                     "Stacked": [["EP_1_FLUX_ERR", "EP_2_FLUX_ERR", "EP_3_FLUX_ERR", "EP_4_FLUX_ERR", "EP_5_FLUX_ERR"],
                              ["EP_1_FLUX_ERR", "EP_2_FLUX_ERR", "EP_3_FLUX_ERR", "EP_4_FLUX_ERR", "EP_5_FLUX_ERR"]],
                     "OM":["UVW2_AB_FLUX_ERR","UVM2_AB_FLUX_ERR","UVW1_AB_FLUX_ERR","U_AB_FLUX_ERR","B_AB_FLUX_ERR","V_AB_FLUX_ERR"],
                     "UVOT":["UVW2_FLUX_ERR","UVM2_FLUX_ERR","UVW1_FLUX_ERR","U_FLUX_ERR","B_FLUX_ERR","V_FLUX_ERR"]
                      }
band_center = {"XMM":[0.35,0.75,1.5,3.25,8.25],
               "NewXMM":[0.35,0.75,1.5,3.25,8.25],
               "Chandra":[0.85,1.6,4.5],
               "Swift":[0.65,1.5,6],
               "Slew":[1.1,7],
               "eRosita":[0.35,0.75,1.5,3.25],
               "RASS":[0.25,0.7,1.45],
               "WGACAT":[0.25,0.7,1.45],
               "Stacked":[0.35,0.75,1.5,3.25,8.25],
               "OM":[2000,2300,2750,3500,4500,5500],
               "UVOT":[2000,2300,2750,3500,4500,5500]}
band_half_width = {"XMM":[0.15,0.25,0.5,1.25,3.75],
              "NewXMM":[0.15,0.25,0.5,1.25,3.75],
              "Chandra":[0.35,0.4,2.5],
              "Swift":[0.35,0.5,4],
              "Slew":[0.9,5],
              "eRosita":[0.15,0.25,0.5,1.25],
              "RASS":[0.15,0.2,0.55],
              "WGACAT":[0.15,0.2,0.55],
              "Stacked":[0.15,0.25,0.5,1.25,3.75],
              "OM":[100,200,250,500,500,500],
              "UVOT":[100,200,250,500,500,500]}

band_edges = {}
frequencies = {}
frequencies_half_width = {}
for cat in band_center.keys():
    band_edges[cat] = [center-width for (center, width) in zip(band_center[cat],band_half_width[cat])]
    band_edges[cat].append(band_center[cat][-1]+band_half_width[cat][-1])
for cat in catalogs:
    frequencies[cat]=[2.41e17*center for center in band_center[cat]]
    frequencies_half_width[cat] = [2.41e17*width for width in band_half_width[cat]]
xband_average_frequency = 2 * 2.41e17 #2keV in Hz, to compute alpha_OX
xband_width = 11.9 * 2.41e17 #11.9 keV in Hz, to compute alpha_OX

hr_bandlimit_index = {"XMM":3,"NewXMM":3,"Chandra":2,"Swift":2,"Slew":1,"eRosita":3,"RASS":3,"WGACAT":3,"Stacked":3}
band_conv_factors_soft = {"XMM":0.35/0.35,
                          "NewXMM":0.35/0.35,
                          "Chandra":0.35/0.28,
                          "Swift":0.35/0.34,
                          "Slew":0.35/0.35,
                          "eRosita":0.35/0.35,
                          "RASS":0.35/0.35,
                          "WGACAT":0.35/0.35,
                          "Stacked":0.35/0.35}
band_conv_factors_hard = {"XMM":0.65/0.65,
                          "NewXMM":0.65/0.65,
                          "Chandra":0.65/0.41,
                          "Swift":0.65/0.56,
                          "Slew":0.65/0.65,
                          "eRosita":0.65/0.25,
                          "RASS":np.nan,
                          "WGACAT":np.nan,
                          "Stacked":0.65/0.65}

hr_track_markers = {"XMM":"o","NewXMM":'o',"Chandra":"v","Swift":"s","Slew":"P","eRosita":"^","RASS":"d","WGACAT":"d","Stacked":"*"}
# endregion

# cleanness='CleanCatalogs/'
# path_to_catalogs = f"/home/erwan/Documents/PhD/LongTermVariability/LongTermVariability_Python/NewMatchMethod/{cleanness}Catalogs/FullData/"

class Source:
    """
    A Source object corresponds to a source from one of the X-ray catalogs. It has several attributes:
    - catalog: the corresponding catalog name, in the same naming convention as the catalog Table defined at the top
    - iau_name: the name of the source, considered as a unique identifier
    - fluxes: a table containing the fluxes extrapolated in the 0.1-12keV band
    - flux_errors: a table containing the 1 sigma flux errors extrapolated in the 0.1-12keV band. It consists of 2 tables,
    one for negative errors and one for positive errors.
    - timesteps: a table containing the MJD dates of detections
    - obsids: a table containing the ObsID for each detection, in the case of XMM and Swift (used for matching with OM & UVOT)
    - band_flux and band_fluxerr: same as fluxes and flux_errors, but is divided in the various detection bands of each instrument.
    Used for spectrum & SED plotting

    In the end, each Source will be associated to a unique MasterSource, each MasterSource having Source objects from several distinct catalogs
    """
    def __init__(self, catalog, iau_name, flux, fluxerr, timesteps, band_flux, band_fluxerr, obsids=[],swift_stacked_flux=[], swift_stacked_flux_err=[[],[]], swift_stacked_times=[[],[]], xmm_offaxis=[], short_term_var=[]):
        self.catalog = catalog
        self.name = iau_name
        self.master_source = []
        self.fluxes = flux
        self.flux_errors=fluxerr
        self.timesteps=[float(elt) for elt in timesteps]
        self.obsids=[int(obsid) for obsid in obsids]

        self.band_flux = band_flux
        self.band_fluxerr = band_fluxerr

        self.soft_dets = [np.sum(det[:hr_bandlimit_index[catalog]])*band_conv_factors_soft[catalog] for det in self.band_flux]
        self.soft_errors = [[np.sum(err_neg[:hr_bandlimit_index[catalog]])*band_conv_factors_soft[catalog] for err_neg in self.band_fluxerr[0]],
                            [np.sum(err_pos[:hr_bandlimit_index[catalog]])*band_conv_factors_soft[catalog] for err_pos in self.band_fluxerr[1]]]
        if catalog!= "RASS" and catalog!="WGACAT":
            self.hard_dets = [np.sum(det[hr_bandlimit_index[catalog]:])*band_conv_factors_hard[catalog] for det in self.band_flux]
            self.hard_errors = [
                [np.sum(err_neg[hr_bandlimit_index[catalog]:]) * band_conv_factors_hard[catalog] for err_neg in
                 self.band_fluxerr[0]],
                [np.sum(err_pos[hr_bandlimit_index[catalog]:]) * band_conv_factors_hard[catalog] for err_pos in
                 self.band_fluxerr[1]]]
        else:
            self.hard_dets = [np.nan for det in self.fluxes]
            self.hard_errors = [[np.nan for det in self.fluxes],[np.nan for det in self.fluxes]]


        self.hardness = [(hard-soft)/(hard+soft) for (soft,hard) in zip(self.soft_dets, self.hard_dets)]
        low_soft = np.where(np.array(self.soft_dets) - np.array(self.soft_errors[0]) < 0, 0,
                            np.array(self.soft_dets) - np.array(self.soft_errors[0]))
        low_hard = np.where(np.array(self.hard_dets) - np.array(self.hard_errors[0]) < 0, 0,
                            np.array(self.hard_dets) - np.array(self.hard_errors[0]))
        up_soft = np.where(np.array(self.soft_dets) + np.array(self.soft_errors[1]) < 0, 0,
                           np.array(self.soft_dets) + np.array(self.soft_errors[1]))
        up_hard = np.where(np.array(self.hard_dets) + np.array(self.hard_errors[1]) < 0, 0,
                           np.array(self.hard_dets) + np.array(self.hard_errors[1]))
        self.hardness_err = [[hr - (hard-soft)/(hard+soft) for (soft,hard,hr) in zip(up_soft, low_hard, self.hardness)],
                             [(hard-soft)/(hard+soft) - hr for (soft,hard,hr) in zip(low_soft, up_hard, self.hardness)]]
        self.swift_stacked_flux = swift_stacked_flux
        self.swift_stacked_flux_err = swift_stacked_flux_err
        self.swift_stacked_times = swift_stacked_times
        self.swift_stacked_variable = False

        self.min_upper = 1
        self.max_lower = 0
        self.var = 1
        if len(flux)>0:
            self.min_upper = min(np.array(flux) + np.array(fluxerr[1]))
            self.max_lower = max(np.array(flux) - np.array(fluxerr[0]))
        if swift_stacked_flux!=[]:
            stacked_min = min(np.array(swift_stacked_flux)+np.array(swift_stacked_flux_err[1]))
            if stacked_min<0.5*self.min_upper:
                self.swift_stacked_variable = True
            self.min_upper = min(self.min_upper, stacked_min)
        if len(flux)+len(swift_stacked_flux) > 1:
            self.var = self.max_lower/self.min_upper

        self.xmm_offaxis = xmm_offaxis
        self.short_term_var = short_term_var

class MasterSource:
    """
    A MasterSource corresponds to a single physical source, built on the association of multiple archival catalogs.
    A MasterSource has several attributes:
    - source: A dictionary which gives access to the underlying catalogs sources, which are Source objects in our framework.
    The keys of this dictionary are the names of the corresponding catalogs.
    - source_fluxes and source_error_bars: compiled flux and flux_errors from all its constituting Source objects, it is
    useful in order to access overall flux properties of the MasterSource (maximum flux, variability,...).
    - tab_hr and tab_hr_err: same thing but for hardness properties instead of flux properties.
    - var_ratio, var_amplitude, var_significance: correspond to different ways of quantifying flux variability of the source
    - hr_var, hr_var_signif: same thing but for hardness properties instead of flux properties

    A MasterSource only has one method, plot_lightcurve(), which produces a multi-panel plot of all relevant information
    """
    def __init__(self, id, tab_sources, ra, dec, poserr):
        self.id = id
        self.sources = {}
        self.sources_fluxes = []
        self.sources_error_bars = [[],[]]
        self.sources_timesteps = []
        self.sources_var = []
        self.tab_hr = []
        self.tab_hr_err = [[],[]]
        self.never_on_axis_xmm = False
        self.has_short_term_var = False
        self.min_time=60000
        self.max_time=0
        for source in tab_sources:
            if ("XMM" in self.sources.keys()) and (source.catalog == "Stacked"):
                #We remove the Stacked detection that correspond to a clean XMM detection
                xmm_obsid = self.sources["XMM"].obsids
                stacked_obsid = source.obsids
                new_det_ind = [i for i in range(len(stacked_obsid)) if stacked_obsid[i] not in xmm_obsid]
                source.fluxes = source.fluxes[new_det_ind]
                source.flux_errors[0] = source.flux_errors[0][new_det_ind]
                source.flux_errors[1] = source.flux_errors[1][new_det_ind]
                source.timesteps = np.array(source.timesteps)[new_det_ind]
                source.obsids = np.array(source.obsids)[new_det_ind]
                source.hardness = np.array(source.hardness)[new_det_ind]
                source.hardness_err[0] = np.array(source.hardness_err[0])[new_det_ind]
                source.hardness_err[1] = np.array(source.hardness_err[1])[new_det_ind]

                source.band_flux = source.band_flux[new_det_ind]
                source.band_fluxerr[0] = source.band_fluxerr[0][new_det_ind]
                source.band_fluxerr[1] = source.band_fluxerr[1][new_det_ind]
            source.master_source = self
            self.sources[source.catalog]=source
            for (flux, fluxerr_neg, fluxerr_pos, timestep) in zip(source.fluxes, source.flux_errors[0], source.flux_errors[1], source.timesteps):
                self.sources_fluxes.append(flux)
                self.sources_error_bars[0].append(fluxerr_neg)
                self.sources_error_bars[1].append(fluxerr_pos)
                self.sources_var.append(source.var)
                self.sources_timesteps.append(timestep)
            self.tab_hr += list(source.hardness)
            self.tab_hr_err[0] += list(source.hardness_err[0])
            self.tab_hr_err[1] += list(source.hardness_err[1])
            for (flux, fluxerr_neg, fluxerr_pos, start, stop) in zip(source.swift_stacked_flux, source.swift_stacked_flux_err[0], source.swift_stacked_flux_err[1], source.swift_stacked_times[0], source.swift_stacked_times[1]):
                self.sources_fluxes.append(flux)
                self.sources_error_bars[0].append(fluxerr_neg)
                self.sources_error_bars[1].append(fluxerr_pos)
                self.min_time = min(start, self.min_time)
                self.max_time = max(stop, self.max_time)
                self.sources_timesteps.append((start+stop)/2)
            if source.xmm_offaxis!=[]:
                if np.nanmin(source.xmm_offaxis)>1:
                    self.never_on_axis_xmm = True
            if source.timesteps!=[]:
                self.min_time = min(min(source.timesteps), self.min_time)
                self.max_time = max(max(source.timesteps), self.max_time)
            for var_flag in source.short_term_var:
                if var_flag>0:
                    self.has_short_term_var=True
        self.sources_fluxes = np.array(self.sources_fluxes)
        self.sources_error_bars = np.array(self.sources_error_bars)


        self.min_upper = 1
        self.max_lower = 0
        self.var_ratio = 1
        self.var_amplitude = 0
        self.var_significance = 0
        if len(self.sources_fluxes)>0 and (not np.isnan(self.sources_fluxes).all()):
            min_upper_ind = np.argmin(self.sources_fluxes + self.sources_error_bars[1])
            self.min_upper = (self.sources_fluxes + self.sources_error_bars[1])[min_upper_ind]
            max_lower_tab = np.where(self.sources_fluxes - self.sources_error_bars[0]>0,
                                     self.sources_fluxes - self.sources_error_bars[0],
                                     self.sources_fluxes)
            max_lower_ind = np.argmax(max_lower_tab)
            self.max_lower = max_lower_tab[max_lower_ind]
            self.var_ratio = self.max_lower/self.min_upper
            self.var_amplitude = self.max_lower - self.min_upper
            self.var_optimistic = self.sources_fluxes[max_lower_ind]/self.sources_fluxes[min_upper_ind]
            self.var_significance = self.var_amplitude/np.sqrt(self.sources_error_bars[1][max_lower_ind]**2 + self.sources_error_bars[0][min_upper_ind]**2)
            #self.frac_var = np.sqrt((np.var(self.sources_fluxes, ddof=1)-np.mean(np.array(self.sources_error_bars)**2))/(np.mean(self.sources_fluxes)**2))

        self.hr_min = np.nan
        self.hr_max = np.nan
        self.hr_var = np.nan
        self.hr_var_signif = np.nan
        if len(self.tab_hr)>1 and (not np.isnan(self.tab_hr).all()) and (not np.isnan(self.tab_hr_err).all()):
            #print(self.tab_hr, self.tab_hr_err)
            index_hr_min = np.nanargmin(np.array(self.tab_hr)+np.array(self.tab_hr_err[1]))
            index_hr_max = np.nanargmax(np.array(self.tab_hr)-np.array(self.tab_hr_err[0]))
            self.hr_min = (np.array(self.tab_hr)+np.array(self.tab_hr_err[1]))[index_hr_min]
            self.hr_max = (np.array(self.tab_hr)-np.array(self.tab_hr_err[0]))[index_hr_max]
            self.hr_var = self.hr_max - self.hr_min
            if self.tab_hr_err[1][index_hr_min]**2 + self.tab_hr_err[0][index_hr_max]**2 >0:
                self.hr_var_signif = self.hr_var/np.sqrt(self.tab_hr_err[1][index_hr_min]**2 + self.tab_hr_err[0][index_hr_max]**2)
            else:
                self.hr_var_signif = np.nan

        self.xmm_ul = []
        self.xmm_ul_dates = []
        self.xmm_ul_obsids = []

        self.slew_ul = []
        self.slew_ul_dates = []
        self.slew_ul_obsids = []

        self.chandra_ul = []
        self.chandra_ul_dates = []

        self.ra = float(ra)
        self.dec = float(dec)
        self.pos_err = float(poserr)

        self.glade_distance=[]

        self.simbad_type=''
        self.has_sdss_widths=False

    def plot_lightcurve(self, limit_dates=[]):
        plt.rcParams.update({'font.size': 8})
        fig, [ax1, ax2, ax3] = plt.subplots(1,3, figsize=(15,5))

        plt.suptitle(f'Simbad type: {self.simbad_type}   -   More details', picker=True, bbox=dict(facecolor=(180 / 256., 204 / 256., 252 / 256.)))
        xmm_name=[]
        swift_name=[]
        if len(self.xmm_ul)!= 0:
            ax1.errorbar(self.xmm_ul_dates, self.xmm_ul, yerr=0.2 * np.array(self.xmm_ul), uplims=True, fmt='none', c=colors["XMM"], label="XMM non-det.")
        if len(self.slew_ul)!= 0:
            ax1.errorbar(self.slew_ul_dates, self.slew_ul, yerr=0.2 * np.array(self.slew_ul), uplims=True, fmt='none', c=colors["Slew"], label="Slew non-det.")
        if len(self.chandra_ul)!= 0:
            ax1.errorbar(self.chandra_ul_dates, self.chandra_ul, yerr=0.2 * np.array(self.chandra_ul), uplims=True, fmt='none', c=colors["Chandra"])


        hardness_track=[]
        hardness_err_track=[[],[]]
        luminosity_track=[]
        luminosity_err_track = [[], []]
        time_track=[]
        catalogs_track=[]
        for cat in catalogs:
            if cat in self.sources.keys():
                source = self.sources[cat]
                ax1.errorbar(np.array(source.timesteps), np.array(source.fluxes),
                            yerr=np.array(source.flux_errors), fmt="o",c=colors[cat],
                            label=source.name, markeredgecolor='gray')
                if cat == "Swift":
                    times=[(stop+start)/2 for (start,stop) in zip(source.swift_stacked_times[0],source.swift_stacked_times[1])]
                    timerange=[(stop-start)/2 for (start,stop) in zip(source.swift_stacked_times[0],source.swift_stacked_times[1])]
                    ax1.errorbar(times, source.swift_stacked_flux,
                                yerr=source.swift_stacked_flux_err, xerr=timerange,
                                 fmt="o", markeredgecolor='gray', c=colors[cat])
                tab_width = 2*np.array(band_half_width[cat])
                for det in range(len(source.band_flux)):
                    ax2.step(band_edges[cat], [source.band_flux[det][0]/tab_width[0]] + list(source.band_flux[det]/tab_width), c=colors[cat], where='pre')
                    ax2.errorbar(band_center[cat], source.band_flux[det]/tab_width,
                                 yerr=[source.band_fluxerr[0][det]/tab_width,source.band_fluxerr[1][det]/tab_width],
                                 fmt="o", markeredgecolor='gray', c=colors[cat], alpha=0.4)
                hardness_track+=list(source.hardness)
                hardness_err_track[0] += list(source.hardness_err[0])
                hardness_err_track[1] += list(source.hardness_err[1])
                luminosity_track+=list(source.fluxes)
                luminosity_err_track[0] += list(source.flux_errors[0])
                luminosity_err_track[1] += list(source.flux_errors[1])
                time_track+=list(source.timesteps)
                catalogs_track+=[cat for elt in source.timesteps]
                if cat=="XMM":
                    xmm_name = source.name
                elif cat=="Chandra":
                    chandra_name = source.name
                elif cat=="Swift":
                    swift_name = source.name

        if limit_dates!=[]:
            ax1.axvline(limit_dates[0], linestyle='--')
            ax1.axvline(limit_dates[1], linestyle='--')
        has_eFeds = "eRosita" in self.sources.keys()
        fig.canvas.mpl_connect('pick_event', lambda event: click_action(self.ra, self.dec, xmm_name, swift_name, has_eFeds, self.has_sdss_widths or self.has_sdss))

        order = np.argsort(time_track)
        hardness_track=np.array(hardness_track)[order]
        luminosity_track = np.array(luminosity_track)[order]
        hardness_err_track=[np.array(hardness_err_track[0])[order],np.array(hardness_err_track[1])[order]]
        luminosity_err_track=[np.array(luminosity_err_track[0])[order],np.array(luminosity_err_track[1])[order]]
        time_track = np.array(time_track)[order]
        catalogs_track = np.array(catalogs_track)[order]
        color_track = np.array([matplotlib.cm.get_cmap('inferno')((time-time_track[0])/(time_track[-1]-time_track[0])) for time in time_track])
        ax3.errorbar(hardness_track,luminosity_track,xerr=hardness_err_track,yerr=luminosity_err_track, alpha=0.2, linestyle="--")
        #ax3.scatter(hardness_track,luminosity_track, c=ax3.lines[-1].get_color())
        for cat in catalogs:
            ax3.scatter(hardness_track[catalogs_track==cat], luminosity_track[catalogs_track==cat], c=ax3.lines[-1].get_color(),marker=hr_track_markers[cat], s=50, label=cat, edgecolors="gray")

        #ax1.tick_params(axis='x', rotation=45)
        ax1.set_title("Long-term lightcurve (0.2-12 keV)")
        ax1.legend(loc="best")
        ax1.set_yscale('log')
        ax1.set_xlabel("Time (MJD)")
        ax1.set_ylabel(r"Flux ($erg.s^{-1}.cm^{-2}$)")
        margin = 0.1*(self.max_time-self.min_time)
        ax1.set_xlim(self.min_time-margin, self.max_time+margin)

        ax2.set_title("Detections X-ray spectra")
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.set_xlabel("Energy (keV)")
        ax2.set_ylabel(r"$F_{\nu}$ ($erg.s^{-1}.cm^{-2}.keV^{-1}$)")

        ax3.set_title("Hardness-Luminosity diagram")
        ax3.legend(loc="best")
        ax3.set_yscale('log')
        ax3.set_xlabel("Hardness")
        ax3.set_xlim((-1.1,1.1))
        ax3.set_ylabel(r"Flux ($erg.s^{-1}.cm^{-2}$)")

        plt.draw()
        plt.tight_layout()
        plt.show()

#plt.rcParams['text.usetex'] = True

def load_relevant_source(cat,file_to_load):
    print(f"Loading {cat}...")
    raw_data = fits.open(file_to_load, memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    sources_raw = sources_raw[np.argsort(sources_raw[src_names[cat]])]

    indices_for_source = [i for i in range(1, len(sources_raw)) if (sources_raw[src_names[cat]][i] != sources_raw[src_names[cat]][i - 1])]

    if cat == "Swift":
        #conv_factors[cat] = np.where(np.isnan(sources_raw["FixedPowECFO"]), 4.2e-11,sources_raw["FixedPowECFO"]) * conv_factors[cat]
        timestartobs = Time(sources_raw["StartTime_UTC"], format="iso").mjd
        timeendobs = Time(sources_raw["StopTime_UTC"], format="iso").mjd
        timestartobs = np.split(timestartobs, indices_for_source)
        timeendobs = np.split(timeendobs, indices_for_source)


    #We divide up the catalog in sub-samples corresponding to each source
    timesteps = np.split(np.array(sources_raw[time_names[cat]]), indices_for_source)
    if cat in ("XMM","Swift","Stacked"):
        obsids = np.split(np.array(sources_raw[obsid_names[cat]]), indices_for_source)
    else:
        obsids = [[] for elt in indices_for_source]
    names = np.split(np.array(sources_raw[src_names[cat]]), indices_for_source)

    band_fluxes = []
    band_flux_errors_neg=[]
    band_flux_errors_pos=[]

    fluxes = np.split(conv_factors[cat]*np.array(sources_raw[flux_names[cat]]), indices_for_source)
    flux_errors_neg = np.split(conv_factors[cat]*np.array(sources_raw[flux_error_names[cat][0]]), indices_for_source)
    flux_errors_pos = np.split(conv_factors[cat]*np.array(sources_raw[flux_error_names[cat][1]]), indices_for_source)
    flux_errors = [[flux_neg, flux_pos] for (flux_neg, flux_pos) in zip(flux_errors_neg, flux_errors_pos)]

    for band_flux_name, band_fluxerr_neg_name, band_fluxerr_pos_name in zip(band_flux_names[cat],
                                                                          band_fluxerr_names[cat][0],band_fluxerr_names[cat][1]):
        band_fluxes.append(np.array(sources_raw[band_flux_name]))
        band_flux_errors_neg.append(np.array(sources_raw[band_fluxerr_neg_name]))
        band_flux_errors_pos.append(np.array(sources_raw[band_fluxerr_pos_name]))
    band_fluxes = np.transpose(np.array(band_fluxes))
    band_flux_errors_neg = np.transpose(np.array(band_flux_errors_neg))
    band_flux_errors_pos = np.transpose(np.array(band_flux_errors_pos))
    band_fluxes = np.split(band_fluxes, indices_for_source)
    band_flux_errors_neg = np.split(band_flux_errors_neg, indices_for_source)
    band_flux_errors_pos = np.split(band_flux_errors_pos, indices_for_source)
    band_fluxerr = [[band_flux_neg, band_flux_pos] for (band_flux_neg, band_flux_pos) in zip(band_flux_errors_neg, band_flux_errors_pos)]


    dic_sources = {}

    #This loops on all sources, to build the Source objects
    for (index, flux, flux_error, time, name, band_flux, band_fluxerr, obsid) in zip(range(len(fluxes)),fluxes, flux_errors, timesteps, names, band_fluxes, band_fluxerr, obsids):
            swift_stacked_flux=[]
            swift_stacked_flux_err=[[],[]]
            swift_stacked_times=[[],[]]
            if cat == "Swift":
                tab_src_timestartobs = timestartobs[index]
                tab_src_timeendobs = timeendobs[index]

                #We select the stacked Swift detections first
                swift_stacked_flux=flux[obsid>1e10]
                swift_stacked_flux_err=[flux_error[0][obsid>1e10],flux_error[1][obsid>1e10]]
                swift_stacked_times=[tab_src_timestartobs[obsid>1e10], tab_src_timeendobs[obsid>1e10]]

                # We then treat the classical, non-stacked Swift detections
                flux = flux[obsid < 1e10]
                flux_error = [flux_error[0][obsid < 1e10], flux_error[1][obsid < 1e10]]
                time = time[np.where(obsid < 1e10)]
                band_flux = band_flux[obsid < 1e10]
                band_fluxerr = [band_fluxerr[0][obsid < 1e10], band_fluxerr[1][obsid < 1e10]]
                obsid = obsid[obsid < 1e10]
            source = Source(cat, name[0].strip(), flux, flux_error, time, band_flux, band_fluxerr, obsid, swift_stacked_flux,swift_stacked_flux_err,swift_stacked_times)
            dic_sources[name[0].strip()] = source
    return dic_sources

def load_XMM_upperlimits(dic_master_sources):
    raw_data = fits.open(f"{path_to_master_sources}Master_source_XMM_UpperLimits.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    xmm_pointed_ul = sources_raw[sources_raw["obstype"]=="pointed"]
    dates = Time(xmm_pointed_ul["MJD_Date"],format="mjd").mjd
    print("Loading XMM pointed upper limits on Master Sources...")
    pbar = tqdm(total=len(dates))
    for line, date in zip(xmm_pointed_ul, dates):
        ms = dic_master_sources[line["MS_ID"]]
        if ("Stacked" not in ms.sources.keys()) or (int(line["obsid"]) not in ms.sources["Stacked"].obsids):
            #if int(line["obsid"]) not in ms.sources["Stacked"].obsids:
                #In the case of a simultaneous Stacked detection / RapidXMM upper limit, we take the Stacked detection
                ms.xmm_ul.append(line["ul_flux8_3sig"])
                ms.xmm_ul_dates.append(date)
                ms.xmm_ul_obsids.append(line["obsid"])
                ms.min_time = min(date, ms.min_time)
                ms.max_time = max(date, ms.max_time)
        """else:
                ms.xmm_ul.append(line["ul_flux8_3sig"])
                ms.xmm_ul_dates.append(date)
                ms.xmm_ul_obsids.append(line["obsid"])
                ms.min_time = min(date, ms.min_time)
                ms.max_time = max(date, ms.max_time)"""
        pbar.update(1)
    pbar.close()

    print("Loading XMM slew upper limits on Master Sources...")
    xmm_slew_ul = sources_raw[sources_raw["obstype"] == 'slew   ']
    dates = Time(xmm_slew_ul["MJD_Date"], format="mjd").mjd
    pbar = tqdm(total=len(dates))
    for line, date in zip(xmm_slew_ul, dates):
        ms = dic_master_sources[line["MS_ID"]]
        ms.slew_ul.append(line["ul_flux8_3sig"])
        ms.slew_ul_dates.append(date)
        ms.slew_ul_obsids.append(line["obsid"])
        ms.min_time = min(date, ms.min_time)
        ms.max_time = max(date, ms.max_time)
        pbar.update(1)
    pbar.close()


    """#Updating variabilities using upper limits
    tab_improvement_withxxm = []
    tab_improvement_withoutxxm = []
    to_plot = []
    improve=[]
    for ms in dic_master_sources.values():
        if len(ms.xmm_ul)>0:
            if ms.min_upper <1 and min(ms.xmm_ul) > 0:
                if "XMM" in ms.sources.keys() or "Stacked" in ms.sources.keys():
                    tab_improvement_withxxm.append(ms.min_upper / min(ms.xmm_ul))
                else:
                    tab_improvement_withoutxxm.append(ms.min_upper / min(ms.xmm_ul))
                    if (ms.min_upper / min(ms.xmm_ul)>10):
                        to_plot.append(ms)
                        improve.append(ms.min_upper / min(ms.xmm_ul))
            #ms.min_upper = min(ms.min_upper, min(ms.xmm_ul))
        #if len(ms.slew_ul)>0:
        #    ms.min_upper = min(ms.min_upper, min(ms.slew_ul))
        ms.var = ms.max_lower/ms.min_upper"""
    """plt.hist(tab_improvement_withxxm, bins=np.geomspace(1e-3, 1e3, 50), color="royalblue", label="Had XMM detection")
    plt.hist(tab_improvement_withoutxxm, bins=np.geomspace(1e-3, 1e3, 50), color="r", label="No XMM detection")
    plt.xscale('log')
    plt.legend()
    plt.xlabel("Minimum detection / XMM-Newton upper limit")
    order = np.argsort(improve)[::-1]
    for ind in order[:10]:
        ms = to_plot[ind]
        ms.plot_lightcurve()"""

def load_Chandra_upperlimits(dic_master_sources):
    print("Load Chandra upper limits...")
    raw_data = fits.open(f"{path_to_master_sources}Chandra_UL.fits", memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)
    for line in tqdm(sources_raw):
        ms = dic_master_sources[line["MS_ID"]]
        ms.chandra_ul.append(line["flux_aper_hilim_b"])
        ms.chandra_ul_dates.append(line["gti_mjd_obs"])


def load_master_sources(file_to_load):
    """Loads the multi-instruments sources in a dictionary"""
    print(f"Loading Master Sources...")
    raw_data = fits.open(os.path.join(file_to_load, 'Master_source_cone.fits'), memmap=True)
    sources_raw = raw_data[1].data
    sources_raw = Table(sources_raw)

    tab_catalog_sources = {}
    for cat in catalogs:
        tab_catalog_sources[cat] = load_relevant_source(cat,os.path.join(file_to_load,cat+'.fits'))

    dic_master_sources = {}
    for line in tqdm(sources_raw):
        tab_sources_for_this_ms = []
        for cat in catalogs:
            if line[cat]!='':
                name=line[cat].strip()
                if name in tab_catalog_sources[cat].keys():
                    tab_sources_for_this_ms.append(tab_catalog_sources[cat][name])
        ms_id = line["MS_ID"]
        ms = MasterSource(ms_id, tab_sources_for_this_ms, line["MS_RA"], line["MS_DEC"], line["MS_POSERR"])
        dic_master_sources[ms_id] = ms
    #load_XMM_upperlimits(dic_master_sources)
    #load_Chandra_upperlimits(dic_master_sources)

    print("Master sources loaded!")
    return dic_master_sources


