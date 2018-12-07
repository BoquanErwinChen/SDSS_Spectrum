from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from matplotlib import pyplot as plt
import numpy as np

class SpectrumCosmology():
	
	def __init__(self, astropy_cosmology, z=None):
		self.cosmology = astropy_cosmology
		self.z = z

	@property
	def age(self):
		return(self.cosmology.age(self.z))

	@property
	def comoving_distance(self):
		return(self.cosmology.comoving_distance(self.z))

	@property
	def luminosity_distance(self):
		return(self.cosmology.luminosity_distance(self.z))


class Spectrum():
	def __init__(self, filename):
		self.filename = filename

		with fits.open(self.filename) as hdu_list:
		
			self.ra   = hdu_list['PRIMARY'].header['RA']  # degrees
			self.dec  = hdu_list['PRIMARY'].header['DEC'] # degrees			
			
			self.z    =  hdu_list['SPALL'].data['Z'][0]   # redshift
			
			flux      = hdu_list['COADD'].data['flux']
			self.flux = flux / np.nanmedian(flux)         # normalise
		
			# wavelength in units of Angstroms	
			self.wavelength = 10**hdu_list['COADD'].data['loglam']

			# standard cosmology
			self._cosmology = None

	@property
	def rest_wavelength(self):
		''' de-redshift the wavelength array '''		
		return(self.wavelength / (1 + self.z))

	@property
	def ra_arcsec(self):
		''' right ascension in arcseconds'''
		return(self.ra * 3600)

	@property
	def dec_arcsec(self):
		''' declination in arcseconds '''
		return(self.dec * 3600)


	def separation(self,s, unit='degree'):
		''' Return the angle on the sky between two objects. defult=degrees '''

		loc1 = SkyCoord(ra=self.ra, dec=self.dec, unit='deg')
		loc2 = SkyCoord(ra=s.ra, dec=s.dec, unit='deg')

		# get the separation in the desired units
		separation = getattr(loc1.separation(loc2), unit)

		return(separation)

	
	@property
	def cosmology(self):
		assert self.z > 0, "Redshift contaimnated by peculiar velocity, not a reliable distance indicator"
		if self._cosmology is None:
			self._cosmology = SpectrumCosmology(FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725), self.z)
			self._cosmology.z = self.z
		return self._cosmology



	def hist_flux(self, bins=10, label='Flux', fontsize=15, figname=None, **kwargs):
		"""
		Plot a histogram for flux.

		bins: int or sequence or str, optional
			Consistent with np.histogram.

		label: str, optional
			A label for the x axis.

		fontsize: int, optional
			The font size of x label.

		figname: str, optional
			If not None, the figure will be saved under the specified name.

		**kwargs:
			Other properties in matplotlib.pyplot.hist.
		"""

		fig, ax = plt.subplots()
		ax.hist(self.flux, bins=bins, **kwargs)
		ax.set_xlabel(label, fontsize=fontsize)
		if figname is not None:
			plt.savefig(f'{figname}.png', dpi=300)
		plt.show()


	def plot_wavelength_flux(self, xlabel='Wavelength', ylabel='Flux', fontsize=15, figname=None, **kwargs):
		"""
		Plot flux over wavelenth.

		xlabel: str, optional
			A label for the x axis.

		ylabel: str, optional
			A label for the y axis.

		fontsize: int, optional
			The font size of labels.

		figname: str, optional
			If not None, the figure will be saved under the specified name.

		**kwargs:
			Other properties in matplotlib.pyplot.plot.
		"""

		fig, ax = plt.subplots()
		ax.plot(self.rest_wavelength, self.flux, **kwargs)
		ax.set_xlabel(xlabel, fontsize=fontsize)
		ax.set_ylabel(ylabel, fontsize=fontsize)
		if figname is not None:
			plt.savefig(f'{figname}.png', dpi=300)
		plt.show()


def scatter_flux(spectra, xlabel='RA (deg)', ylabel='Dec (deg)', cblabel='Total Flux', fontsize=15, figname=None, **kwargs):
    """
    Plot sources on the sky colored by total flux.

    spectra: list of objects
        A list of Spectrum objects.

    xlabel: str, optional
        A label for the x axis.

    ylabel: str, optional
        A label for the y axis.

    cblabel: str, optional
        A label for color bar.

    fontsize: int, optional
        The font size of labels.

    figname: str, optional
        If not None, the figure will be saved under the specified name.

    **kwargs:
        Other properties in matplotlib.pyplot.plot.
    """

    fig, ax = plt.subplots()
    ra = [s.ra for s in spectra]
    dec = [s.dec for s in spectra]
    c = [sum(s.flux) for s in spectra]
    sca = ax.scatter(ra, dec, c=c, **kwargs)
    ax.invert_xaxis()
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    cb = fig.colorbar(sca, ax=ax)
    cb.set_label(cblabel)

    if figname is not None:
        plt.savefig(f'{figname}.png', dpi=300)
    plt.show()

def scatter_redshift(spectra, xlabel='RA (deg)', ylabel='Dec (deg)', cblabel='Redshift', fontsize=15, figname=None, **kwargs):
    """
    Plot sources on the sky colored by redshift.

    spectra: list of objects
        A list of Spectrum objects.

    xlabel: str, optional
        A label for the x axis.

    ylabel: str, optional
        A label for the y axis.

    cblabel: str, optional
        A label for color bar.

    fontsize: int, optional
        The font size of labels.

    figname: str, optional
        If not None, the figure will be saved under the specified name.

    **kwargs:
        Other properties in matplotlib.pyplot.plot.
    """

    fig, ax = plt.subplots()
    ra = [s.ra for s in spectra]
    dec = [s.dec for s in spectra]
    c = [s.z for s in spectra]
    sca = ax.scatter(ra, dec, c=c, **kwargs)
    ax.invert_xaxis()
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    cb = fig.colorbar(sca, ax=ax)
    cb.set_label(cblabel)
    if figname is not None:
        plt.savefig(f'{figname}.png', dpi=300)
    plt.show()







