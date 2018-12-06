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
		return(self.ra * 36000)

	@property
	def dec_arcsec(self):
		''' declination in arcseconds '''
		return(self.dec * 36000)


	def separation(self, s, unit='degree'):
		''' Return the angle on the sky between two objects. defult=degrees '''

		loc1 = SkyCoord(ra=self.ra, dec=self.dec, unit='deg')
		loc2 = SkyCoord(ra=s.ra, dec=s.dec, unit='deg')

		# get the separation in the desired units
		separation = getattr(loc1.separation(loc2), unit)

		return(separation)

	
	@property
	def cosmology(self):
		if self._cosmology is None:
			self._cosmology = SpectrumCosmology(FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725), self.z)
			self._cosmology.z = self.z
		return self._cosmology


	def hist_flux(self, bins=10, label='Flux', fontsize=15, **kwargs):
		"""
		Plot a histogram for flux.

        bins: int or sequence or str, optional
            Consistent with np.histogram.

        label: str, optional
            A label for the x axis.

        fontsize: int, optional
            The font size of x label.

        **kwargs:
            Other properties in matplotlib.pyplot.hist.
        """

		fig, ax = plt.subplots()
		ax.hist(self.flux, bins=bins, **kwargs)
		ax.set_xlabel(label, fontsize=fontsize)
		plt.show()


	def plot_wavelength_flux(self, xlabel='Wavelength', ylabel='Flux', fontsize=15, **kwargs):
		"""
		Plot flux over wavelenth.

        xlabel: str, optional
            A label for the x axis.

        ylabel: str, optional
            A label for the y axis.

        fontsize: int, optional
            The font size of labels.

        **kwargs:
            Other properties in matplotlib.pyplot.plot.
        """

		fig, ax = plt.subplots()
		ax.plot(self.rest_wavelength, self.flux, **kwargs)
		ax.set_xlabel(xlabel, fontsize=fontsize)
		ax.set_ylabel(ylabel, fontsize=fontsize)
		plt.show()



		







