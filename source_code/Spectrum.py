from astropy.io import fits
from astropy.coordinates import SkyCoord
import numpy as np


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


	@property
	def rest_wavelength(self):
		''' de-redshift the wavelength array '''		
		return(self.wavelength / (1 + self.z))


	def separation(self, s):
		''' Return the angle on the sky between two objects '''
		loc1 = SkyCoord(ra=self.ra, dec=self.dec, unit='deg')
		loc2 = SkyCoord(ra=s.ra, dec=s.dec, unit='deg')

		separation = loc1.separation(loc2)
		return(separation)







