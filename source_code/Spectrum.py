from astropy.io import fits
import numpy as np

class Spectrum():
	def __init__(self, filename):
		self.filename = filename

		with fits.open(self.filename) as hdu_list:
		
			self.ra   = hdu_list['PRIMARY'].header['RA'] 
			self.dec  = hdu_list['PRIMARY'].header['DEC']
			
			self.z    =  hdu_list[2].data['Z'][0] # redshift

			flux      = hdu_list['COADD'].data['flux']
			self.flux = flux / np.nanmedian(flux) # normalise
		
			# wavelength in units of Angstroms	
			self.wavelength = 10**hdu_list['COADD'].data['loglam']


	@property
	def rest_wavelength(self):
		return(self.wavelength / (1 + self.z)) # de-redshift the wavelength array


	


