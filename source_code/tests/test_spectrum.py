#from astropy.io import fits as pyfits
from source_code.Spectrum import Spectrum
import numpy as np

#



def test_spectrum_open(spectrumfile):
	''' Check that we can open a spectrum file and the file exists '''
	s1 = Spectrum(spectrumfile)
	assert s1 is not None, "File does is empty or is not read properly" 

#def test_ra(spectrumfile):
	''' Check that we can read header keyword RA from an SDSS file. '''
	np.testing.assert_approx_equal(spectrumfile.ra, 29.815837)

def test_ra_bound(spectrumfile):
	'''Check that RA is within acceptable bounds'''
	assert 0.0 <= spectrumfile.ra <= 360, "RA is out of bounds"

def test_dec_bound(spectrumfile):
	'''Check that RA is within acceptable bounds'''
	assert -90 <= spectrumfile.dec <= 90, "DEC is out of bounds"

def test_redshift_bounds(spectrumfile):
	'''Check that redshift is not negative since big peciliar velocities preclude use as a distance indicator'''
	assert spectrumfile.z > 0.0, "Redshift is negative!"

def test_flux_data(spectrumfile):
	'''Check that flux data is being read'''
	assert len(spectrumfile.flux)>0, "No flux data read"

def test_flux_real(spectrumfile):
	'''Check that flux array has meaningful data'''
	assert spectrumfile.flux[0] is not None, "Flux data is invalid: Nones are present"

def test_flux_real(spectrumfile):
	'''Check that flux array has meaningful data'''
	assert np.nanmedian(spectrumfile.flux) != 0, "Flux data is invalid: only zeroes"

def test_arcsec_transform(spectrumfile):
	'''Check ra transformation'''
	np.testing.assert_approx_equal(spectrumfile.ra_arcsec, 29.815837*3600)

#test that distance is true
def test_separation(spectrumfile, spec2file):
	'''Check sky point separation'''
	np.testing.assert_approx_equal(spectrumfile.separation(spec2file), 141.234249521342)

def test_separation(spectrumfile, spec2file):
	'''Check sky point separation is symmetric'''
	np.testing.assert_approx_equal(spectrumfile.separation(spec2file),spec2file.separation(spectrumfile))


#test that a plot object is returned

#

