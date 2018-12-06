#from astropy.io import fits as pyfits
from source_code.Spectrum import Spectrum
import numpy as np

#------------Checks on multiple spectra-------------
def test_spectrum_open(spectrumfile):
	#Check that we can open a spectrum file and the file exists 
	s1 = Spectrum(spectrumfile)
	assert s1 is not None, "File does is empty or is not read properly" 

def test_ra_bound(spectrumfile):
	#Check that RA is within acceptable bounds
	s1 = Spectrum(spectrumfile)
	assert 0.0 <= s1.ra <= 360, "RA is out of bounds"

def test_dec_bound(spectrumfile):
	#Check that RA is within acceptable bounds
	s1 = Spectrum(spectrumfile)
	assert -90 <= s1.dec <= 90, "DEC is out of bounds"

def test_redshift_bounds(spectrumfile):
	#Check that redshift is not negative since big peciliar velocities preclude use as a distance indicator
	s1 = Spectrum(spectrumfile)
	assert s1.z > 0.0, "Redshift is negative!"

def test_flux_data(spectrumfile):
	#Check that flux data is being read
	s1 = Spectrum(spectrumfile)
	assert len(s1.flux)>0, "No flux data read"

def test_flux_real(spectrumfile):
	#Check that flux array has meaningful data
	s1 = Spectrum(spectrumfile)
	assert s1.flux[0] is not None, "Flux data is invalid: Nones are present"

def test_flux_real(spectrumfile):
	#Check that flux array has meaningful data
	s1 = Spectrum(spectrumfile)
	assert np.nanmedian(s1.flux) != 0, "Flux data is invalid: only zeroes"


#-----------Test on specific spectra-------
def test_ra(spec1file):
	#Check that we can read header keyword RA from an SDSS file. 
	np.testing.assert_approx_equal(spec1file.ra, 29.815837)

def test_arcsec_transform(spec1file):
	#Check ra transformation
	np.testing.assert_approx_equal(spec1file.ra_arcsec, 29.815837*3600)

def test_separation(spec1file, spec2file):
	#Check sky point separation
	np.testing.assert_approx_equal(spec1file.separation(spec2file), 141.234249521342)

def test_separation_symmetry(spec1file, spec2file):
	#Check sky point separation is symmetric
	np.testing.assert_approx_equal(spec1file.separation(spec2file),spec2file.separation(spec1file))

def test_flux_sum(spec1file):
	#Check flux sum is accurate; important since used in making flux plot
	np.testing.assert_approx_equal(sum(spec1file.flux), 4539.678025157074)

