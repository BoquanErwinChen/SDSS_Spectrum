#from astropy.io import fits as pyfits
from source_code.Spectrum import Spectrum
import numpy as np

sample_1 = "/home/osboxes/Scicoder/SDSS_Spectrum/spectra/spec-10000-57346-0002.fits"
sample_2 = "/home/osboxes/Scicoder/SDSS_Spectrum/spectra/spec-10000-57346-0003.fits"

s = Spectrum(sample_1)


'''Some kind ob object I make that has attributes RA and dec, user specifiable, for testing edge cases'''


def test_spectrum_open():
	''' Check that we can open a spectrum file and the file exists '''
	assert s is not None, "File does is empty or is not read properly" 

def test_ra():
	''' Check that we can read header keyword RA from an SDSS file. '''
	np.testing.assert_approx_equal(s.ra, 29.815837)

#def test_ra_bounds():
#	'''Check that RA is within acceptable bounds'''


#test that ra and dec are within given bounds

#test to open a spectrum file

#test to read header keywords

#test to read data

#test for coordinate transform

#test that distance is true

#test that distance works for edge cases(?)

#test that distance is symmetric (i.e. dist 1 of 2 is the same as dist 2 of 2

#test that redshift is not negative (means peculiar velocity is large, and not usable as a redshift/distance indicator)

#test that a plot object is returned

#

