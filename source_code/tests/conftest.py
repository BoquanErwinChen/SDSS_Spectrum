import pytest
from source_code.Spectrum import Spectrum
import glob

spectrum_paths = glob.glob('/home/osboxes/Scicoder/SDSS_Spectrum/spectra/*.fits')
#print(spectrum_paths)

@pytest.fixture(scope="module",params=spectrum_paths)
def spectrumfile(request):
 	return request.param

@pytest.fixture(scope="module")
def spec2file():
	sample_2 = "/home/osboxes/Scicoder/SDSS_Spectrum/spectra/spec-4055-55359-0001.fits"
	s2 = Spectrum(sample_2)
	return s2 

@pytest.fixture(scope="module")
def spec1file():
	sample_1 = "/home/osboxes/Scicoder/SDSS_Spectrum/spectra/spec-10000-57346-0002.fits"
	s1 = Spectrum(sample_1)
	return s1 

