# Example code #

from source_code.Spectrum import Spectrum, scatter_flux, scatter_redshift
import os

s1 = Spectrum('spectra/spec-10000-57346-0002.fits')
s2 = Spectrum('spectra/spec-6055-56102-0008.fits')

print("Right ascension and Declination (degrees) \n" 
	  "----------------------------------------  \n"    
	  "{0:.2f}, {1:.2f} \n{2:.2f}, {2:.2f} \n".format(s1.ra, s1.dec,s2.ra, s2.dec))


print("Seperation in arcseconds \n"
	  "-----------------------  \n"
	  "{0:.2f} arcsec".format(s1.separation(s2, 'arcsec')))
print()

# load all spectra under a directory
directory_name = 'spectra/'
directory = os.fsencode(directory_name)
list_spectrum = []
for file in os.listdir(directory):
	filename = os.fsdecode(file)
	if filename.endswith(".fits"):
		print(filename)
		list_spectrum.append(Spectrum(directory_name+filename))
		continue
	else:
		continue

# each spectrum has its own plotting functions
list_spectrum[0].hist_flux(bins=50)
list_spectrum[0].plot_wavelength_flux()

# other parameters in matplotlib can be set through a dictionary of keyword arguments
param = {'cumulative': True, 'histtype': 'step'}
list_spectrum[0].hist_flux(bins=50, label='Flux', fontsize=15, figname=None, **param)
param = {'c': 'r', 'linestyle': '--'}
list_spectrum[0].plot_wavelength_flux(xlabel='Wavelength', ylabel='Flux', fontsize=15, figname=None, **param)

# plot multiple spectra on the sky color coded in attributes
scatter_redshift(list_spectrum)

# other parameters in matplotlib can be set through a dictionary of keyword arguments
param = {'s': 3, 'linewidth': None, 'alpha': 0.5}
scatter_flux(list_spectrum, xlabel='RA (deg)', ylabel='Dec (deg)', cblabel='Total Flux', fontsize=15, figname=None,
			 **param)