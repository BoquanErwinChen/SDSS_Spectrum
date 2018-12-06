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

# plot the sources on the sky color coded in attributes
scatter_flux(list_spectrum)
scatter_redshift(list_spectrum)
