# Example code #

from source_code.Spectrum import Spectrum

s1 = Spectrum('spectra/spec-10000-57346-0002.fits')
s2 = Spectrum('spectra/spec-6055-56102-0008.fits')

print("Right ascension and Declination (degrees) \n" 
	  "----------------------------------------  \n"    
	  "{0:.2f}, {1:.2f} \n{2:.2f}, {2:.2f} \n".format(s1.ra, s1.dec,s2.ra, s2.dec))


print("Seperation in arcseconds \n"
	  "-----------------------  \n"
	  "{0:.2f} arcsec \n".format(s1.separation(s2, 'arcsec')))

print("Age of Universe at that redshift \n"
	  "-------------------------------- \n"
	  "redshift : {0:.2f}, age : {1:.2f} \n".format(s1.z, s1.cosmology.age))



