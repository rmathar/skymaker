#!/usr/bin/env python3

# Richard J. Mathar, 2025-06-26

'''
This python program has two command line arguments:
an existing, read+writable FITS file, and an existing
 readable ASCII file with cfitsio template syntax
 of introducing new FITS header cards into the primary HDU.

This works similar to a repeated call of modhead
of the HEASARC software in 
https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html 
looping over the lines of the header card  files.
'''

from astropy.io import fits
import sys
from pathlib import Path

if len(sys.argv) < 2 :
	print("Usage:", sys.argv[0], "fitsfile.fits hdrfile.asc") 
	sys.exit(1)

# collect name of the FITS file
# and try to open it.
ffile=sys.argv[1]

try:
	hdu = fits.open(ffile,mode='update')
except:
	print("Cannot open ",ffile) 
	sys.exit(1)

# collect th ename of the file with the header cards
# and try to open it.
hfile=Path(sys.argv[2])
if not hfile.is_file():
	print(hfile,"does not exist") 
	sys.exit(1)

# hduHdr = hdu.PrimaryHDU().header
hduHdr = hdu[0].header
afile=open(hfile,'r')
# loop over all lines in the file with header cards
# (this does not yet remove comment lines that start with the hash)
for aline in afile :
	# split the header line at the equal sign
	asplit = aline.split("=")
	# split the value and comment at the slash
	a2split = asplit[1].split("/")
	# test whether the value is a string or a number
	# to remove the quotes for numbers.
	try:
		val = float(a2split[0])
	except:
		# need to remove quotes from a string value because astropy adds them again
		val = a2split[0].replace('\'','')
	hduHdr.set( asplit[0], val,a2split[1].rstrip() )

afile.close()
hdu.close()
