#!/usr/bin/env python

import argparse, sys
import datetime, time
import numpy
import ppgplot
import generalUtils, configHelper
import curses
from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.table import Table

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Loads an IFAS reduced image. Displays it with PGPLOT.')
	parser.add_argument('filename', type=str, help='The FITS image file.')
	parser.add_argument('--save', action="store_true", help='Write the input parameters to the config file as default values.')
	args = parser.parse_args()
	print args
	
	config = configHelper.configClass("inspectIPHASImage")

	"""exposureTime = config.getProperty("ExposureTime")
	if args.exposureTime!=None:
		exposureTime = args.exposureTime
		config.ExposureTime = exposureTime
	if exposureTime == None:
		print "Please specify an exposure time or save one in the config file."
		sys.exit()
	"""
	
	if args.save:
		config.save()
	
	hdulist = fits.open(args.filename)
	
	print hdulist.info()
	
	
	for card in hdulist:
		print card.header.keys()
		print repr(card.header)
	
	imageData =  hdulist[1].data
	
	wcsSolution = WCS(hdulist[1].header)
	
	hdulist.close()
	
	(height, width) = numpy.shape(imageData)
	
	aspectRatio = float(height)/float(width)
	print aspectRatio
	
	""" Set up the PGPLOT windows """
	imagePlot = {}
	imagePlot['pgplotHandle'] = ppgplot.pgopen('/xs')
	ppgplot.pgpap(6, aspectRatio)
	ppgplot.pgsvp(0.0, 1.0, 0.0, 1.0)
	ppgplot.pgswin(0, width, 0, height)
	
	# ppgplot.pgenv(0., width,0., height, 1, -2)
	imagePlot['pgPlotTransform'] = [0, 1, 0, 0, 0, 1]
	
	boostedImage = generalUtils.percentiles(imageData, 20, 99)
	ppgplot.pggray(boostedImage, 0, width-1, 0, height-1, 0, 255, imagePlot['pgPlotTransform'])
	
	# Determine the RA, DEC of the centre of the image, using the WCS solution found in the FITS header
	imageCentre = [ width/2, height/2]
	
	
	ra, dec = wcsSolution.all_pix2world([imageCentre], 1)[0]
	positionString = generalUtils.toSexagesimal((ra, dec))
	print ra, dec
	print "RA, DEC of image centre is: ", positionString
	
	cached = True
	if not cached:
		search = conesearch(center=(ra, dec),
	                    radius=0.12,
	                    verb=3,
	                    catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
		dr2nearby = search.to_table()
		dr2nearby.write('iphas-data.fits', format='fits', overwrite=True)
	else:
		dr2nearby = Table.read('iphas-data.fits')
	
	print dr2nearby
	print dr2nearby.colnames
	
	raArray = dr2nearby['RAJ2000']
	decArray = dr2nearby['DEJ2000']
	xArray, yArray = wcsSolution.all_world2pix(raArray, decArray, 1)

	ppgplot.pgsci(3)
	ppgplot.pgsfs(2)
	index = 0
	for ra, dec, x, y in zip(raArray, decArray, xArray, yArray):
		print index, ra, dec, x, y
		index+= 1

		ppgplot.pgcirc(x, y, 10)
		
	ppgplot.pgclos()
	
	
	
