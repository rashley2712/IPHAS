#!/usr/bin/env python

import argparse, sys
import datetime, time
import numpy
import ppgplot
import generalUtils, configHelper, os
import curses
from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.table import Table, vstack
from astropy.utils import data
import astropy.table

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Loads an IFAS reduced image. Displays it with PGPLOT.')
	parser.add_argument('filename', type=str, help='The FITS image file.')
	parser.add_argument('--save', action="store_true", help='Write the input parameters to the config file as default values.')
	parser.add_argument('--ignorecache', action="store_true", help='Ignore the cached DR2 catalogue. Will overwrite one if it already exists.')
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
	imageCentre = [ width/2, height/3]
	
	
	ra, dec = wcsSolution.all_pix2world([imageCentre], 1)[0]
	positionString = generalUtils.toSexagesimal((ra, dec))
	print ra, dec
	print "RA, DEC of image centre is: ", positionString
	
	filenameParts = args.filename.split('.')
	dr2Filename = filenameParts[0] + "_dr2_cache.fits"
	cached = False
	if not args.ignorecache:
		print "Looking for a cached copy of the DR2 catalogue:", dr2Filename
		if os.path.exists(dr2Filename):
			cached = True



	if not cached:
		for index, yCentre in enumerate(numpy.arange(height/4, height, height/4)):
			print yCentre
			imageCentre = [ width/2, yCentre]
			ra, dec = wcsSolution.all_pix2world([imageCentre], 1)[0]
			print index, imageCentre, ra, dec
			with data.conf.set_temp('remote_timeout', 30):
				search = conesearch(center=(ra, dec),
	                    radius=0.1,
	                    verb=3,
	                    catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
			dr2nearbyTemp = search.to_table()
			if index==0: 
				dr2nearby = dr2nearbyTemp
			else:
				dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
				dr2nearby = dr2nearbyWhole
			print dr2nearby
			
			
		dr2nearby.write(dr2Filename, format='fits', overwrite=True)
	
	else:
		dr2nearby = Table.read(dr2Filename)
	
	print dr2nearby.colnames
	
	dr2nearby.remove_column('errBits2')
	# print "Table length:", len(dr2nearby)
	# print "(Unique) table length:", len(astropy.table.unique(dr2nearby))
	
	
	raArray = dr2nearby['RAJ2000']
	decArray = dr2nearby['DEJ2000']
	classArray = dr2nearby['mergedClass']
	xArray, yArray = wcsSolution.all_world2pix(raArray, decArray, 1)
	pStarArray = dr2nearby['pStar']
	ppgplot.pgsci(3)
	ppgplot.pgsfs(2)
	index = 0
	for ra, dec, x, y, c, p in zip(raArray, decArray, xArray, yArray, classArray, pStarArray):
		# print index, ra, dec, x, y, c
		index+= 1
		colour = 1
		if c==-9: colour = 0   # Black  = Saturated
		if c==1: colour = 4    # Blue   = Galaxy
		if c==-3: colour = 5   # Cyan   = Probable Galaxy
		if c==-1: colour = 7   # Yellow = Star
		if c==-2: colour = 8   # Orange = Probable Star
		if c==0: colour = 2    # Red    = Noise
		ppgplot.pgsci(colour)
		ppgplot.pgcirc(x, y, 5 + (5*p))
		
	screen = curses.initscr()
	screen.keypad(1)
	
	try:
		while True: 
			event = screen.getch()
			if event == ord("q"): break
			elif event == ord(" "): 
				screen.addstr("The User Pressed The Space Bar")
				x = 0
				y= 0
				ch = None
				ppgplot.pgcurs(x, y)
				print x, y, ch
	except KeyboardInterrupt:
		curses.endwin()
		print "Ctrl-C pressed, but I dealt with it. "

	curses.endwin()

	ppgplot.pgclos()
	
	
	
