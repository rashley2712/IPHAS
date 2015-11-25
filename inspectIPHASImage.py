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

def plotCircles(objectTable, margins):
	raArray = objectTable['RAJ2000']
	decArray = objectTable['DEJ2000']
	classArray = objectTable['mergedClass']
	xArray, yArray = wcsSolution.all_world2pix(raArray, decArray, 1)
	pStarArray = objectTable['pStar']
	ppgplot.pgsci(3)
	ppgplot.pgsfs(2)
	index = 0
	print "Margins:", margins
	for ra, dec, x, y, c, p in zip(raArray, decArray, xArray, yArray, classArray, pStarArray):
		if ra > margins[1][0] and ra < margins[0][0] and dec<margins[0][1] and dec>margins[1][1]:
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
	print "%d sources plotted."%(index+1)
		
def withinMargins(table, namedColumns):
	print namedColumns
	return True
	

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
	
	paperSize = 6  # Paper size in inches
	imageMinMax = (0, 255)
	
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
	ppgplot.pgpap(paperSize, aspectRatio)
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
	print "RA, DEC of image centre is: ", positionString, ra, dec
	margins = wcsSolution.all_pix2world([[0, 0], [width, height]], 1)
	print "ra, dec limits:", margins
	
	
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
	xlimits = (0, width)
	ylimits = (0, height)
	
	plotCircles(dr2nearby, margins)
	
		
	# try: 
	x=width/2
	y=height/2
	ppgplot.pgsci(3)
	keyPressed = None
	while keyPressed != 'q':
		ch = ppgplot.pgcurs(x, y)
		x=ch[0]
		y=ch[1]
		keyPressed = ch[2]
		print "Key pressed:", ch[2]
		if keyPressed=='i':
			print "Zoom requested at (%0.0f, %0.0f)"%(x, y)
			zoomFactor = 1.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			if (newWidth<100) or (newHeight<100):
				print "Maximum zoom reached."
				continue
			else:
				xlimits = (int(x - newWidth/2), int(x + newWidth/2)) 
				if xlimits[0] < 0:
					xlimits = (0, xlimits[1] + abs(xlimits[0]))
				if xlimits[1] > width:
					xlimits = (width - newWidth, width)
				ylimits = (int(y - newHeight/2), int(y + newHeight/2)) 
				if ylimits[0] < 0:
					ylimits = (0, ylimits[1] + abs(ylimits[0]))
				if ylimits[1] > height:
					ylimits = (height - newHeight, height)
				
			xlimits = (int(xlimits[0]), int(xlimits[1]))
			ylimits = (int(ylimits[0]), int(ylimits[1]))
			print "new limits:", xlimits, ylimits
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			print "new limits (world)", ra_limits, dec_limits
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			ppgplot.pggray(boostedImage, xlimits[0], xlimits[1]-1, ylimits[0], ylimits[1]-1, 0, 255, imagePlot['pgPlotTransform'])
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]
			plotCircles(dr2nearby, margins)

		if keyPressed=='o':
			if newWidth == width: continue
			print "Zoom out requested at (%0.0f, %0.0f)"%(x, y)
			zoomFactor = 0.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			
			if (newWidth >= width) or (newHeight>=height):
				print "Back to 1:1 scale."
				newWidth = width
				newHeight = height
				
			xlimits = (int(x - newWidth/2), int(x + newWidth/2)) 
			if xlimits[0] < 0:
				xlimits = (0, xlimits[1] + abs(xlimits[0]))
			if xlimits[1] > width:
				xlimits = (width - newWidth, width)
			ylimits = (int(y - newHeight/2), int(y + newHeight/2)) 
			if ylimits[0] < 0:
				ylimits = (0, ylimits[1] + abs(ylimits[0]))
			if ylimits[1] > height:
				ylimits = (height - newHeight, height)
				
			xlimits = (int(xlimits[0]), int(xlimits[1]))
			ylimits = (int(ylimits[0]), int(ylimits[1]))
			print "new limits:", xlimits, ylimits
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			print "new limits (world)", ra_limits, dec_limits
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			ppgplot.pggray(boostedImage, xlimits[0], xlimits[1]-1, ylimits[0], ylimits[1]-1, 0, 255, imagePlot['pgPlotTransform'])
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]
			plotCircles(dr2nearby, margins)
			
			
			
	# except KeyboardInterrupt:
	#	print "Ctrl-C pressed, but I dealt with it. "


	ppgplot.pgclos()
	
	
	