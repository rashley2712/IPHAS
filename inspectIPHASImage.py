#!/usr/bin/env python

import argparse, sys
import datetime, time
import numpy, math
import ppgplot
import generalUtils, configHelper, os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.table import Table, vstack
from astropy.utils import data
import astropy.table

def createPGplotWindow(handle, width, height):
	""" Set up the PGPLOT windows """
	
	newPlot = {}
	newPlot['pgplotHandle'] = ppgplot.pgopen('/xs')
	ppgplot.pgpap(2, 1)
	ppgplot.pgsvp(0.0, 1.0, 0.0, 1.0)
	ppgplot.pgswin(0, width, 0, height)
	
	return newPlot

def getBrightStars(ra, dec, radius):
	print(ra, dec, radius)
	maglimit = 30
	
	# First look for a cached copy of this data
	filenameParts = args.filename.split('.')
	usnoCache = filenameParts[0] + "_usno_cache.fits"
	usnoCached = False
	if not args.ignorecache:
		print("Looking for a cached copy of the USNO catalogue:", usnoCache)
		if os.path.exists(usnoCache):
			usnoCached = True
	
	if usnoCached:
		brightStarsTable = Table.read(usnoCache)
	else:		
		with data.conf.set_temp('remote_timeout', 60):
			try: 
				usno = 'The USNO-A2.0 Catalogue (Monet+ 1998) 1'
				search = conesearch(center=(ra, dec),
               		radius=radius,
                	verb=3,
					cache=True, 
                	catalog_db=usno)
			except: 
				print("Failed to retrieve any results from Vizier.")
				return None
			brightStarsTable = search.to_table()
			print("Found %d bright stars in %f degree radius."%(len(brightStarsTable), radius))
			brightStarsTable.write(usnoCache, format='fits', overwrite=True)
			
	brightStarsArray = []
	
	for row in brightStarsTable:
		if row['Rmag']<maglimit:
			star = {}
			star['ra'] = row['RAJ2000']
			star['dec'] = row['DEJ2000']
			star['Rmag'] = row['Rmag']
			x, y = wcsSolution.all_world2pix([star['ra']], [star['dec']], 1)
			star['x'] = x
			star['y'] = y
			brightStarsArray.append(star)
			if star['Rmag']>12:
				star['radius'] = 40*math.exp((-star['Rmag']+12)/4)
			else: 
				star['radius'] = 40
	print ("%d are brighter than R=%.1f"%(len(brightStarsArray), maglimit))
	return brightStarsArray

def distance(pos1, pos2):
	return math.sqrt( (pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2)

def gridCircle(x0, y0, radius, grid):
	x = radius;
	y = 0;
	decisionOver2 = 1 - x;   

	while( y <= x ):
		grid[ x0: x + x0,  y0: y + y0] = 132 # Octant 1
		grid[ x0: y + x0,  y0: x + y0] = 132 # Octant 2
		grid[ -x + x0:x0, y0:y + y0] = 132 # Octant 4
		grid[ -y + x0: x0, y0:x + y0] = 132 # Octant 3
		grid[-x + x0:x0, -y + y0:y0] = 132 # Octant 5
		grid[-y + x0:x0, -x + y0:y0] = 132 # Octant 6
		grid[ x0:x + x0, -y + y0:y0] = 132 # Octant 7
		grid[ x0:y + x0, -x + y0:y0] = 132 # Octant 8
		y+= 1
		if (decisionOver2<=0):
			decisionOver2 += 2 * y + 1;  # Change in decision criterion for y -> y+1
		else:
			x-= 1;
			decisionOver2 += 2 * (y - x) + 1;   # Change for y -> y+1, x -> x-1
	return grid

def getVizierObjects(ra, dec, radius):
	with data.conf.set_temp('remote_timeout', 60):
		try: 
			search = conesearch(center=(ra, dec),
                radius=0.05,
                verb=3,
                catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
		except: 
			print("Failed to retrieve any results from Vizier.")
	return search.to_table()

def plotCircles(objectTable, margins):
	margins = checkMargins(margins)
	ppgplot.pgsci(3)
	ppgplot.pgsfs(2)
	index = 0
	print("Margins:", margins)
	for obj in objectTable:
		ra = obj['ra']
		dec = obj['dec']
		c = obj['class']
		if ra > margins[1][0] and ra < margins[0][0] and dec<margins[0][1] and dec>margins[1][1]:
			# print index, ra, dec, x, y, c
			index+= 1
			colour = 1
			if c==-9: colour = 2   # Red  = Saturated
			if c==1: colour = 4    # Blue   = Galaxy
			if c==-3: colour = 5   # Cyan   = Probable Galaxy
			if c==-1: colour = 3   # Green = Star
			if c==-2: colour = 8   # Orange = Probable Star
			if c==0: colour = 2    # Red    = Noise
			ppgplot.pgsci(colour)
			ppgplot.pgcirc(obj['x'], obj['y'], 5 + (5* obj['pStar']))
	return (index+1)
	
def redraw():
	ppgplot.pgslct(imagePlot['pgplotHandle'])
	ppgplot.pgslw(3)
	ppgplot.pggray(boostedImage, xlimits[0], xlimits[1]-1, ylimits[0], ylimits[1]-1, imageMinMax[0], imageMinMax[1], imagePlot['pgPlotTransform'])
	if plotSources: plotCircles(dr2Objects, margins)
	if plotHa:
		reduceddr2cat = []
		for selected in extendedHaSources:
			reduceddr2cat.append(dr2Objects[selected])
		plotCircles(reduceddr2cat, margins)
	if plotGrid: 
		print("Plotting grid")
		ppgplot.pgsci(6)
		xVals = [p[0] for p in pixelGrid] 
		yVals = [p[1] for p in pixelGrid]
		ppgplot.pgpt(xVals, yVals, 2)
	if plotPointings:
		ppgplot.pgsci(2)
		ppgplot.pgsfs(2)
		ppgplot.pgslw(10)
		for p in pointings:
			ppgplot.pgcirc(p['x'], p['y'], 30)
	if plotBrightStars:
		ppgplot.pgsci(3)
		ppgplot.pgsfs(2)
		ppgplot.pgslw(10)
		for b in brightStars:
			ppgplot.pgcirc(b['x'], b['y'], 40)
			
			
def makeGrid(width, height, size):
	grid = []
	for x in numpy.arange(size, width, size):
		for y in numpy.arange(size, height, size):
			grid.append((x,y))
	return grid	
		
def makeMask():
	bitmap = numpy.zeros(numpy.shape(imageData))
	border = 50
	# First mask off the border
	bitmap[0:border, 0:width] = 132
	bitmap[height-border:height, 0:width] = 132
	bitmap[0:height, 0:border] = 132
	bitmap[0:height, width-border:width] = 132
	
	for index, object in enumerate(dr2Objects):
		if object['class'] != -1: continue   # Skip objects that are not stars  
		radius = object['pixelFWHM'] * 4
		x = object['x']  
		y = object['y'] 
		if (x<border) or (x>(width-border)): continue
		if (y<border) or (y>(height-border)): continue
		bitmap = gridCircle(y, x, radius, bitmap)
		sys.stdout.write("\rMasking: %d of %d."%(index, len(dr2Objects)))
		sys.stdout.flush()
	sys.stdout.write("\n")
	sys.stdout.flush()
	
	# Also mask out the really bright stars
	for index, object in enumerate(brightStars):
		radius = object['radius']
		x = object['x']  
		y = object['y'] 
		if (x<border) or (x>(width-border)): continue
		if (y<border) or (y>(height-border)): continue
		bitmap = gridCircle(y, x, radius, bitmap)
		
		
	return bitmap
	
def drawSuperPixel(superPixel):
	
	if 'pgplotHandle' not in maskPlot.keys():
		maskPlot['pgplotHandle'] = ppgplot.pgopen('/xs')
		maskPlot['pgPlotTransform'] = [0, 1, 0, 0, 0, 1]
	else:
		ppgplot.pgslct(maskPlot['pgplotHandle'])
	ppgplot.pgpap(paperSize, aspectRatio)
	ppgplot.pgsvp(0.0, 1.0, 0.0, 1.0)
	ppgplot.pgswin(0, width, 0, height)	
	ppgplot.pggray(mask, 0, width-1, 0, height-1, 0, 255, maskPlot['pgPlotTransform'])
	ppgplot.pgslct(imagePlot['pgplotHandle'])
	
def drawMask(mask):
	print("Drawing the mask.")
	if 'pgplotHandle' not in maskPlot.keys():
		maskPlot['pgplotHandle'] = ppgplot.pgopen('/xs')
		maskPlot['pgPlotTransform'] = [0, 1, 0, 0, 0, 1]
	else:
		ppgplot.pgslct(maskPlot['pgplotHandle'])
	ppgplot.pgpap(paperSize, aspectRatio)
	ppgplot.pgsvp(0.0, 1.0, 0.0, 1.0)
	ppgplot.pgswin(0, width, 0, height)	
	ppgplot.pggray(mask, 0, width-1, 0, height-1, 0, 255, maskPlot['pgPlotTransform'])
	ppgplot.pgslct(imagePlot['pgplotHandle'])
	
	

def filterGrid(pixelGrid):
	print("Old grid length:%d"%len(pixelGrid))
	newGrid = []
	radius = 7
	rejectXarray = [object['x'] for object in dr2Objects]
	rejectYarray = [object['y'] for object in dr2Objects]
	for position in pixelGrid:
		found = False
		for x,y in zip(rejectXarray, rejectYarray):
			if distance((x,y), position) < radius:
				print (position,  " is rejected, distance was %f"%distance((x,y), position))
				found=True
				continue
		if not found: newGrid.append(position)
	print("New grid length:%d"%len(newGrid))
	return newGrid
		
		
def checkMargins(margins):
	if margins[0][0] < margins[1][0]:
		temp = margins[0][0]
		margins[0][0] = margins[1][0]
		margins[1][0] = temp
	if margins[0][1] < margins[1][1]:
		temp = margins[0][1]
		margins[0][1] = margins[1][1]
		margins[1][1] = temp
	return margins
		
def withinMargins(table, namedColumns):
	print(namedColumns)
	return True
	

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Loads an IFAS reduced image. Displays it with PGPLOT.')
	parser.add_argument('filename', type=str, help='The FITS image file.')
	parser.add_argument('--save', action="store_true", help='Write the input parameters to the config file as default values.')
	parser.add_argument('--ignorecache', action="store_true", help='Ignore the cached DR2 catalogue. Will overwrite one if it already exists.')
	args = parser.parse_args()
	print(args)
	
	config = configHelper.configClass("inspectIPHASImage")
	if args.save:
		config.save()
	
	paperSize = 5.5  # Paper size in inches
	imageMinMax = (0, 255)
	plotSources = False
	plotHa = False
	plotGrid = False
	plotBrightStars = False
	brightStars = []
	invertedColours = True
	pixelGrid = []
	maskPlot = {}
	mask = []
	pixelScale = 0.333 # arcseconds per pixel
	pointings = []
	plotPointings = False
	
	print ("Astropy cache dir %s."%astropy.config.get_cache_dir())
		
	hdulist = fits.open(args.filename)
	
	print(hdulist.info())
	
	filter = None
	for card in hdulist:
		print(card.header.keys())
		print(repr(card.header))
		for key in card.header.keys():
			if 'WFFBAND' in key:
				filter = card.header[key]
	
	imageData =  hdulist[1].data
	savedImageData = numpy.copy(imageData)

	
	wcsSolution = WCS(hdulist[1].header)
	
	hdulist.close()
	
	(height, width) = numpy.shape(imageData)
	
	aspectRatio = float(height)/float(width)
	print(aspectRatio)
	
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
	print("RA, DEC of image centre is: ", positionString, ra, dec)
	margins = wcsSolution.all_pix2world([[0, 0], [width, height]], 1)
	margins = checkMargins(margins)
	print("ra, dec limits:", margins)
	
	
	print("Looking for bright stars")
	brightStars = getBrightStars(ra, dec, 0.5)
	
	filenameParts = args.filename.split('.')
	dr2Filename = filenameParts[0] + "_dr2_cache.fits"
	cached = False
	if not args.ignorecache:
		print("Looking for a cached copy of the DR2 catalogue:", dr2Filename)
		if os.path.exists(dr2Filename):
			cached = True

	if not cached:
		for index, yCentre in enumerate(numpy.arange(height/4, height, height/4)):
			print(yCentre)
			imageCentre = [ width/2, yCentre]
			ra, dec = wcsSolution.all_pix2world([imageCentre], 1)[0]
			print(index, imageCentre, ra, dec)
			with data.conf.set_temp('remote_timeout', 60):
				try: 
					search = conesearch(center=(ra, dec),
	                    radius=0.1,
	                    verb=3,
	                    catalog_db="http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&")
				except: 
					print("Failed to retrieve any results from Vizier.")
			dr2nearbyTemp = search.to_table()
			if index==0: 
				dr2nearby = dr2nearbyTemp
			else:
				print("Table was %d rows. Additional data is %d rows."%(len(dr2nearby), len(dr2nearbyTemp)))
				dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
				dr2nearby = dr2nearbyWhole
				print("New table is %d rows."%len(dr2nearby))
			print(dr2nearby)
		
		radius = 0.05	
		# Get bottom left corner from vizier
		ra, dec = wcsSolution.all_pix2world([[0, 0]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		# Get bottom right corner from vizier
		ra, dec = wcsSolution.all_pix2world([[width, 0]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		# Get top left corner from vizier
		ra, dec = wcsSolution.all_pix2world([[0, height]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		# Get top right corner from vizier
		ra, dec = wcsSolution.all_pix2world([[width, height]], 1)[0]
		print ("Getting corner: %f, %f"%(ra, dec))
		dr2nearbyTemp = getVizierObjects(ra, dec, radius) 
		dr2nearbyWhole = vstack([dr2nearby, dr2nearbyTemp])
		dr2nearby = dr2nearbyWhole
		
		print("Checking for duplicate objects.")
		IPHASnames = []
		for index, row in enumerate(dr2nearby):
			IPHASname = row['IPHAS2']
			if IPHASname in IPHASnames:
				dr2nearby.remove_row(index)
			else:
				IPHASnames.append(IPHASname)
			if (index%100) == 0:
				sys.stdout.write("\rChecking row: %d    "%index)
				sys.stdout.flush()
		
		sys.stdout.write("\n")
		sys.stdout.flush()
		print("dr2nearby table is %d rows long."%len(dr2nearby))
		print("Trimming out objects that lie outside the limits of this image.")
		margins = checkMargins(margins)
		print("Margins: " + str(margins))
		rejectedRows = 0
		for index, row in enumerate(dr2nearby):
			ra = row['RAJ2000']
			dec = row['DEJ2000']
			if ra>margins[0][0] or ra<margins[1][0] or dec>margins[0][1] or dec<margins[1][1]:
				sys.stdout.write("\rRejecting this row: %d %f, %f is outside image borders."%(index, ra, dec))
				sys.stdout.flush()
				dr2nearby.remove_row(index)
				rejectedRows+= 1
		
		sys.stdout.write("\n%d rows rejected."%rejectedRows)
		sys.stdout.flush()
		print("dr2nearby table is %d rows long."%len(dr2nearby))
				
		dr2nearby.write(dr2Filename, format='fits', overwrite=True)
	
	else:
		dr2nearby = Table.read(dr2Filename)
	
	print("Data columns found in the DR2 catalogue: " + str(dr2nearby.colnames))
	
	dr2nearby.remove_column('errBits2')
	
	print("Length of the dr2 table: %d"%len(dr2nearby))
	
	# Move table into a dictionary object
	dr2Objects = []
	for row in dr2nearby:
		IPHAS2name = row['IPHAS2']
		dr2Object={}
		dr2Object['name'] = IPHAS2name
		dr2Object['ra'] = row['RAJ2000']
		dr2Object['dec'] = row['DEJ2000']
		dr2Object['class'] = row['mergedClass']
		dr2Object['pStar'] = row['pStar']
		dr2Object['iClass'] = row['iClass']
		dr2Object['haClass'] = row['haClass']
		dr2Object['pixelFWHM'] = row['haSeeing'] / pixelScale
		x, y = wcsSolution.all_world2pix([dr2Object['ra']], [dr2Object['dec']], 1)
		dr2Object['x'] = x[0]
		dr2Object['y'] = y[0]
		dr2Objects.append(dr2Object)
		
	
	
	# Run through all objects and find the ones that are haClass = "+1" but iClass != "+1" and overall class = "+1"
	extendedHaSources = []
	for index, d in enumerate(dr2Objects):
		if (d['class'] ==1) and (d['haClass'] == 1) and (d['iClass'] != 1):
			extendedHaSources.append(index)
	
	print("Found %d extended Ha sources out of %d total objects in DR2."%(len(extendedHaSources), len(dr2Objects)))
	
	
	xlimits = (0, width)
	ylimits = (0, height)
	
	# plotCircles(dr2Objects, margins)
	redraw()

	help = []
	helpItem = {'key': "p", 'text': "Toggle the plotting of DR2 sources on/off."}
	help.append(helpItem)
	helpItem = {'key': "h", 'text': "Plot only the Ha sources."}
	help.append(helpItem)
	helpItem = {'key': "w", 'text': "Invert the gray-scale."}
	help.append(helpItem)
	helpItem = {'key': "i/o", 'text': "Zoom in/out."}
	help.append(helpItem)
	helpItem = {'key': "q", 'text': "Quit."}
	help.append(helpItem)
	helpItem = {'key': "f", 'text': "Find local Ha maxima."}
	help.append(helpItem)	
	helpItem = {'key': "l", 'text': "Show location of cursor."}
	help.append(helpItem)
	# try: 
	x=width/2
	y=height/2
	newWidth = width
	newHeight = height
	pixelGrid = makeGrid(width, height, 100)
	
	ppgplot.pgsci(3)
	keyPressed = None
	while keyPressed != 'q':
		ch = ppgplot.pgcurs(x, y)
		x=ch[0]
		y=ch[1]
		keyPressed = ch[2]
		# print "Key pressed:", ch[2]
		if keyPressed== '?':
			print("Help:")
			for helpItem in help:
				print("\t" + helpItem['key'] + " : " +  helpItem['text'])
			print("")
		if keyPressed== 'p':
			if plotSources == True:
				plotSources = False
			else:
				plotSources = True
				plotHa = False
			redraw()
		
		if keyPressed== 'h':
			if not plotHa:
				plotHa = True
				plotSources = False
			else:
				plotHa = False
			redraw() 
		
			
		if keyPressed== 'w':
			imageMinMax = (imageMinMax[1], imageMinMax[0])
			if invertedColours: invertedColours= False
			else: invertedColours=True
			redraw()
			
		if keyPressed== 'l':
			ra, dec = wcsSolution.all_pix2world(x, y, 1)
			positionString = generalUtils.toSexagesimal((ra, dec))
			print "Cursor location: [%d, %d]  %s (%f, %f)"%(x, y, positionString, ra, dec)
			
			
		if keyPressed=='i':
			print("Zoom requested at (%0.0f, %0.0f)"%(x, y))
			zoomFactor = 1.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			if (newWidth<100) or (newHeight<100):
				print("Maximum zoom reached.")
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
			print("new limits:", xlimits, ylimits)
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			print("new limits (world)", ra_limits, dec_limits)
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]
			redraw()
			
		if keyPressed=='g':
			if not plotGrid:
				plotGrid = True
				redraw()
			else:
				plotGrid = False
				redraw()
				
		if keyPressed=='m':
			mask = makeMask()
			drawMask(mask)
		
		if keyPressed=='c':
			boostedImage+= mask
			redraw()
	
		if keyPressed=='a':
			booleanMask = numpy.ma.make_mask(mask)
			maskedImageData = numpy.ma.masked_array(imageData, booleanMask)
			maskedBoostedImage = numpy.ma.masked_array(boostedImage, booleanMask)
			imageData = numpy.ma.filled(maskedImageData, 0)
			boostedImage = numpy.ma.filled(maskedBoostedImage, 0)
			redraw()
	
		if keyPressed=='f':
			radius = 180
			superPixelSize = 15
			superPixelArea = superPixelSize*superPixelSize
			spPreview = createPGplotWindow("preview", superPixelSize, superPixelSize)
			imageCopy = numpy.copy(imageData)
			maskedImageCopy = numpy.ma.masked_array(imageCopy, mask = False)
			for index in range(50):
				maximum = numpy.amax(imageCopy)
				flatPosition =  numpy.argmax(imageCopy)
				position = numpy.unravel_index(flatPosition, numpy.shape(imageData))
				print(maximum, position[1], position[0])
				startX = position[1]-superPixelSize/2
				startY = position[0]-superPixelSize/2
				superPixel = imageCopy[startY:startY+superPixelSize, startX:startX+superPixelSize]
				ppgplot.pgslct(spPreview['pgplotHandle'])
				boostedImage = generalUtils.percentiles(superPixel, 20, 99)
				ppgplot.pggray(boostedImage, 0, superPixelSize-1, 0, superPixelSize-1, 0, 255, imagePlot['pgPlotTransform'])
				ppgplot.pgslct(imagePlot['pgplotHandle'])
				pointingObject = { 'x': position[1], 'y': position[0]}
				pointings.append(pointingObject)
				ppgplot.pgsfs(2)
				ppgplot.pgsci(5)
				ppgplot.pgcirc(position[1], position[0], radius)
				xpts = [startX, startX, startX+superPixelSize, startX+superPixelSize]
				ypts = [startY, startY+superPixelSize, startY+superPixelSize, startY]
				ppgplot.pgsfs(2)
				ppgplot.pgsci(4)
				ppgplot.pgpoly(xpts, ypts)
				tempBitmap = numpy.zeros(numpy.shape(imageCopy))
				tempBitmap = gridCircle(position[0], position[1], radius, tempBitmap)
				booleanMask = numpy.ma.make_mask(tempBitmap)
				maskedImageCopy = numpy.ma.masked_array(imageCopy, booleanMask)
				imageCopy = numpy.ma.filled(maskedImageCopy, 0)
				time.sleep(3)
				
			print pointings	
				
				
		if keyPressed=='z':
			plotPointings = True
			imageData = savedImageData
			boostedImage = generalUtils.percentiles(imageData, 20, 99)
			redraw()
		
		if keyPressed=='b':
			if plotBrightStars:
				plotBrightStars = False
			else:
				plotBrightStars = True
			redraw()
			
		if keyPressed=='r':
			pixelGrid = filterGrid(pixelGrid)
			redraw() 
			
		if keyPressed=='o':
			if newWidth == width: continue
			print("Zoom out requested at (%0.0f, %0.0f)"%(x, y))
			zoomFactor = 0.5
			currentWidth = (xlimits[1] - xlimits[0])
			newWidth = currentWidth / zoomFactor
			currentHeight = (ylimits[1] - ylimits[0])
			newHeight = currentHeight / zoomFactor
			
			if (newWidth >= width) or (newHeight>=height):
				print("Back to 1:1 scale.")
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
			print("new limits:", xlimits, ylimits)
			ra_limits, dec_limits = wcsSolution.all_pix2world(numpy.array(xlimits), numpy.array(ylimits), 1)
			margins = [[ra_limits[0], dec_limits[0]], [ra_limits[1], dec_limits[1]]]	
			print("new limits (world)", ra_limits, dec_limits)
			
			ppgplot.pgswin(xlimits[0], xlimits[1], ylimits[0], ylimits[1])
			redraw()
		
		if plotSources: sourceStatus = "ON"
		else: sourceStatus = "OFF" 
		if plotHa: HaStatus = "ON"
		else: HaStatus = "OFF" 
		if invertedColours: invertStatus = "ON"
		else: invertStatus = "OFF" 
		print("Filter: %s PlotSources[%s]  Plot Ha extended[%s]  Invert Greyscale[%s]"%(filter, sourceStatus, HaStatus, invertStatus))
			
			
	# except KeyboardInterrupt:
	#	print "Ctrl-C pressed, but I dealt with it. "


	ppgplot.pgclos()
	
	
	
