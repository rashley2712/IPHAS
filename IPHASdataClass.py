from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.vo.client.conesearch import list_catalogs
from astropy.table import Table, vstack
from astropy.utils import data
from matplotlib.path import Path

import numpy, math, os, sys
import generalUtils
import astroquery
import matplotlib.pyplot

def distance(p1, p2):
	return math.sqrt( (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 )

catalogMetadata = {
	'tycho': {
		'columns': {
			'ra': 'RA_ICRS_',
			'dec': 'DE_ICRS_',
			'B': 'BTmag',
			'V': 'VTmag',
			'mag': 'VTmag' },
		'catalog_db': "http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=I/259/tyc2&-out.all&",
		'colour': 'blue',
		'VizierName': 'I/259/tyc2',
		'VizierLookup': 'tycho'
		}, 
	'usno': {
		'columns': {
					'ra': 'RAJ2000',
					'dec': 'DEJ2000',
					'B': 'Bmag',
					'R': 'Rmag',
					'mag': 'Rmag' },
				'colour': 'blue',
				'VizierName': 'I/252/out',
				'VizierLookup': 'usno'
			},
	'dr2': {
		'columns': {
			'ra': 'RAJ2000', 
			'dec': 'DEJ2000',
			'i': 'i',
			'r': 'r',
			'H': 'ha',
			'mag': 'r',
		    'class': 'mergedClass',
			'pStar': 'pStar', 
		    'iclass': 'iClass', 
		    'haClass': 'haClass',
		    'pixelFWHM': 'haSeeing'},
		'catalog_db': "http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&",
		'VizierLookup': 'DR2',
		'VizierName': 'II/321/iphas2',
		'colour': 'green'
	}
}

class IPHASdataClass:
	def __init__(self):
		print "Initialising an empty IPHAS data class"
		self.originalImageData = None
		self.boostedImage = None
		self.FITSHeaders = {}
		self.filter = None
		self.pixelScale = None
		self.centre = None
		self.filename = None
		self.ignorecache = False
		self.catalogs = {}
		self.figSize = 12.
		self.magLimit = 18
		self.mask = None
		return None
		
	def setProperty(self, property, value):
		truths = ["true", "yes", "on", "1", "Y", "y", "True"]
		falses = ["false", "no", "off", "0", "N", "n", "False"]
		if property=='maglimit':
			self.__dict__['magLimit'] = float(value)
		if property=="ignorecache":
			if value in truths:
				self.ignorecache = True
			if value in falses:
				self.ignorecache = False
			
			
		
	def loadFITSFile(self, filename):
		hdulist = fits.open(filename)
		self.filename = filename
		FITSHeaders = []
		for card in hdulist:
			# print(card.header.keys())
			# print(repr(card.header))
			for key in card.header.keys():
				self.FITSHeaders[key] = card.header[key]
				if 'WFFBAND' in key:
					self.filter = card.header[key]
		self.originalImageData =  hdulist[1].data
		self.height, self.width = numpy.shape(self.originalImageData)
		self.wcsSolution = WCS(hdulist[1].header)
		print "width, height", self.width, self.height, "shape:", numpy.shape(self.originalImageData)
		self.getRADECmargins()
		imageCentre = (self.width/2, self.height/2)
		ra, dec = self.wcsSolution.all_pix2world([imageCentre], 1)[0]
		self.centre = (ra, dec)
		positionString = generalUtils.toSexagesimal((ra, dec))
		print "RA, DEC of image centre is: ", positionString, ra, dec
		
		hdulist.close()
		
	def showVizierCatalogs(self):
		(ra, dec) = self.centre
		from astroquery.vizier import Vizier
		Vizier.ROW_LIMIT = 50
		from astropy import coordinates
		from astropy import units as u
		c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
		skyHeight= coordinates.Angle(self.raRange, unit = u.deg)
		results = Vizier.query_region(coordinates = c, radius= 1.0 * u.deg)
		print results
		
		
	def getVizierObjects(self, catalogName):
		""" Make a request to Vizier to get an Astropy Table of catalog object for this field. """
		(ra, dec) = self.centre
		
		availableCatalogs = catalogMetadata.keys()
		if catalogName not in availableCatalogs:
			print "The definitions for this catalogue are unknown. Available catalogues are:", availableCatalogs
			return
		
		# First look for a cached copy of this data
		filenameParts = self.filename.split('.')
		catalogCache = filenameParts[0] + "_" + catalogName + "_cache.fits"
		cached = False
		if not self.ignorecache:
			print "Looking for a cached copy of the catalogue:", catalogCache, 
			if os.path.exists(catalogCache):
				print "FOUND"
				cached = True
			else: print "NOT FOUND"
	
		if cached:
			newCatalog = Table.read(catalogCache)
		else:			
			print "Going online to fetch %s results from Vizier with mag limit %f."%(catalogName, self.magLimit)
			from astroquery.vizier import Vizier
			Vizier.ROW_LIMIT = 1E5
			Vizier.column_filters={"r":"<%d"%self.magLimit}
			from astropy import coordinates
			from astropy import units as u
			c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
			skyHeight= coordinates.Angle(self.raRange, unit = u.deg)
			skyWidth = coordinates.Angle(self.decRange, unit = u.deg)
			print "Sky width, height:", skyWidth, skyHeight
			result = Vizier.query_region(coordinates = c, width = skyWidth, height = skyHeight, catalog = catalogMetadata[catalogName]['VizierLookup'], verbose=True)
			newCatalog = result[catalogMetadata[catalogName]['VizierName']]
			newCatalog.pprint()
			
			# Write the new catalog to the cache file
			newCatalog.write(catalogCache, format='fits', overwrite=True)
		
		self.addCatalog(newCatalog, catalogName)
		
		return
		
		
	def printCatalog(self, catalogName):
		catalog = self.catalogs[catalogName]
		for b in catalog:
			print b
		print "%d rows printed."%len(catalog)
	
	def addCatalog(self, catTable, catalogName):
		newCatalog = []
		columnMapper = catalogMetadata[catalogName]['columns']
		for index, row in enumerate(catTable):
			object={}
			skipRow = False
			for key in columnMapper.keys():
				object[key] = row[columnMapper[key]]
				if numpy.isnan(row[columnMapper[key]]): skipRow = True
			if skipRow: continue		
			x, y = self.wcsSolution.all_world2pix([object['ra']], [object['dec']], 1)
			object['x'] = x[0]
			object['y'] = y[0]
			newCatalog.append(object)
			if  ((index+1)%100) == 0:
				sys.stdout.write("\rCopying: %d of %d."%(index+1, len(catTable)))
				sys.stdout.flush()
		sys.stdout.write("\rCopying: %d of %d.\n"%(index+1, len(catTable)))
		sys.stdout.flush()
		
		print "Adding catalog %s to list of stored catalogs."%catalogName
		self.catalogs[catalogName] =  newCatalog
		return
				
	def getRADECmargins(self):
		boundingBox = self.wcsSolution.all_pix2world([[0, 0], [0, self.width], [self.height, self.width], [self.height, 0]], 1, ra_dec_order = True)
		# boundingBox = self.wcsSolution.all_pix2world([[0, 0], [0, self.height], [self.width, self.height], [self.width, 0]], 1, ra_dec_order = True)
		print "Bounding box:", boundingBox
		pixelDiagonal = math.sqrt(self.height**2 + self.width**2)
		pixel1 = boundingBox[0]
		pixel2 = boundingBox[2]
		skyDiagonal = distance(pixel1, pixel2)
		print "Diagonal size:", pixelDiagonal, skyDiagonal
		self.pixelScale = (skyDiagonal / pixelDiagonal) * 3600.
		print "Pixel scale: %3.2f \"/pixel"%self.pixelScale
		raMin = numpy.min([r[0] for r in boundingBox])
		raMax = numpy.max([r[0] for r in boundingBox])
		decMin = numpy.min([r[1] for r in boundingBox])
		decMax = numpy.max([r[1] for r in boundingBox])
		print "RA, DEC min/max:", raMin, raMax, decMin, decMax
		raRange = raMax - raMin
		decRange = decMax - decMin
		print "RA range, DEC range", raRange, decRange, raRange*60, decRange*60
		self.raRange = raRange
		self.decRange = decRange
		print "Pixel scale: %6.4f \"/pixel"%self.pixelScale
		self.boundingBox = boundingBox
		
	def showFITSHeaders(self):
		headersString = ""
		for key in self.FITSHeaders.keys():
			print key + " : " + str(self.FITSHeaders[key])
			headersString+= str(key) + " : " + str(self.FITSHeaders[key]) + "\n"
		return headersString
			
	def getFITSHeader(self, key):
		try:
			print key + " : " + str(self.FITSHeaders[key]) 
			return self.FITSHeaders[key]
		except KeyError:
			print "Could not find a header with the name:", key
			return None 
			
	def plotCatalog(self, catalogName):
		catalog = self.catalogs[catalogName]
		catalogColour = catalogMetadata[catalogName]['colour']
		try:
			fig = self.figure
			
			xArray = []
			yArray = []
			rArray = []
			for o in catalog:
				# Check that the catalog has a class flag
				if 'class' in o.keys():
					if o['class'] != -1: continue   # Skip objects that are not stars  
				xArray.append(o['x'])
				yArray.append(self.height - o['y'])
				if catalogName=='dr2':
					r = o['pixelFWHM']*5.
				else:
					if o['mag']>12:
						r = 40*math.exp((-o['mag']+12)/4)
					else: 
						r = 40
				rArray.append(r)
	
			# Nick Wright 
			# R / pixels = 8192/M^2 + 1000/M + 100 
			
			patches = [matplotlib.pyplot.Circle((x_, y_), s_, fill=False, linewidth=1) for x_, y_, s_ in numpy.broadcast(xArray, yArray, rArray)]
			collection = matplotlib.collections.PatchCollection(patches, alpha = 0.25, color = catalogColour)
			ax = matplotlib.pyplot.gca()
			ax.add_collection(collection)
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			# matplotlib.pyplot.savefig("test.png", bbox_inches='tight')
		except AttributeError as e:
			print "There is no drawing surface defined yet. Please use the 'draw' command first."
			print e
		except Exception as e:
			print e
			
	def maskCatalog(self, catalogName):
		try:
			catalog = self.catalogs[catalogName]
		except KeyError:
			print "Could not find a catalog called %s."%catalogName
			return
		
		if self.mask==None:
			self.mask = numpy.zeros(numpy.shape(self.originalImageData))
			print "Creating a new blank mask of size:", numpy.shape(self.mask)

		xArray = []
		yArray = []
		rArray = []
		for o in catalog:
			# Check that the catalog has a class flag
			if 'class' in o.keys():
				if o['class'] != -1: continue   # Skip objects that are not stars  
			xArray.append(o['x'])
			yArray.append(self.height - o['y'])
			if catalogName=='dr2':
				r = o['pixelFWHM']*5.
			else:
				if o['mag']>12:
					r = 40*math.exp((-o['mag']+12)/4)
				else: 
					r = 40
			rArray.append(r)
			
		#for x, y, r in zip(xArray, yArray, rArray):
		#	print x, y, r
		
			
	def drawBitmap(self):
		if self.boostedImage is None:
			print "Boosting the image"
			self.boostedImage = generalUtils.percentiles(self.originalImageData, 20, 99)
		matplotlib.pyplot.ion()
		# mplFrame = numpy.rot90(self.boostedImage)
		mplFrame = self.boostedImage
		mplFrame = numpy.flipud(mplFrame)
		self.figure = matplotlib.pyplot.figure(self.filename, figsize=(self.figSize/1.618, self.figSize))
		self.figure.frameon = False
		self.figure.set_tight_layout(True)
		axes = matplotlib.pyplot.gca()
		axes.set_axis_off()
		self.figure.add_axes(axes)
		imgplot = matplotlib.pyplot.imshow(mplFrame, cmap="gray_r", interpolation='nearest')
		
		verts = []
		for b in self.boundingBox:
			print b
			y, x = self.wcsSolution.all_world2pix(b[0], b[1], 1, ra_dec_order=True)
			coord  = (float(x), float(y))
			print coord
			verts.append(coord)	
			
		verts.append((0, 0))
			
		print verts
		codes = [Path.MOVETO,
		         Path.LINETO,
		         Path.LINETO,
		         Path.LINETO,
		         Path.CLOSEPOLY,
		         ]

		path = Path(verts, codes)

		patch = matplotlib.patches.PathPatch(path, fill=None, lw=2)
		axes.add_patch(patch)
		matplotlib.pyplot.show()
		
		# matplotlib.pyplot.savefig("test.png",bbox_inches='tight')
		
			
			
