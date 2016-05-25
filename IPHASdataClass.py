from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.vo.client.conesearch import list_catalogs
from astropy.table import Table, vstack
from astropy.utils import data

import numpy, math, os, sys
import generalUtils
import astroquery
import matplotlib.pyplot


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
	'dr2': {
		'columns': {
			'ra': 'RAJ2000', 
			'dec': 'DEJ2000',
			'i': 'i',
			'r': 'r',
			'H': 'ha',
			'mag': 'ha',
		    'class': 'mergedClass',
			'pStar': 'pStar', 
		    'iclass': 'iClass', 
		    'haClass': 'haClass',
		    'pixelFWHM': 'haSeeing'},
		'catalog_db': "http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=IPHAS2&-out.all&",
		'VizierLookup': 'DR2',
		'VizierName': 'II/321/iphas2',
		'color': 'green'
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
		self.catalog = []
		self.figSize = 12.
		self.magLimit = 18
		return None
		
	def setProperty(self, property, value):
		if property=='magLimit':
			self.__dict__[property] = float(value)
		if property==''
		
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
		print "width, height", self.width, self.height
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
			print "Going online to fetch %s results from Vizier."%catalogName
			from astroquery.vizier import Vizier
			Vizier.ROW_LIMIT = 1E5
			Vizier.column_filters={"r":"<%d"%self.magLimit}
			from astropy import coordinates
			from astropy import units as u
			c = coordinates.SkyCoord(ra,dec,unit=('deg','deg'),frame='icrs')
			skyHeight= coordinates.Angle(self.raRange, unit = u.deg)
			skyWidth = coordinates.Angle(self.decRange, unit = u.deg)
			print "Sky width, height:", skyWidth, skyHeight
			result = Vizier.query_region(coordinates = c, width = skyWidth, height = skyHeight, catalog = catalogMetadata[catalogName]['VizierLookup'])
			newCatalog = result[catalogMetadata[catalogName]['VizierName']]
			newCatalog.pprint()
			newCatalog.write(catalogCache, format='fits', overwrite=True)
		
		self.addCatalog(newCatalog, catalogName)
		
		return
		
		
	def printCatalog(self):
		for b in self.catalog:
			print b
		print "%d rows printed."%len(self.catalog)
	
	def addCatalog(self, catTable, catalogName):
		newCatalog = []
		columnMapper = catalogMetadata[catalogName]['columns']
		for index, row in enumerate(catTable):
			object={}
			for key in columnMapper.keys():
				object[key] = row[columnMapper[key]]
			x, y = self.wcsSolution.all_world2pix([object['ra']], [object['dec']], 1)
			object['x'] = x[0]
			object['y'] = y[0]
			newCatalog.append(object)
			if  (index%100) == 0:
				sys.stdout.write("\rCopying: %d of %d."%(index+1, len(catTable)))
				sys.stdout.flush()
		sys.stdout.write("\rCopying: %d of %d.\n"%(index+1, len(catTable)))
		sys.stdout.flush()
		self.catalog = newCatalog
				
	def addIPHASCatalog(self, catTable):
		for index, row in enumerate(catTable):
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
			if  (index%100) == 0:
				sys.stdout.write("\rCopying: %d of %d."%(index, len(dr2nearby)))
				sys.stdout.flush()
		sys.stdout.write("\n")
		sys.stdout.flush()
				
	def getRADECmargins(self):
		margins = self.wcsSolution.all_pix2world([[0, 0], [self.width, self.height]], 1)
		print "Margins:", margins
		if margins[0][0] > margins[1][0]:
			temp = margins[0][0]
			margins[0][0] = margins[1][0]
			margins[1][0] = temp
		if margins[0][1] > margins[1][1]:
			temp = margins[0][1]
			margins[0][1] = margins[1][1]
			margins[1][1] = temp
		print "Repaired margins:", margins
		raRange = abs(margins[0][0] - margins[1][0])
		decRange = abs(margins[0][1] - margins[1][1])
		print "RA range, DEC range", raRange, decRange, raRange*60, decRange*60
		longRange = raRange
		if decRange>longRange: longRange = decRange
		self.pixelScale = (longRange / self.height) * 3600.
		self.raRange = raRange
		self.decRange = decRange
		print "Pixel scale: %6.4f \"/pixel"%self.pixelScale
		
		return margins
		
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
			
	def plotCatalog(self):
		try:
			fig = self.figure
			for index, object in enumerate(self.catalog):
				x = object['x'] 
				y = self.height - object['y'] 
				radius = 10
				fig.gca().add_artist(matplotlib.pyplot.Circle((x,y), radius, color='green', fill=False, linewidth=1.0))
				if  (index%100) == 0:
						sys.stdout.write("\rPlotting: %d of %d."%(index+1, len(self.catalog)))
						sys.stdout.flush()
			sys.stdout.write("\rPlotting: %d of %d.\n"%(index+1, len(self.catalog)))
			sys.stdout.flush()
			matplotlib.pyplot.show()
			matplotlib.pyplot.savefig("test.png", bbox_inches='tight')
		except AttributeError:
			print "There is no drawing surface defined yet. Please use the 'draw' command first."
		except Exception as e:
			print e
			
			
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
		matplotlib.pyplot.savefig("test.png",bbox_inches='tight')
		
			
			
