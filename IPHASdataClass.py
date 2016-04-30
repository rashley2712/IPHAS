from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.table import Table, vstack
from astropy.utils import data

import numpy, math, os, sys
import generalUtils


class IPHASdataClass:
	def __init__(self):
		print "Initiliasing an empty data class"
		self.originalImageData = None
		self.FITSHeaders = {}
		self.filter = None
		self.pixelScale = None
		self.centre = None
		self.filename = None
		self.ignorecache = False
		return None
		
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
		imageCentre = (self.width/2, self.height/2)
		ra, dec = self.wcsSolution.all_pix2world([imageCentre], 1)[0]
		self.centre = (ra, dec)
		positionString = generalUtils.toSexagesimal((ra, dec))
		print "RA, DEC of image centre is: ", positionString, ra, dec
		
		hdulist.close()
		
	def getVizierObjects(self):
		""" Make a request to Vizier to get an Astropy Table of catalog object for this field. """
		fullRadius = math.sqrt((self.width/2)**2 + (self.height/2)**2) * self.pixelScale
		(ra, dec) = self.centre
		radius = fullRadius / 3600.
		print(ra, dec, radius)
		maglimit = 30
	
		# First look for a cached copy of this data
		filenameParts = self.filename.split('.')
		usnoCache = filenameParts[0] + "_usno_cache.fits"
		usnoCached = False
		if not self.ignorecache:
			print("Looking for a cached copy of the USNO catalogue:", usnoCache)
			if os.path.exists(usnoCache):
				usnoCached = True
	
		print conesearch.list_catalogs()
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
				
		self.addBrightCatalog(brightStarsTable)
	
	def addBrightCatalog(self, catTable):
		newCatalog = []
		print "Available columns:", catTable.columns
		for index, row in enumerate(catTable):
			object={}
			#object['name'] = IPHAS2name
			object['ra'] = row['RAJ2000']
			object['dec'] = row['DEJ2000']
			object['R'] = row['Rmag']
			object['B'] = row['Bmag']
			#dr2Object['class'] = row['mergedClass']
			#dr2Object['pStar'] = row['pStar']
			#dr2Object['iClass'] = row['iClass']
			#dr2Object['haClass'] = row['haClass']
			#dr2Object['pixelFWHM'] = row['haSeeing'] / pixelScale
			x, y = self.wcsSolution.all_world2pix([object['ra']], [object['dec']], 1)
			object['x'] = x[0]
			object['y'] = y[0]
			newCatalog.append(object)
			if  (index%100) == 0:
				sys.stdout.write("\rCopying: %d of %d."%(index+1, len(catTable)))
				sys.stdout.flush()
		sys.stdout.write("\rCopying: %d of %d.\n"%(index+1, len(catTable)))
		sys.stdout.flush()
				
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
		if margins[0][0] < margins[1][0]:
			temp = margins[0][0]
			margins[0][0] = margins[1][0]
			margins[1][0] = temp
		if margins[0][1] < margins[1][1]:
			temp = margins[0][1]
			margins[0][1] = margins[1][1]
			margins[1][1] = temp
		self.pixelScale = (abs(margins[0][0] - margins[1][0]) / self.height) * 3600.
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
			
			