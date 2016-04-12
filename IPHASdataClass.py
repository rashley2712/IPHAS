from astropy.io import fits
from astropy.wcs import WCS
from astropy.vo.client.conesearch import conesearch
from astropy.table import Table, vstack
from astropy.utils import data

import numpy, math, os
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
		fullRadius = math.sqrt((self.width/2)**2 + (self.height/2)**2) * self.pixelScale
		(ra, dec) = self.centre
		radius = fullRadius
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
			
			