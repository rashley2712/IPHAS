from astropy.io import fits
import numpy

class IPHASdataClass:
	def __init__(self):
		print "Initiliasing an empty data class"
		self.originalImageData = None
		self.FITSHeaders = {}
		self.filter = None
		return None
		
	def loadFITSFile(self, filename):
		hdulist = fits.open(filename)
	
		FITSHeaders = []
		for card in hdulist:
			# print(card.header.keys())
			# print(repr(card.header))
			for key in card.header.keys():
				self.FITSHeaders[key] = card.header[key]
				if 'WFFBAND' in key:
					self.filter = card.header[key]
		self.originalImageData =  hdulist[1].data
		self.width, self.height = numpy.shape(self.originalImageData)
		print "width, height", self.width, self.height
		
	def showFITSHeaders(self):
		headersString = ""
		for key in self.FITSHeaders.keys():
			print key + " : " + str(self.FITSHeaders[key])
			headersString+= str(key) + " : " + str(self.FITSHeaders[key]) + "\n"
		print "... headers 6"
			
	def getFITSHeader(self, key):
		try:
			print key + " : " + str(self.FITSHeaders[key]) 
			return self.FITSHeaders[key]
		except KeyError:
			print "Could not find a header with the name:", key
			return None 
			
			