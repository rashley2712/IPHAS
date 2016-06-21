import cmd, sys
import IPHASdataClass
import os

class commandClass(cmd.Cmd):
	""" The command processor """
	
	
	"""Simple command processor example."""
	prompt = 'iphas: '
	# Disable rawinput module use
	use_rawinput = False
	# Do not show a prompt after each command read
	prompt = ''
	
	IPHASdata = IPHASdataClass.IPHASdataClass()
	echo = False
	
	def default(self, line):
		command, arg, line = self.parseline(line)
		func = [getattr(self, n) for n in self.get_names() if n.startswith('do_' + command)]
		if func:
			func[0](arg)
		return 
	
	def do_versions(self, line):
		import sys
		print "Sys: %s"%sys.version
		import astropy
		print "Astropy: %s"%astropy.__version__
		import numpy
		print "Numpy: %s"%numpy.__version__
		import matplotlib
		print "Matplotlib: %s"%matplotlib.__version__
		from PIL import Image,ImageDraw,ImageFont
		print "Python Image Library (PIL): %s"%Image.VERSION
		import astroquery
		print "Astroquery: %s"%astroquery.__version__
		
		return
	
	def precmd(self, line):
		if len(line)==0: return line
		if line[0] == '#':
			# This line is a comment, therefore send a 'NOP' to the cmd processor to ignore.
			print "Comment:", line
			return "NOP"
			
		if line[0] == "@":
			# A script needs to be run.
			scriptname = line[1:]
			print "Running the script:", scriptname
			self.run_script(scriptname)
			return "NOP"
			
		if self.echo: print "iphas> ", line
		return line
		
	def do_echo(self, line):
		""" Toggle the echo of commands (useful if you are running from a script)."""
		if self.echo: self.echo=False
		else: self.echo=True
		return 
	
	def run_script(self, scriptname):
		try:
			scriptFile = open(scriptname, 'rt')
			for line in scriptFile:
				self.onecmd(line)
			return "NOP"
		except IOError as e:
			print "Could not find the script:", scriptname
			return 
			
	def do_reload(self, line):
		if 'IPHASdataClass' in sys.modules:  
		    del sys.modules["IPHASdataClass"]
		del IPHASdata
		import IPHASdataClass
		IPHASdata = IPHASdataClass.IPHASdataClass()
		return 
	
	def do_quit(self, line):
		""" quit 
		Leave and exit to the shell. """
		print "Leaving iphas. Goodbye."
		sys.exit()
		return True
		
	def do_type(self, line):
		""" type
		Print out to the terminal the contents of an object.
		The object could be a catalog downloaded from Vizier, etc. """
		
		self.IPHASdata.printCatalog(line)
			
		
	def do_plot(self, line):
		""" plot
		Plot a catalog over an existing image """
		
		self.IPHASdata.plotCatalog(line)
		
		return
		
	def do_draw(self, line):
		""" draw
		Draw the image of the CCD """
		
		self.IPHASdata.drawBitmap()
		
		return
		
	def do_apply(self, line):
		"""
		Apply mask to the bitmap """
		
		self.IPHASdata.applyMask()
		
		return
		
		
	def do_load(self, line):
		""" load
		Load a FITS file. """
		print "Loading...", line
		if not os.path.exists(line):
			print "Could not find a file named:", line
			return
		
		self.IPHASdata.loadFITSFile(line)
		return 
		
	def do_show(self, line):
		""" show
		Output some information about a stored object """
		if line=="headers":
			self.IPHASdata.showFITSHeaders()
			return
		if line=="margins":
			print "Margins:", self.IPHASdata.getRADECmargins()
			return
		if line=="catalogs":
			self.IPHASdata.showVizierCatalogs()
			return
		self.IPHASdata.getFITSHeader(line)
		return
	
	def do_shell(self, line):
		"Run a shell command"
		print "running shell command:", line
		output = os.popen(line).read()
		print output
		self.last_output = output
		
	def do_get(self, line):
		""" Get an additional object for the image. 
		eg get cat : Get a catalogue of sources from Vizier."""
		catalogName = line.split()[0]
		print "Getting ", catalogName
		self.IPHASdata.getVizierObjects(catalogName)
		
		return
		
	def do_set(self, line):
		""" set
		Set a property of the IPHAS object."""
		# print self.IPHASdata.__dict__
		params = line.split()
		if len(params)< 2:
			print "Usage 'set [property] [value]'"
			return
		propertyToSet = params[0]
		value = params[1]
		self.IPHASdata.setProperty(propertyToSet, value)
		return
		
	def do_mask(self, line):
		self.IPHASdata.maskCatalog(line)
		return
		
	def do_cats(self, line):
		self.IPHASdata.showVizierCatalogs()
		return
		
	def emptyline(self):
		return
		
	def do_EOF(self, line):
		return True
	
	def postloop(self):
		return True
	
	def do_NOP(self, line):
		""" NOP
		Do nothing."""
		return 
	
