import cmd, sys
import IPHASdataClass


class commandClass(cmd.Cmd):
	""" The command processor """
	IPHASdata = IPHASdataClass.IPHASdataClass()
	
	"""Simple command processor example."""
	prompt = 'iphas: '
	# Disable rawinput module use
	use_rawinput = False
	# Do not show a prompt after each command read
	prompt = ''
	
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
		return line
	
	def run_script(self, scriptname):
		try:
			scriptFile = open(scriptname, 'rt')
			for line in scriptFile:
				self.onecmd(line)
		except:
			print "Could not find the script:", scriptname
			return "NOP"
			
	def do_reload(self, line):
		if 'IPHASdataClass' in sys.modules:  
		    del sys.modules["IPHASdataClass"]
		import IPHASdataClass
		IPHASdata = IPHASdataClass.IPHASdataClass()
		return 
	
	def do_quit(self, line):
		""" quit 
		Leave and exit to the shell. """
		print "Leaving iphas. Goodbye."
		sys.exit()
		return True
		
	def do_load(self, line):
		""" load
		Load a FITS file. """
		print "Loading...", line
		
		self.IPHASdata.loadFITSFile(line)
		return 
		
	def do_show(self, line):
		""" show
		Output some information about an stored object """
		if line=="headers":
			self.IPHASdata.showFITSHeaders()
			return
		self.IPHASdata.getFITSHeader(line)
	
	def do_shell(self, line):
		"Run a shell command"
		print "running shell command:", line
		output = os.popen(line).read()
		print output
		self.last_output = output
		
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
	