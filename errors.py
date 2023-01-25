class NotPairedReads(Exception):
	"""Exception raised for error where paired reads are not found 
		using the nomenclautre of the files."""

	def __init__(self, ID, message):
		self.ID = ID
		self.message = message
		super().__init__(self.message)

#ID = 2342323
#errmsg= f"Paired reads for FastQ files with ID {ID} have not been identified"
#raise NotPairedReads(ID,errmsg)

class moreThanTwoReads(Exception):
	"""Exception raised when an ID has more than 2 reads"""

	def __init__(self,ID, message):
		self.ID = ID
		self.message = message
		super().__init__(self.message)

class TruncatedBamError(Exception):
	"""Exception raised when a truncated BAM file is found """

	def __init__(self, bam_file, message):
		self.bam_file = bam_file
		self.message = message
		super().__init__(self.message)

class NotTypicalIDNomenclature(Exception):
	"""Raised when nomenclature"""

	def __init__(self, bam_file, message):
		self.bam_file = bam_file
		self.message = message
		super().__init__(self.message)

