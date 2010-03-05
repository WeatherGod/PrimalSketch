


class ScaleSpace_Blob :
	def __init__(self) :
		self.isBright = None
		self.signifVal = None
		self.appearance = None
		self.disappearance = None
		self.grey_blobs = []
		self.approp_scalelvl = None
		self.approp_greyblob = None
		self.atBoundary = None


class Bifurcation_Event :
	def __init__(self, event_type = None, scaleblobs_above = [], scaleblobs_below = [], position = None, scaleVal = None) :
		self.event_type = event_type
		self.scaleblobs_above = scaleblobs_above
		self.scaleblobs_below = scaleblobs_below
		self.position = position
		self.scaleVal = scaleVal


class Scale_Level :
	def __init__(self, scaleVal, image, coarserLevel = None, finerLevel = None) :
		self.scaleVal = scaleVal
		self.image = image
		self.greyBlobs_bright = []
		self.greyBlobs_dark = []
		self.coarserLevel = coarserLevel
		self.finerLevel = finerLevel



class Primal_Sketch :
	def __init__(self) :
		self.scale_levels = []

		self.scaleBlobs_bright = []
		self.bifurcation_bright = []

		self.scaleBlobs_dark = []
		self.bifurcation_dark = []



