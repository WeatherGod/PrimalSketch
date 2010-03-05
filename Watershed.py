


class GreyLevel_Blob :
	def __init__(self) :
		self.polarity = True
		self.scale_level = None
		self.extremum = None
		self.saddle = None
		self.support_region = None



class Extremum_Region :
	def __init__(self, pos = [], val = None, greyblob = None) :
		self.positions = pos
		self.grey_val = val
		self.grey_blob = greyblob

	def isNeighbor(pos) :
		for aPoint in self.positions :
			if ((aPoint[0] - pos[0])**2 + (aPoint[1] - pos[1])**2 <= 1.0) :
				return True
		return False

	def addPoint(pos) :
		self.positions.append(pos)

class Saddle_Region :
	def __init__(self, pos = [], val = None, greyblobs = []) :
		self.position = pos
		self.grey_val = val
		self.grey_blobs = greyblobs


class Support_Region :
	def __init__(self) :
		self.blob_area = None
		self.first_moment = None
		self.second_moment = None
		self.pixels = None
		self.atBoundary = None



def Watershed_Transform(image) :
	# Assume that image is a normalized, quantized image...
	pixels = [[]] * (image.max() + 1)

	# loading the 'pixels' array with coordinates of pixels with associated values
	for x in range(image.shape[0]) :
		for y in range(image.shape[1]) :
			pixels[image[x, y]].append((x, y))

	centers = [[]] * (image.max() + 1)
	for q in range(image.max(), -1, -1) :
		for p in pixels[q] :
			isCenter = True
			for aCenter in centers[q] :
				if aCenter.isNeighbor(p) :
					aCenter.addPoint(p)
					isCenter = False
					break
			
			if isCenter : centers[q].append(Extremum_Region([p], q))
			
	# Zero for background, -1 for unchecked, positive values for blob number
	basins = -1 * numpy.ones_like(image)

	for depth in range(



"""
int starting_delta = (SMALLEST_VALID)? 0 : alg->myDelta;
  for (int delta = starting_delta; delta <= alg->myDelta; ++delta){
     std::vector<Center> deferredFromLast, deferredToNext;
     for (int bin=maxbin; bin >= 0; --bin){
        int bin_lower = bin - delta;
        deferredFromLast = deferredToNext; deferredToNext.clear();
        std::vector<Glob> foothills;
        size_t n_centers = centers[bin].size();
        size_t tot_centers = n_centers + deferredFromLast.size();
        for (size_t i = 0; i < tot_centers; ++i){
           // done this way to minimize memory overhead of maintaining two lists
           const Center& center = (i < n_centers)? centers[bin][i] : deferredFromLast[i-n_centers];
           if ( bin_lower < 0 ) bin_lower = 0;
           if ( marked[center.x][center.y] == UNMARKED ){
             bool captured=setMaximum( data, marked, center, bin_lower , alg->myMinSize, centers, foothills );
             if (!captured){
               // decrement to lower value to see if it'll get big enough
               Center defer( center.x, center.y, center.bin-1 );
               deferredToNext.push_back(defer);
             }
           }
        }// all centers
        ErrorLogInfo("Finished processing " << tot_centers << " potential maxima at bin=" << bin << " and delta=" << delta << "\n");
        // this is the last one for this bin
        removeFoothills( data, marked, bin, bin_lower, centers, foothills );
     } // all bins
  } // all deltas
"""
