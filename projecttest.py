from LoadRastRadar import *
import ScaleSpace as ss
import numpy

import pylab

import mpl_toolkits.mplot3d.axes3d as p3

def foo() :

	scales = [0, 1, 2, 4, 6, 8, 16, 32, 64]
	ps = ss.Primal_Sketch()
	Sketch_Radar(ps, "6500KTLX20050514_052255.nc", scales)
	Create_Images(ps, "KTLX")

	# Need to re-assign scales because the array gets emptied by Sketch_Radar
	scales = [0, 1, 2, 4, 6, 8, 16, 32, 64]
	ps = ss.Primal_Sketch()
	Sketch_Radar(ps, "PAR_20090210_202625_rast.nc", scales)
	Create_Images(ps, "MPAR")



def Sketch_Radar(ps, dataFile, scales) :
	radarData = LoadRastRadar(dataFile)
	radarData['vals'] = numpy.nan_to_num(numpy.squeeze(radarData['vals']))

	# Anchoring to zero and type-casting as integer.
	imageData = (radarData['vals'] - radarData['vals'].min()).astype(int)

	ps.CreateSketch(imageData, scales)



def Create_Images(ps, namestem) :
	finalScales = list(ps.scale_levels)
	finalScales.sort()

	if len(finalScales) == 0 :
		print "No Scales!"
		return

	scaleVols = numpy.array([aScaleBlob.volume(finalScales) for aScaleBlob in ps.scaleBlobs_bright])

	# Sort for greatest to least
	sortIndices = numpy.argsort(scaleVols)[::-1]


	#whichBlobs = scaleVols > numpy.mean(scaleVols)
	whichBlobs = scaleVols > 0

	theBlobs = [ps.scaleBlobs_bright[anIndex] for anIndex in sortIndices[whichBlobs]]
	print "Plotting the significant blobs"
	displayRegion = set([])
	figureIndex = 0
	fig = pylab.figure()
	ax = fig.gca()
	xLims = [0, ps.scale_levels[0].image.shape[1]]
	yLims = [0, ps.scale_levels[0].image.shape[0]]

	ax.set_title("%s - Significant Blobs #%d" % (namestem, figureIndex)) 


	for blobIndex in sortIndices :
		if len(displayRegion.intersection(ps.scaleBlobs_bright[blobIndex].approp_support)) > 0 :

			ax.set_xlim(xLims)
			ax.set_ylim(yLims)
			fig.savefig("SignifBlobs%0.3d_%s.png" % (figureIndex, namestem))
			pylab.close()
			figureIndex += 1

			# time for the next figure!
			fig = pylab.figure()
			ax = fig.gca()
			ax.set_title("%s - Significant Blobs #%d" % (namestem, figureIndex))
			displayRegion = set([])

		Plot_ScaleBlob(ax, ps.scaleBlobs_bright[blobIndex])
		displayRegion.update(ps.scaleBlobs_bright[blobIndex].approp_support)

	ax.set_xlim(xLims)
	ax.set_ylim(yLims)
	fig.savefig("SignifBlobs%0.3d_%s.png" % (figureIndex, namestem))
	pylab.close()

	print "Plotting all blobs"
	fig = pylab.figure()
	ax = fig.gca()

	for aBlob in ps.scaleBlobs_bright :
		Plot_ScaleBlob(ax, aBlob)

	ax.set_xlim(xLims)
	ax.set_xlabel("X")
	ax.set_ylim(yLims)
	ax.set_ylabel("Y")
	ax.set_title("%s - All Scale Blobs" % namestem)
	fig.savefig("AllBlobs_%s.png" % namestem)
	pylab.close()

	
	print "Plotting extrema paths"
	fig = Plot_ExtremaPaths(theBlobs, ps.events_bright)
	ax = fig.gca()
	#ax.set_xlim(xLims)
	#ax.set_ylim(yLims)
	ax.set_zlabel("Scale Space")
	ax.set_xlabel("X")
	ax.set_ylabel("Y")
	ax.set_title("%s - Extrema Paths" % namestem)

	fig.savefig("ExtremaPaths_%s.png" % namestem)
	pylab.close()



	



def Plot_ScaleBlob(ax, scaleBlob) :
        x, y, z = Blob2Coords(scaleBlob)
	ax.scatter(x, y, c=z, marker='s', s = 2, edgecolors='none', alpha = 0.05, vmin=0, vmax=64)

def Plot_ExtremaPaths(scaleBlobs, scaleEvents) :
	fig = pylab.figure()
	ax = p3.Axes3D(fig)

	scaleBlobs = set(scaleBlobs)

	for aScaleBlob in scaleBlobs :
		ys, xs = zip(*[aRegion.first_moment() for aRegion in aScaleBlob.support_regions])
		scaleCoords = [aScaleLevel.scaleVal for aScaleLevel in aScaleBlob.scale_levels]

		ax.plot(xs, ys, zs=scaleCoords, color='k')

	# Plot the extrema paths for splits, merges and complex events
	for aScaleEvent in scaleEvents :
		if aScaleEvent.event_type in [ss.Scale_Event.SPLIT,
					      ss.Scale_Event.MERGE,
					      ss.Scale_Event.COMPLEX] :

			# Only plot the merge/split/complex if the scale blobs above and below are already plotted.
			if scaleBlobs.issuperset(aScaleEvent.scaleblobs_above) and scaleBlobs.issuperset(aScaleEvent.scaleblobs_below) :
				# Get all permutations of the paths between the scale blobs above and below the event
				paths = [(scaleBlobAbove.support_regions[-1].first_moment(),
					  scaleBlobBelow.support_regions[0].first_moment()) for scaleBlobAbove in aScaleEvent.scaleblobs_above
										    for scaleBlobBelow in aScaleEvent.scaleblobs_below]
				scaleCoords = [(scaleBlobAbove.scale_levels[-1].scaleVal,
						scaleBlobBelow.scale_levels[0].scaleVal) for scaleBlobAbove in aScaleEvent.scaleblobs_above
										 for scaleBlobBelow in aScaleEvent.scaleblobs_below]

				# Now plot each path
				for aPath, scales in zip(paths, scaleCoords) :
					# The fancy zipping is so that the y coordinates go second while the x coordinates go first
					ax.plot(*zip(*aPath)[::-1], zs=scales, color='r')

	return fig
			

"""
def Make_ScaleSpace_Plot(ps) :
	someColors = ['r', 'g', 'm', 'b', 'c', 'k']
	X, Y = numpy.meshgrid(numpy.arange(1200), numpy.arange(925))


	fig = pylab.figure()
	ax = p3.Axes3D(fig)
	#pylab.hold(True)
	fakeImage = numpy.zeros(X.shape)
	fakeImage[0, 0] = ps.scale_levels[-1].scaleVal

	had_data = ax.has_data()
	
	for aLevel in ps.scale_levels[-2:] :
		print aLevel.scaleVal
		cset = matplotlib.axes.Axes.contourf(ax, X, Y,
						     numpy.ma.masked_array(aLevel.image, aLevel.image < 0),
						     range(aLevel.image.max() + 1))
		levels = cset.levels
		colls = cset.collections
		for z1, linec in zip(levels, colls) :
			art3d.poly_collection_2d_to_3d(linec, aLevel.scaleVal)
			linec.set_sort_zpos(aLevel.scaleVal)
#			linec.set_sort_zpos(z1)

	ax.auto_scale_xyz(X, Y, fakeImage, had_data)

	#ax.set_xlim((X.min(), X.max()))
	#ax.set_ylim((Y.min(), Y.max()))
	#ax.set_zlim3d((0, 10.0))
	pylab.show()
"""

def Blob2Coords(blob) :
        x = []
        y = []
        z = []

	for aRegion, scale_lev in zip(blob.support_regions, blob.scale_levels) :
                ytemp, xtemp = zip(*list(aRegion))
                x += xtemp
                y += ytemp
                z += [scale_lev.scaleVal] * len(xtemp)

        return x, y, z
"""		
def testwork(blob) :
	fig = pylab.figure()
	ax = fig.gca()
	x, y, z = Blob2Coords(blob)

	levelVal = max(z)
	x1, y1 = zip(*[(x[index], y[index]) for index in range(len(z)) if z[index] == levelVal])
	x1 = numpy.asarray(x1)
	y1 = numpy.asarray(y1)

	print x1.shape, y1.shape, levelVal

	X, Y = numpy.meshgrid(range(x1.min(), x1.max() + 2), range(y1.min(), y1.max() + 2))
	C = numpy.empty((X.shape[0] - 1, X.shape[1] - 1))
	C.fill(numpy.nan)
	C[y1 - y1.min(), x1 - x1.min()] = levelVal

	X = numpy.ma.asarray(X)
	Y = numpy.ma.asarray(Y)
	C_masked = numpy.ma.masked_array(C, numpy.isnan(C))

	Ny, Nx = X.shape

	mask = numpy.ma.getmaskarray(C_masked)[0:Ny - 1, 0:Nx - 1]

	ravelmask = (mask==0).ravel()
	X1 = numpy.compress(ravelmask, numpy.ma.filled(X[0:-1,0:-1]).ravel())
	Y1 = numpy.compress(ravelmask, numpy.ma.filled(Y[0:-1,0:-1]).ravel())
	X2 = numpy.compress(ravelmask, numpy.ma.filled(X[1:,0:-1]).ravel())
	Y2 = numpy.compress(ravelmask, numpy.ma.filled(Y[1:,0:-1]).ravel())
	X3 = numpy.compress(ravelmask, numpy.ma.filled(X[1:,1:]).ravel())
	Y3 = numpy.compress(ravelmask, numpy.ma.filled(Y[1:,1:]).ravel())
	X4 = numpy.compress(ravelmask, numpy.ma.filled(X[0:-1,1:]).ravel())
	Y4 = numpy.compress(ravelmask, numpy.ma.filled(Y[0:-1,1:]).ravel())
        npoly = len(X1)

	xy = numpy.concatenate((X1[:,numpy.newaxis], Y1[:,numpy.newaxis],
                             X2[:,numpy.newaxis], Y2[:,numpy.newaxis],
                             X3[:,numpy.newaxis], Y3[:,numpy.newaxis],
                             X4[:,numpy.newaxis], Y4[:,numpy.newaxis],
                             X1[:,numpy.newaxis], Y1[:,numpy.newaxis]),
                             axis=1)

	verts = xy.reshape((npoly, 5, 2))
	C = numpy.compress(ravelmask, numpy.ma.filled(C[0:Ny-1, 0:Nx-1]).ravel())

	collection = mcolls.PolyCollection(verts)
	collection.set_color('r')
	#collection.set_alpha(1.0)
	#collection.set_array(C)
	#collection.set_cmap(None)
	#collection.set_norm(None)
	#collection.set_clim(levelVal - 1.0, levelVal + 1.0)

	#ax.grid(False)
	#corners = (x1.min(), y1.min()), (x1.max(), y1.max())
	ax.set_xlim((x1.min(), x1.max()))
	ax.set_ylim((y1.min(), y1.max()))
	#ax.autoscale_view()
	ax.add_collection(collection)

	pylab.show()
"""
