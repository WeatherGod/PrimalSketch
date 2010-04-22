#!/usr/bin/env python

from LoadRastRadar import *
import ScaleSpace as ss
import numpy

import pylab
import mpl_toolkits.mplot3d.axes3d as p3

scales = [1, 5, 11, 21, 41]

def foo(ps) :
	radarData = LoadRastRadar("6500KTLX20050514_052255.nc")
	radarData['vals'] = numpy.nan_to_num(numpy.squeeze(radarData['vals']))

	print "Data min:", radarData['vals'].min(), "   Data max:", radarData['vals'].max()

	imageData = (radarData['vals'] - radarData['vals'].min()).astype(int)

	print "Image max: ", imageData.max()


	ps.CreateSketch(imageData, scales)

def Make_ScaleSpace_Plot(ps) :
	someColors = ['r', 'g', 'm', 'b', 'c', 'k']
	#X, Y = numpy.meshgrid(numpy.arange(1200), numpy.arange(925))


	fig = pylab.figure()
	ax = p3.Axes3D(fig)
	pylab.hold(True)

	for index, aBlob in enumerate(ps.scaleBlobs_bright) :
		print index
		#Z = numpy.zeros(Y.shape)
		x, y, z = Blob2Coords(aBlob)
		#Z[y, x] = z
		if len(x) > 300 :
			ax.scatter3D(x, y, z, color = someColors[index % 6])


	pylab.show()


def Blob2Coords(blob) :
        x = []
        y = []
        z = []

	for aRegion, scale_lev in zip(blob.support_regions, blob.scale_levels) :
                ytemp, xtemp = zip(*list(aRegion.pixels))
                x += xtemp
                y += ytemp
                z += [scale_lev.scaleVal] * len(xtemp)

        return x, y, z
		
def testwork(blob) :
	pylab.figure()
	x, y, z = Blob2Coords(blob)

	levelVal = max(z)
	x1, y1 = zip(*[(x[index], y[index]) for index in range(len(z)) if z[index] == levelVal])
	x1 = numpy.asarray(x1)
	y1 = numpy.asarray(y1)

	print x1.shape, y1.shape, levelVal

	X, Y = numpy.meshgrid(range(min(x1), max(x1) + 2), range(min(y1), max(y1) + 2))
	C = numpy.empty((X.shape[0] - 1, X.shape[1] - 1)).fill(numpy.nan)
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

	xy = np.concatenate((X1[:,numpy.newaxis], Y1[:,numpy.newaxis],
                             X2[:,numpy.newaxis], Y2[:,numpy.newaxis],
                             X3[:,numpy.newaxis], Y3[:,numpy.newaxis],
                             X4[:,numpy.newaxis], Y4[:,numpy.newaxis],
                             X1[:,numpy.newaxis], Y1[:,numpy.newaxis]),
                             axis=1)

	verts = xy.reshape((npoly, 5, 2))
	C = numpy.compress(ravelmask, numpy.ma.filled(C[0:Ny-1, 0:Nx-1]).ravel())

	collection = pylab.collections.PolyCollection(verts, {'edgecolors': (0, 0, 0, 1), 'antialiaseds': (0,), 'linewidths': (0.25,)})
	collection.set_alpha(1.0)
	collection.set_array(C)
	collection.set_cmap(None)
	collection.set_norm(None)
	collection.set_clim(levelVal - 1, levelVal + 1)

	pylab.grid(False)
	corners = (x1.min(), y1.min()), (x1.max(), y1.max())
	pylab.update_datalim(corners)
	pylab.autoscale_view()
	pylab.add_collection(collection)

	pylab.show()
