#!/usr/bin/env python

from LoadRastRadar import *
import ScaleSpace as ss
import numpy

import pylab
import matplotlib.collections as mcolls
import matplotlib.axes
import mpl_toolkits.mplot3d.axes3d as p3
import mpl_toolkits.mplot3d.art3d as art3d

scales = [41, 21, 11, 5, 1, 0]

def foo(ps) :
	radarData = LoadRastRadar("6500KTLX20050514_052255.nc")
	radarData['vals'] = numpy.nan_to_num(numpy.squeeze(radarData['vals']))

	print "Data min:", radarData['vals'].min(), "   Data max:", radarData['vals'].max()

	imageData = (radarData['vals'] - radarData['vals'].min()).astype(int)

	print "Image max: ", imageData.max()


	ps.CreateSketch(imageData, scales)

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
