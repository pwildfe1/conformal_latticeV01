FREECADPATH = 'C:/Users/wildf/Documents/Software/FreeCAD_0.19.22474-Win-Conda_vc14.x-x86_64/bin' # path to your FreeCAD.so or FreeCAD.dll file
import sys
sys.path.append(FREECADPATH)
import FreeCAD
import Part
import numpy as np
import scipy


def readSurface(path):

	obj = Part.Shape()
	obj.read(path)

	srf = obj.Faces[0]

	return srf


"""
SURFACE reads a surface .stp or .iges file for reading and grid extraction
Input:
path (str) = file location of the .iges or .stp file
"""

class surface:

	def __init__(self, path):

		self.srf = readSurface(path)

	####
	# evalSurface (u, v)
	# evaluates the surface at a normalized parameter to extract the point location and normal
	# Input:
	# u (float) = normalized value between 0 - 1 indicating the relative location of evaluation in the u value
	# v (float) = normalized value between 0 - 1 indicating the relative location of evaluation in the v value
	####

	def evalSurface(self, u, v):

		minU, maxU, minV, maxV = self.srf.ParameterRange

		U = u*(maxU-minU) + minU
		V = v*(maxV-minV) + minV

		point = self.srf.valueAt(U, V)
		normal = self.srf.normalAt(U, V)

		return [point, normal]


	####
	# extractGrid (dimU, dimV)
	# extracts a grid of points and normals that is dimU by dimV
	# Input:
	# dimU (int) = number of points in the U direction
	# dimV (int) = number of points in the V direction
	####

	def extractGrid(self, dimU, dimV):

		U = []
		V = []
		normals = []

		for i in range(dimU+1):

			col = []
			n_col = []

			for j in range(dimV+1):

				u, v = i/dimU, j/dimV

				info = self.evalSurface(u, v)

				col.append(info[0])
				n_col.append(info[1])

			U.append(col)
			normals.append(n_col)

		self.grid = np.array(U)
		self.norm = np.array(normals)

		self.grid = np.swapaxes(self.grid, 0, 1)
		self.norm = np.swapaxes(self.norm, 0, 1)

		return self.grid