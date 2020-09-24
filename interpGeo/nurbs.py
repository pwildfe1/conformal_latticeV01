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


class surface:

	def __init__(self, path):

		self.srf = readSurface(path)


	def evalSurface(self, u, v):

		minU, maxU, minV, maxV = self.srf.ParameterRange

		U = u*(maxU-minU) + minU
		V = v*(maxV-minV) + minV

		point = self.srf.valueAt(U, V)

		return point

	def extractGrid(self, dimU, dimV):

		U = []
		V = []
		normals = []

		for i in range(dimU+1):

			col = []
			n_col = []

			for j in range(dimV+1):

				u, v = i/dimU, j/dimV

				col.append(self.evalSurface(u, v))
				n_col.append(self.normalAt(u,v))

			U.append(col)
			normals.append(n_col)

		self.grid = np.array(U)
		self.norm = np.array(normals)

		self.grid = np.swapaxes(self.grid, 0, 1)
		self.norm = np.swapaxes(self.norm, 0, 1)

		return self.grid