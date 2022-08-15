import math as m
import numpy as np
import scipy.interpolate
import open3d as o3d
import copy

import re
import sys
import os
import json


class Warp:

	def __init__(self, mesh):

		self.mesh = mesh
		self.v = np.asarray(mesh.vertices)
		self.f = np.asarray(mesh.triangles)
		self.mesh.compute_vertex_normals()
		self.vn = np.asarray(mesh.vertex_normals)
		self.attractors = []

	def define_valid_vertices(self, zRange):

		maxY = np.max(self.v[:, 1])

		flattened = np.zeros(self.v.shape)
		flattened[:] = self.v[:]
		flattened[:, 2] = 0

		self.valid_indexes = np.where(np.linalg.norm(flattened, axis = 1) > maxY - .5)[0]
		self.v_valid = self.v[self.valid_indexes]
		self.valid_indexes = self.valid_indexes[np.where(self.v_valid[:, 2] < zRange[1])[0]]
		self.v_valid = self.v[self.valid_indexes]
		self.valid_indexes = self.valid_indexes[np.where(self.v_valid[:, 2] > zRange[0])[0]]
		self.v_valid = self.v[self.valid_indexes]

		print(self.valid_indexes.shape)
		print(self.v_valid.shape)

	def addAttractor(self, curve):

		self.attractors.append(curve)

	def imprint(self, thres = .5, depth = .5):

		attPts = []

		for i in range(len(self.attractors)):
			attPts.extend(self.attractors[i])

		attPts = np.array(attPts)

		print(attPts.shape)

		for i in range(self.valid_indexes.shape[0]):
			index = self.valid_indexes[i]
			distances = np.linalg.norm(attPts[:] - self.v[index], axis = 1)
			if np.min(distances) < thres:
				self.v[index] = self.v[index] - self.vn[index] * (1 - np.min(distances)/thres) * depth

	def write_mesh(self, output):

		self.output = copy.deepcopy(self.mesh)
		self.output.vertices = o3d.utility.Vector3dVector(self.v)
		self.output.compute_vertex_normals()

		o3d.io.write_triangle_mesh(output, self.output)



def main(path, radius = 12):

	mesh = o3d.io.read_triangle_mesh(path)
	myWarp = Warp(mesh)
	myWarp.define_valid_vertices([1, 3.75])

	spacing = 2
	resoU = int(2 * m.pi * radius/spacing)

	for i in range(resoU):
		pts = []
		pts02 = []
		for j in range(10):
			pts.append([radius * m.sin(i/resoU * 2 * m.pi), radius * m.cos(i/resoU * 2 * m.pi), 1 + j/9 * 2.75 * m.sin(i/resoU * 2 * m.pi * 3)])
			# pts02.append([radius * m.sin(i/resoU * 2 * m.pi), radius * m.cos(i/resoU * 2 * m.pi), 3.75 - j/9 * 2 * m.cos(i/resoU * 2 * m.pi * 3)])
		# curve02 = iG.interpCrv(np.array(pts02))
		myWarp.addAttractor(pts)
		# myWarp.addAttractor(curve02)
	
	myWarp.imprint()
	myWarp.write_mesh("test_imprint.obj")


main("ring_300x100.obj")