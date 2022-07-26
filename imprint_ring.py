import interpGeo.interpolateGeo as iG
import interpGeo.exportPolyline as eP
import interpGeo.nurbs as NURBS
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

	def imprint(self, curve, thres = 2, depth = .5):

		pts = np.array(curve.pts)
		crv_cnt = np.mean(pts, 0)
		z_size = np.max(pts[:, 2]) - np.min(pts[:, 2])

		if z_size == 0:
			z_size = 1000

		for i in range(self.valid_indexes.shape[0]):
			index = self.valid_indexes[i]
			distances = np.linalg.norm(pts[:] - self.v[index], axis = 1)
			if np.min(distances) < thres:
				self.v[index] = self.v[index] - self.vn[index] * (1 - np.min(distances)/thres) * depth

	def write_mesh(self, output):

		self.output = copy.deepcopy(self.mesh)
		self.output.vertices = o3d.utility.Vector3dVector(self.v)
		self.output.compute_vertex_normals()

		o3d.io.write_triangle_mesh(output, self.output)



def main(path):

	mesh = o3d.io.read_triangle_mesh(path)
	myWarp = Warp(mesh)
	myWarp.define_valid_vertices([.75, 4.25])

	pts = []

	for i in range(10):
		pts.append([12, 0, .75 + i/9*3])

	curve = iG.interpCrv(np.array(pts))
	myWarp.imprint(curve)
	myWarp.write_mesh("test_imprint.obj")


main("ring_construction.obj")