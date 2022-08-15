import math as m
import numpy as np
import scipy.interpolate
import open3d as o3d
import copy

import re
import sys
import os
import json


class Attractor:

	def __init__(self, source, threshold):
		
		self.source = source
		self.thres = threshold

	def calculate_factors(self, locations):

		factors = np.zeros(locations.shape[0])
		valid_indexes = []

		for i in range(self.source.shape[0]):
			distances = np.linalg.norm(locations[:] - self.source[i], axis = 1)
			if np.min(distances) < self.thres:
				factors[np.where(distances < self.thres)[0]] = 1

		valid_indexes = np.where(factors > 0)[0]

		for i in range(valid_indexes.shape[0]):
			loc = locations[valid_indexes[i]]
			distances = np.linalg.norm(self.source[:] - loc, axis = 1)
			factors[valid_indexes[i]] = np.min(distances)/self.thres

		return factors




class Blend:

	def __init__(self, mesh, blend, warped = True):

		self.mesh = mesh
		self.blend = blend
		self.v = np.asarray(self.mesh.vertices)
		self.blend_v = np.asarray(self.blend.vertices)
		self.output_v = np.asarray(self.mesh.vertices)
		self.v_relative = np.zeros(self.v.shape)
		self.blend_relative = np.zeros(self.blend_v.shape)
		self.attractors = []
		self.warped = warped

	def calculate_blend_vectors(self, blend_box = np.array([])):

		if blend_box.shape[0] > 0:
			self.calculate_box(blend_box)
		else:
			self.blend_vectors = self.blend_v[:] - self.v[:]

	def addAttractor(self, pts, threshold):

		self.attractors.append(Attractor(pts, threshold))


	def blend_mesh(self, thres):

		attPts = []

		for i in range(len(self.attractors)):
			attPts.extend(list(self.attractors[i].source))

		attPts = np.array(attPts)

		for i in range(self.blend_vectors.shape[0]):
			distances = np.linalg.norm(attPts[:] - self.v[i], axis = 1)
			factor = np.min(distances)/thres
			if factor < 1:
				self.output_v[i] = self.v[i] + self.blend_vectors[i]*(1 - factor)


	def write_mesh(self, output):

		self.output = copy.deepcopy(self.mesh)
		self.output.vertices = o3d.utility.Vector3dVector(self.output_v)
		self.output.compute_vertex_normals()

		o3d.io.write_triangle_mesh(output, self.output)

def main(path, blend_path, thres = 10, radius = 12):

	mesh = o3d.io.read_triangle_mesh(path)
	blend = o3d.io.read_triangle_mesh(blend_path)
	myBlend = Blend(mesh, blend)
	myBlend.calculate_blend_vectors()

	resoU = 4

	for i in range(resoU):
		pts = []
		pts.append([radius * m.sin(i/resoU * 2 * m.pi), radius * m.cos(i/resoU * 2 * m.pi), 3])
		myBlend.addAttractor(np.array(pts), thres)
	
	myBlend.blend_mesh(thres)
	myBlend.write_mesh("test_blend.obj")


main("base_woven_ring.obj", "blend_woven_ring.obj")