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

	def cymatic_imprint(self, attPts, frequency, amp =.5, thres_angle = 90):

		if(thres_angle > 2*m.pi) : thres_angle = thres_angle * m.pi/180

		bottom = .5 #np.min(self.v[self.valid_indexes[:], 2])
		top = 4.5 #np.max(self.v[self.valid_indexes[:], 2])
		thres_rim = 2

		print(bottom)
		print(top)

		for i in range(self.valid_indexes.shape[0]):
			index = self.valid_indexes[i]
			distances = np.linalg.norm(attPts[:] - self.v[index], axis = 1)
			attPt = attPts[np.argmin(distances)]
			attVector = attPt/np.linalg.norm(attPt)
			currVector = self.v[index]/np.linalg.norm(self.v[index])
			ang = m.acos(np.dot(attVector, currVector))
			if ang < thres_angle:
				
				wave = (m.cos(2 * m.pi/frequency * ang) + 1)/2
				factor = (1 - ang/thres_angle) * amp * wave

				if self.v[index, 2] - bottom < thres_rim: 
					amplitude = (self.v[index, 2] - bottom)/thres_rim * factor
				elif top - self.v[index, 2] < thres_rim:
					amplitude = (top - self.v[index, 2])/thres_rim * factor
				else:
					amplitude = factor
				
				self.v[index] = self.v[index] - self.vn[index] * amplitude

	def write_mesh(self, output):

		self.output = copy.deepcopy(self.mesh)
		self.output.vertices = o3d.utility.Vector3dVector(self.v)
		self.output.compute_vertex_normals()

		o3d.io.write_triangle_mesh(output, self.output)



def main(path, freq, thres_angle, radius = 12):

	mesh = o3d.io.read_triangle_mesh(path)
	myWarp = Warp(mesh)
	myWarp.define_valid_vertices([1, 3.75])

	# for i in range(resoU):
	# 	pts = []
	# 	pts02 = []
	# 	for j in range(10):
	# 		pts.append([radius * m.sin(i/resoU * 2 * m.pi), radius * m.cos(i/resoU * 2 * m.pi), 1 + j/9 * 2.75 * m.sin(i/resoU * 2 * m.pi * 3)])
	# 	myWarp.addAttractor(pts)

	attractors = []
	resoU = 6
	for i in range(resoU): 
		z = 0
		if i%2 == 0: z = 5
		attractors.append([radius * m.sin(i/(resoU-1) * 2 * m.pi), radius * m.cos(i/(resoU-1) * 2 * m.pi), z])
	attractors = np.array(attractors)
	
	myWarp.cymatic_imprint(attractors, freq, thres_angle = thres_angle)
	myWarp.write_mesh("test_cymatic.obj")


main("ring_300x100.obj", freq = .1, thres_angle = 180)