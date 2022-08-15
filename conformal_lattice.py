import interpGeo.interpolateGeo as iG
import interpGeo.exportPolyline as eP
import interpGeo.nurbs as NURBS
import math as m
import numpy as np
import scipy.interpolate

import re
import sys
import os
import json


"""
MESH CLASS takes in a list of vertices and faces (ALL QUADS OR ALL TRIANGLES) and calculates normals
and face vertices based on the given information

Input: 
verts = numpy array of 3D points
faces = numpy array of 3 or 4 Indicies

"""


class Mesh:

	def __init__(self, verts, faces):

		self.v = verts
		self.f = faces
		self.fv = [] # face vertices
		self.fn = [] # face normals
		self.vn = [] # vertex normals

		self.calculateNormals()

	#### 
	# CACLULATE NORMALS
	# generates the fn, fv and vn values for the mesh
	####

	def calculateNormals(self):

		self.fn = []
		self.fv = []
		self.vn = []

		# 
		# CALCULATE THE FACE NORMALS AND FACE VERTICES
		#
		for i in range(self.f.shape[0]):
			
			face = self.v[self.f[i, :]]
			
			# The face normals are calculated by taking the first two edges and creating a cross product 
			# this is an approximation for quad faces, but accurate for triangular faces
			f_x = np.subtract(face[1], face[0])
			f_y = np.subtract(face[2], face[0])
			n = np.cross(f_x, f_y)
			n = n/np.linalg.norm(n)

			self.fn.append(n)
			self.fv.append(face)

		self.fn = np.array(self.fn)
		self.fv = np.array(self.fv)

		vn = np.zeros(self.v.shape) # initialize vertex normals list

		#
		# CALCULATE THE VERTEX NORMALS
		#
		for i in range(self.v.shape[0]):
			
			count = 0

			locations = np.where(self.f == i)[0] # identify all the faces that this vertex is a part of

			vn[i] = np.sum(self.fn[locations[:]]) # add up all the face normals that vertex is attached to

			self.vn.append(np.multiply(vn[i], 1/len(locations))) # divide by the number of faces to get the average

		self.vn = np.array(self.vn)

	####
	# joinMesh joins another mesh with this mesh
	# Input:
	# mesh = another mesh object
	####
	
	def joinMesh(self, mesh):

		# THIS STRATEGY IGNORES DUPLICATE VERTICES AND DUPLICATE FACES TO JUST APPEND A MESH

		faces = mesh.f[:, :] + len(self.v)

		self.f = np.concatenate((self.f, faces), axis = 0)
		self.v = np.concatenate((self.v, mesh.v), axis = 0)
		self.vn = np.concatenate((self.vn, mesh.vn), axis = 0)

		# THE COMMENTED OUT SECTION IS DESIGNED TO REMOVE DUPLICATE FACES AND VERTICES

		# f = np.zeros(mesh.f.shape)
		# f.fill(-1)
		# fv = np.zeros(mesh.fv.shape)
		# fv.fill(-1)
		# v = []
		# vn = []
		# f_join = []
		# fv_join = []

		# for i in range(self.v.shape[0]):

		# 	distances = np.linalg.norm(mesh.v[:] - self.v[i], axis = 1)
		# 	duplicates = np.where(distances < .01)[0]

		# 	for j in range(duplicates.shape[0]):
		# 		locations = np.where(f == duplicates[j])[0]
		# 		f[locations[:]] = i
		# 		fv[locations[:]] = self.v[i]

		# for i in range(f.shape[0]):
		# 	for j in range(f.shape[1]):
		# 		if f[i, j] == -1:
		# 			v.append(mesh.v[mesh.f[i, j]])
		# 			vn.append(mesh.vn[mesh.f[i, j]])
		# 			f[i, j] = self.v.shape[0] + len(v) - 1
		# 			fv[i, j] = mesh.fv[i, j]

		# for i in range(f.shape[0]):
		# 	duplicates = np.where(self.f == f[i])[0]
		# 	if duplicates.shape[0] == 0:
		# 		f_join.append(f[i])
		# 		fv_join.append(fv[i])

		# v = np.array(v)
		# vn = np.array(vn)
		# f_join = np.array(f_join)
		# fv_join = np.array(fv_join)

		# self.f = np.concatenate((self.f, f_join), axis = 0)
		# self.fv = np.concatenate((self.fv, fv_join), axis = 0)
		# self.v = np.concatenate((self.v, v), axis = 0)
		# self.vn = np.concatenate((self.vn, vn), axis = 0)

		return [self.v, self.f]


		
"""

UNIT CLASS is very similar to the mesh class but it is used for curves
Input
verts = numpy array of 3D points
members = list of 2D indicies connecting 2 points

"""


class unit:

	def __init__(self, verts = np.array([]), members = []):

		self.v = verts
		self.members = members
		self.els = []

	def genMembers(self):

		self.els = []

		for i in range(len(self.members)):

			member = []

			for j in range(len(self.members[i])):

				member.append(self.v[self.members[i][j]])

			self.els.append(member)

		return self.els


class Attractor:

	def __init__(self, lattice, cnt, threshold, unit):

		self.lattice = lattice
		self.blend = unit
		self.cnt = cnt
		self.thres = threshold



"""
CONFORMAL LATTICE CLASS uses up to 2-6 surfaces to form a lattice based on user defined boundary surfaces
Input:
guides (list[str]) = file locations of two guiding surfaces
formV (list[str]) = file locations of up to two surfaces that warp the V direction to meet the V boundary
formW (list[str]) = file locations of up to two surfaces that warp the W direction to meet the W boundary
divU (int) = the number of units in the U direction
divV (int) = the number of units in the V direction
divW (int) = the number of units in the W direction

"""


class conformalLattice:


	def __init__(self, guides, formV, formW, divU, divV, divW):

		self.dU = divU
		self.dV = divV
		self.dW = divW

		self.grid = []
		self.cells = []
		self.v_pts = []
		self.w_pts = []

		self.formV = []
		self.formW = []

		self.attractors = []

		self.rails = [NURBS.surface(guides[0]), NURBS.surface(guides[1])]  # the surfaces that generate the base grid

		for i in range(len(formV)):
			self.formV.append(NURBS.surface(formV[i]))  # the surfaces that warp the base grid in the V direction
		
		for i in range(len(formW)):
			self.formW.append(NURBS.surface(formW[i]))  # the surfaces that warp the base grid in the W direction


	####
	# createGrid(self)
	# takes the two rail surfaces and "lofts" between their UV coordinates to form the W axis and the base grid
	####

	def createGrid(self):

		self.rails[0].extractGrid(self.dU, self.dV)
		self.rails[1].extractGrid(self.dU, self.dV)

		w_vecs = np.subtract(self.rails[1].grid[:], self.rails[0].grid[:]) # create vectors between UV rail 1 and rail 2
		w_vecs = w_vecs[:,:]/self.dW  # divide the connecting vectors by the W dimension

		self.grid = []

		for i in range(self.dW+1):

			self.grid.append(np.add(self.rails[0].grid[:,:], w_vecs[:,:]*i))

		self.grid = np.array(self.grid)
		self.grid = np.swapaxes(self.grid, 0, 2)  # swap axes to make sure that the shape remains U, V, W


	####
	# warpGrid(self, factor)
	# warps the grid based on the give form V surfaces and the form W surfaces
	# INPUT:
	# factor (float) = optional float to determine how much the interior grid warps relative to the boundary surface
	####


	def warpGrid(self, factor=1):

		# USES THE V SURFACES TO WARP THE GRID

		for i in range(len(self.formV)):

			self.v_pts.append(self.formV[i].extractGrid(self.dW, self.dV))

			if i==0:
				self.grid[0, :, :] = self.v_pts[i] # moves the first U layer to the first boundary surface
				print('FIRST U LAYER WARPED')
			if i==1:
				self.grid[-1, :, :] = self.v_pts[i] # moves the last U layer to the second boundary surface
				print('LAST U LAYER WARPED')

			# CREATES POS AND NEG VECTORS FOR EACH LAYER OF THE GRID IN THE U DIRECTION (ONE LAYER TO THE NEXT)

			for j in range(self.grid.shape[0]-2):

				pos_vec = (self.grid[j+2, :, :] - self.grid[j+1, :, :])*abs(j/self.grid.shape[0])*factor
				neg_vec = (self.grid[j, :, :] - self.grid[j+1, :, :])*abs(1 - j/self.grid.shape[0])*factor

				# as the layers move from the first to the last the pos vector gets weaker and the neg vector gets stronger

				self.grid[j+1, :, :] = self.grid[j+1, :, :] + pos_vec + neg_vec


		# USES THE W SURFACES TO WARP THE GRID

		for i in range(len(self.formW)):

			self.w_pts.append(self.formW[i].extractGrid(self.dW, self.dU))

			if i==0:
				self.grid[:, 0, :] = self.w_pts[0] # moves the first V layer to the first boundary surface
				print('FIRST V LAYER WARPED')
			if i==1:
				self.grid[:, -1, :] = self.w_pts[1] # moves the first V layer to the second boundary surface
				print('LAST V LAYER WARPED')

			# CREATES POS AND NEG VECTORS FOR EACH LAYER OF THE GRID IN THE V DIRECTION (ONE LAYER TO THE NEXT)

			for j in range(self.grid.shape[1]-2):

				pos_vec = (self.grid[:, j+2, :] - self.grid[:, j+1, :])*abs(j/self.grid.shape[1])*factor
				neg_vec = (self.grid[:, j, :] - self.grid[:, j+1, :])*abs(1 -j/self.grid.shape[1])*factor

				# as the layers move from the first to the last the pos vector gets weaker and the neg vector gets stronger

				self.grid[:, j+1, :] = self.grid[:, j+1, :] + pos_vec + neg_vec


	####
	# addMerges(self, factor)
	# merges a specified number of rows and columns in 2D along an axis for a specified range of units along that axis 
	# INPUT:
	# travel_axi (int) = a number indicating the direction of travel (U = 0, V = 1, W = 2)
	# travel_center (int) = the location of the row and column where surrounding rows and columns will be merged towards
	# travel_range ([float, float]) = the start and end parameters along the travel axis
	# merge_range ([int, int]) = the number of rows-columns that will merge together towards the travel_center
	####

	def addMerges(self, travel_axi, travel_center, travel_range, merge_range = [2,2]):

		travel_st = int(travel_range[0]*self.grid.shape[travel_axi])  # calculate start index along travel
		travel_en = int(travel_range[1]*self.grid.shape[travel_axi])  # calculate end index along travel

		cnt = travel_center

		# SCROLL THROUGH EVERY GRID "ROW" IN THE TRAVEL-AXIS

		for i in range(self.grid.shape[travel_axi]):

			if i>=travel_st and i<=travel_en:

				if travel_axi == 0:  # if travel axis is in the U direction

					# SCROLL THROUGH MERGE RANGE IN ROWS AND COLUMNS RIGHT AND LEFT TO MERGE

					for j in range(merge_range[0]):

						merge_center = cnt[0] + j

						if merge_center < self.grid.shape[1]:

							for k in range(merge_range[1]):

								if cnt[1] + k < self.grid.shape[2]:

									self.grid[i, merge_center, cnt[1] + k] = self.grid[i, cnt[0], cnt[1] + k]

								if cnt[1] - k > 0:

									self.grid[i, merge_center, cnt[1] - k] = self.grid[i, cnt[0], cnt[1] - k]

						
						merge_center = cnt[0] - j

						if merge_center > 0:

							for k in range(merge_range[1]):

								if cnt[1] + k < self.grid.shape[2]:

									self.grid[i, merge_center, cnt[1] + k] = self.grid[i, cnt[0], cnt[1] + k]

								if cnt[1] - k > 0:

									self.grid[i, merge_center, cnt[1] - k] = self.grid[i, cnt[0], cnt[1] - k]




				if travel_axi == 1:

					for j in range(merge_range[0]):

						merge_center = cnt[0] + j

						if merge_center < self.grid.shape[0]:

							for k in range(merge_range[1]):

								if cnt[1] + k < self.grid.shape[2]:

									self.grid[merge_center, i, cnt[1] + k] = self.grid[cnt[0], i, cnt[1] + k]

								if cnt[1] - k > 0:

									self.grid[merge_center, i, cnt[1] - k] = self.grid[cnt[0], i, cnt[1] - k]

						
						merge_center = cnt[0] - j

						if merge_center > 0:

							for k in range(merge_range[1]):

								if cnt[1] + k < self.grid.shape[2]:

									self.grid[merge_center, i, cnt[1] + k] = self.grid[cnt[0], i, cnt[1] + k]

								if cnt[1] - k > 0:

									self.grid[merge_center, i, cnt[1] - k] = self.grid[cnt[0], i, cnt[1] - k]



				if travel_axi == 2:

					for j in range(merge_range[0]):

						merge_center = cnt[0] + j

						if merge_center < self.grid.shape[0]:

							for k in range(merge_range[1]):

								if cnt[1] + k < self.grid.shape[1]:

									self.grid[merge_center, cnt[1] + k, i] = self.grid[cnt[0], cnt[1] + k, i]

								if cnt[1] - k > 0:

									self.grid[merge_center, cnt[1] - k, i] = self.grid[cnt[0], cnt[1] - k, i]

						
						merge_center = cnt[0] - j

						if merge_center > 0:

							for k in range(merge_range[1]):

								if cnt[1] + k < self.grid.shape[1]:

									self.grid[merge_center, cnt[1] + k, i] = self.grid[cnt[0], cnt[1] + k, i]

								if cnt[1] - k > 0:

									self.grid[merge_center, cnt[1] - k, i] = self.grid[cnt[0], cnt[1] - k, i]


	def addAttractor(self, attPt, thres, unit):

		self.attractors.append(Attractor(self, attPt, thres, unit))

	def applyAttractors(self):

		attCenters = []

		for i in range(len(self.attractors)):
			attCenters.append(self.attractors[i].cnt)

		attCenters = np.array(attCenters)
		
		limiter = 3*self.attractors[0].thres
		if (limiter < np.max(self.cell_radii)): limiter = np.max(self.cell_radii)

		affected = 0
		total = 0

		for i in range(self.cells.shape[0]):
			for j in range(self.cells.shape[1]):
				for k in range(self.cells.shape[2]):
					if np.min(np.linalg.norm(attCenters[:] - self.units[i, j, k].v[0], axis = 1)) < limiter:
						self.units[i, j, k].v = self.cells[i, j, k].blendSrc(self.attractors[0].blend, attCenters, self.attractors[0].thres)
						affected = affected + 1
					total = total + 1

	####
	# createCells(self, unit, genBoxes)
	# creates basic boxes and their corresponding units
	# unit (mesh or unit object) = a mesh/unit object that contains vertices and faces/members
	# genBoxes (boolean) = Generate prievew boxes to preview the cells
	####

	def createCells(self, unit, genBoxes = True):

		self.unit = unit
		self.units = []

		self.cell_radii = []

		points = unit.v
		if len(unit.f) != 0:
			faces = self.unit.f
		else:
			faces = []
		boxes = []

		# SCROLLS THROUGH EVERY PART OF THE GRID AND GENERATES BOXES FROM 8-PTS

		for i in range(self.grid.shape[0] - 1):

			cell_layers = []
			unit_layers = []

			for j in range(self.grid.shape[1] - 1):

				cell_cols = []
				cell_units = []

				for k in range(self.grid.shape[2] - 1):

					pt01 = self.grid[i, j, k]
					pt02 = self.grid[i+1, j, k]
					pt03 = self.grid[i+1, j+1, k]
					pt04 = self.grid[i, j+1, k]

					pt05 = self.grid[i, j, k+1]
					pt06 = self.grid[i+1, j, k+1]
					pt07 = self.grid[i+1, j+1, k+1]
					pt08 = self.grid[i, j+1, k+1]

					crns = np.array([pt01, pt02, pt03, pt04, pt05, pt06, pt07, pt08])

					cell_radius = np.max(np.linalg.norm(crns[:] - np.mean(crns, 0), axis = 1))*2
					self.cell_radii.append(cell_radius)

					# create a box morph object based on the corners and points from the vertices in the base unit
					cell_cols.append(iG.boxMorph(crns, points, faces))
					
					# if genBoxes create preview object cells
					if genBoxes:
						boxes.extend(cell_cols[-1].genDefault())
					
					if len(faces) > 0:
						cell_units.append(Mesh(cell_cols[-1].v, np.array(cell_cols[-1].f)))

				cell_layers.append(cell_cols)
				
				if len(faces)>0:
					unit_layers.append(cell_units)

			self.cells.append(cell_layers)
			self.units.append(unit_layers)

		self.cells = np.array(self.cells)
		self.units = np.array(self.units)

		print('all units read')

		if genBoxes:
			eP.exportPolyline(boxes, 'test_boxes.obj')


	####
	# genMeshUnits(self, loc)
	# joins all the individual meshes created for each box into a single mesh 
	# loc (str) = path location of mesh file to be written (MUST BE .OBJ)
	####

	def genMeshUnits(self, loc = 'test_lattice.obj'):

		self.mesh = Mesh(self.units[0, 0, 0].v, self.units[0, 0, 0].f)

		if self.units.shape[0] != 0:

			for i in range(self.units.shape[0]):
				for j in range(self.units.shape[1]):
					for k in range(self.units.shape[2]):

						if i>0 or j>0 or k>0:
							if self.cells[i, j, k].W > 0 and self.cells[i, j, k].L > 0:
								self.mesh.joinMesh(self.units[i, j, k])

			print(len(self.mesh.f))
			print(self.mesh.v.shape)
			print(self.mesh.vn.shape)
			iG.exportMeshOBJ(self.mesh.f, self.mesh.v, self.mesh.vn, loc)




def Main(base_unit, guides, formV = [], formW = [], blend_unit = ""):

	# guide = ['NURBS/srf_01.stp', 'NURBS/srf_02.stp']
	# formV = ['NURBS/srf_03.stp', 'NURBS/srf_04.stp']
	# formW = ['NURBS/srf_05.stp', 'NURBS/srf_06.stp']

	lattice = conformalLattice(guides, formV, formW, 56, 2, 1)
	lattice.createGrid()
	lattice.warpGrid(.25)

	imported = iG.importMeshOBJ(base_unit)
	baseUnit = Mesh(imported[0], np.array(imported[1]))

	lattice.createCells(baseUnit, False)

	if (blend_unit != ""):

		imported_blend = iG.importMeshOBJ(blend_unit)
		blendUnit = Mesh(imported_blend[0], np.array(imported_blend[1]))

		lattice.addAttractor([25, 0, 1], 25, blendUnit.v)
		lattice.addAttractor([-25, 0, 4], 25, blendUnit.v)
		lattice.applyAttractors()

	# addMerges(self, travel_axi, travel_center, travel_range, merge_axi, merge_range):

	# axis moved across U, starts at 4V,3W location, moves from 4/14 to 10/14 in U axis, merging 4V (rows) and 3W (columns)

	lattice.addMerges(0, [4, 1], [4/14, 10/14], [1,0])

	lattice.genMeshUnits("woven_blend_ring.obj")
	

Main("unit_tri.obj", ["inner_scaffold.igs", "outer_scaffold.igs"], blend_unit = "unit_blend.obj")