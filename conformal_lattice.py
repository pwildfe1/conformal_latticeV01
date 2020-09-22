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

		self.rails = [NURBS.surface(guides[0]), NURBS.surface(guides[1])]

		if len(formV) > 0:

			if len(formV) != 2:
				print("YOU NEED EXACTLY TWO V FORMING GUIDES")

			self.formV = [NURBS.surface(formV[0]), NURBS.surface(formV[1])]
		
		if len(formW) > 0:

			if len(formW) != 2:
				print("YOU NEED EXACTLY TWO W FORMING GUIDES")

			self.formW = [NURBS.surface(formW[0]), NURBS.surface(formW[1])]



	def createGrid(self):

		self.rails[0].extractGrid(self.dU, self.dV)
		self.rails[1].extractGrid(self.dU, self.dV)

		w_vecs = np.subtract(self.rails[1].grid[:], self.rails[0].grid[:])
		w_vecs = w_vecs[:,:]/self.dW

		self.grid = []

		for i in range(self.dW+1):

			self.grid.append(np.add(self.rails[0].grid[:,:], w_vecs[:,:]*i))

		self.grid = np.array(self.grid)
		self.grid = np.swapaxes(self.grid, 0, 2)


	def warpPts(self):

		if len(self.formV) != 0:

			if len(self.formV) != 2:
				print("YOU NEED EXACTLY TWO V FORMING GUIDES")
				return -1

			self.v_pts = [self.formV[0].extractGrid(self.dW, self.dV)]
			self.v_pts.append(self.formV[1].extractGrid(self.dW, self.dV))

		if len(self.formW) !=0:

			if len(self.formW) != 2:
				print("YOU NEED EXACTLY TWO W FORMING GUIDES")
				return -1

			self.w_pts = [self.formW[0].extractGrid(self.dW, self.dU)]
			self.w_pts.append(self.formW[1].extractGrid(self.dW, self.dU))


	def warpGrid(self, factor=1):

		if len(self.formV)>0:

			self.warpPts()

			# tan = [self.v_pts[0] - self.grid[0,:,:]]
			# layer_tangents = [tan]

			# for i in range(self.grid.shape[0]-1):

			# 	if i < self.grid.shape[0]/2 and i>0:
			# 		tan = (self.grid[i-1, :, :] - self.grid[i, :, :]) * abs(self.grid.shape[0]/2 - i)
			# 		layer_tangents.append(tan)
			# 	elif i!=0:
			# 		tan = (self.grid[i+1, :, :] - self.grid[i, :, :]) * abs(self.grid.shape[0]/2 - i)
			# 		layer_tangents.append(tan)

			# layer_tangents.append(self.v_pts[1] - self.grid[-1,:,:])

			# self.v_pts[0] = np.flip(self.v_pts[0],axis=0)
			self.v_pts[0] = np.flip(self.v_pts[0],axis=1)

			self.grid[0, :, :] = self.v_pts[0]
			self.grid[-1, :, :] = self.v_pts[1]


		# for i in range(self.grid.shape[0]):

		# 	self.grid[i, : ,:] = self.grid[i, :, :] + layer_tangents[i]*factor



	def tileUnit(self, unit):

		members = []

		for i in range(self.grid.shape[0] - 1):

			cell_layers = []

			for j in range(self.grid.shape[1] - 1):

				cell_cols = []

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

					# cell_cols.append(iG.boxMorph(crns, unit.pts))
					cell_cols.append(iG.boxMorph(crns))
					members.extend(cell_cols[-1].genDefault())

				cell_layers.append(cell_cols)

			self.cells.append(cell_layers)

		eP.exportPolyline(members, 'test_boxes.obj')



def Main():

	guide = ['NURBS/srf_01.stp', 'NURBS/srf_02.stp']
	formV = []
	formW = []
	formV = ['NURBS/srf_03.stp', 'NURBS/srf_04.stp']
	formW = ['NURBS/srf_05.stp', 'NURBS/srf_06.stp']

	lattice = conformalLattice(guide, formV, formW, 14, 6, 4)

	lattice.createGrid()

	lattice.warpGrid(1)

	lattice.tileUnit(np.array([[0,0,0]]))


Main()