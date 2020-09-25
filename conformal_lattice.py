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

		for i in range(len(formV)):
			self.formV.append(NURBS.surface(formV[i]))
		
		for i in range(len(formW)):
			self.formW.append(NURBS.surface(formW[i]))



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
			


	def warpGrid(self, factor=1):

		for i in range(len(self.formV)):

			self.v_pts.append(self.formV[i].extractGrid(self.dW, self.dV))

			if i==0:
				self.grid[0, :, :] = self.v_pts[i]
				print('FIRST U LAYER WARPED')
			if i==1:
				self.grid[-1, :, :] = self.v_pts[i]
				print('LAST U LAYER WARPED')


		for i in range(len(self.formW)):

			self.w_pts.append(self.formW[i].extractGrid(self.dW, self.dU))

			if i==0:
				self.grid[:, 0, :] = self.w_pts[0]
				print('FIRST V LAYER WARPED')
			if i==1:
				self.grid[:, -1, :] = self.w_pts[1]
				print('LAST V LAYER WARPED')



	def addMerges(self, travel_axi, travel_center, travel_range, merge_range = [1,1]):

		travel_st = int(travel_range[0]*self.grid.shape[travel_axi])
		travel_en = int(travel_range[1]*self.grid.shape[travel_axi])

		cnt = travel_center

		for i in range(self.grid.shape[travel_axi]):

			if i>=travel_st and i<=travel_en:

				if travel_axi == 0:

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



	def createCells(self, unit):

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

	lattice = conformalLattice(guide, formV, formW, 14, 8, 6)

	lattice.createGrid()

	lattice.warpGrid(1)

	# def addMerges(self, travel_axi, travel_center, travel_range, merge_axi, merge_range = 2):

	lattice.addMerges(1, [7, 3], [.33,.66], [3,2])

	lattice.createCells(np.array([[0,0,0]]))


Main()