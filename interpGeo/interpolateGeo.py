import numpy as np
import os
import sys
import math as m
from scipy.interpolate import CubicSpline


"""
interpCrv takes a series of points and forms separate equations for the x,y,z coordinates based on index numbers
parameters:
- pts: list of points for function generation
- reso: the number of points evaluated using the equations generated
"""

"""
interpSrf takes a series of contours (interpolated curves) and uses their points to form a grid that defines the UV coordinates of an interpolated srf
parameters:
- contours: interpolated curves
- resolution U: number of points between contours that define the surface grid
- resolution V: number of points along the contours that define the surface grid
- distDivide (optional, default = False): whether or not the surface grid is normalized based on distance 

"""

class interpCrv:

	def __init__ (self, pts, reso = 200):

		self.data = pts
		self.pts = pts
		self.tans = []
		self.reso = reso
		self.crvLen = self.updateLength()
		self.updateCrv()


	def updateCrv(self):

		indexes = np.arange(self.pts.shape[0])

		if np.linalg.norm(self.pts[0]-self.pts[-1]) == 0:

			self.fx = CubicSpline(indexes,self.pts[:,0],bc_type='periodic')
			self.fy = CubicSpline(indexes,self.pts[:,1],bc_type='periodic')
			self.fz = CubicSpline(indexes,self.pts[:,2],bc_type='periodic')

		else:

			self.fx = CubicSpline(indexes,self.pts[:,0])
			self.fy = CubicSpline(indexes,self.pts[:,1])
			self.fz = CubicSpline(indexes,self.pts[:,2])

		self.cntPoint()


	def updateLength(self):

		crvLen = 0

		for i in range(len(self.pts)-1):

			vec = np.subtract(np.array(self.pts[i]),np.array(self.pts[i+1]))
			crvLen = crvLen + np.linalg.norm(vec)
		
		self.crvLen = crvLen

		return self.crvLen

	# updateTans() recalculates the tangents derivatives at the spline control points

	def updateTans(self):

		#self.updateCrv()

		indexes = np.arange(self.pts.shape[0])

		tans = []
		gaps = []

		for i in range(len(indexes)):

			tan = [self.fx(indexes[i],1),self.fy(indexes[i],1),self.fz(indexes[i],1)]
			mag = np.linalg.norm(np.array(tan))
			if mag == 0 and i < len(indexes)-1:
				tan = np.subtract(self.pts[i + 1], self.pts[i])
			elif mag ==0 and i == len(indexes)-1:
				tan = np.subtract(self.pts[i], self.pts[i-1])
			else:
				tan = [tan[0]/mag, tan[1]/mag, tan[2]/mag]

			tans.append(tan)

		self.tans = np.array(tans)


	# genPts(param) function generates the curve based a number of points (reso). Each dimension gets its own spline equation [f(x),f(y),f(z)] 


	def genPts(self,reso):

		step = (len(self.data)-1)/reso
		indexes = np.arange(0,step*reso+step,step)

		pts = []
		tans = []
		gaps = []

		for i in range(len(indexes)):

			pts.append([self.fx(indexes[i]),self.fy(indexes[i]),self.fz(indexes[i])])
			tan = [self.fx(indexes[i],1),self.fy(indexes[i],1),self.fz(indexes[i],1)]
			mag = np.linalg.norm(np.array(tan))
			if mag == 0:
				gaps.append(i)
			tans.append([tan[0]/mag, tan[1]/mag, tan[2]/mag])

		for i in range(len(gaps)):

			if gaps[i] < len(pts)-1:
				vec = np.subtract(pts[gaps[i]+1], pts[gaps[i]])
			else:
				vec = np.subtract(pts[gaps[i]], pts[gaps[i]-1])

			mag = np.linalg.norm(vec)
			tans[gaps[i]] = list(vec/mag)

		self.pts = np.array(pts)
		self.tans = np.array(tans)

		return pts


	# genUniform(param) function generates the curve based a number of points (reso) that are equidistant from one another on the overall curve (approximately)
	# before this function can be called the self.pts array must be created.
	# parameters:
	# param = value from 0-1 determining where along the curve you want to evaluate a point

	def genUniform(self,segments):

		pts = []

		print(segments)

		for i in range(segments+1):

			param = self.evalDistParam(i/segments)
			pts.append(self.evalParam(param))

		self.pts = np.array(pts)

		self.updateTans()

		return pts


	# evalDistParam(normParam) function gets point closest to normalized parameter on curve
	# before this function can be called the self.pts array must be created.
	# parameters:
	# param = value from 0-1 determining where along the curve you want to evaluate a point


	def evalDistParam(self,normParam):

		self.updateLength()

		dist = self.crvLen*normParam

		param = 0

		total = 0

		for i in range(len(self.pts)-1):

			vec = np.subtract(self.pts[i],self.pts[i+1])

			total = np.linalg.norm(vec) + total
			
			if total > dist:

				diff = total - dist
				mag = np.linalg.norm(vec)
				factor = diff/mag
				if factor>1:
					print('shit')
				param = (i + factor)/len(self.pts)
				break

			elif total == dist:

				param = i/len(self.pts)
				break


		return param

	#evalParam points retrieves a point on the spline based on the curves existing x y and z functions. If you want to evaluate based on a unitized
	#see genDistPts

	def evalParam(self,param):

		index = param*(len(self.data)-1)

		result = [self.fx(index),self.fy(index),self.fz(index)]

		return result

	#cntPoint finds the approximate center point of the curve based on the cntrl points

	def cntPoint(self):

		points = self.genPts(self.reso)
		sum = [0,0,0]

		for i in range(len(points)):

			for j in range(len(points[i])):

				sum[j] = sum[j]+points[i][j]

		self.cnt = np.array(sum)*1/len(points)

		return self.cnt

	#closestPt(pt) finds the approx location and index of the closest point on the curve to a point given by the user
	#pt = test point to find closest point to.

	def closestPt(self,pt):

		points = self.genPts(self.reso)
		index = 0
		minimum = np.linalg.norm(np.array(points[0]) - pt)

		#print(np.full((self.pts.shape),pt))

		for i in range(len(points)):

			dist = np.linalg.norm(np.array(points[i]) - pt)
			if dist<minimum:
				minimum = dist
				index = i

		return [index,np.array(points[index])]

	#recordCrv(dir) writes the point information that defines the curve to a .csv file at a specified directory
	#dir = output directory

	def recordCrv(self,dir):

		#self.genPts(200)
		f = open(dir,'w')

		for i in range(len(self.pts)):
			for j in range(len(self.pts[i])):
				f.write(str(self.pts[i][j]))
				if j<len(self.pts[i])-1:
					f.write(',')
			if i<len(self.pts)-1:
				f.write('\n')

		f.close()

	#mvCrv(vec) moves the entire curve by a user specified vector
	#vec = vector you want to move the entire curve by

	def mvCrv(self,vec):

		for i in range(len(self.pts)):
			self.pts[i] = np.array([self.pts[i][0]+vec[0],self.pts[i][1]+vec[1],self.pts[i][2]+vec[2]])


	#rotCrv(ang,origin,axi) rotates an entire curve around an origin point and axis by a specified angle
	#ang = angle to rotate curve by
	#origin = center of rotation
	#axi = axis you are rotating the curve around


	def rotCrv(self,ang,origin,axi):
		
		for i in range(len(self.pts)):
			vec = np.subtract(self.pts[i],origin)
			vec = np.array(vecRotate(vec,ang,axi))
			self.pts[i] = np.add(vec,origin)


	#offsetCrv(plane,value) offsets a curve in a specified plane by a certain value
	#plane = a numpy vector defining the plane normal
	#value = the magnitude the curve will be offset by (can be - or + to offset in and out)


	def offsetCrv(self,plane,value):

		oldVec = np.cross(self.tans[0],plane)
		oldVec = oldVec/np.linalg.norm(oldVec)

		for i in range(len(self.pts)):

			vec = np.cross(self.tans[i],plane)
			vec = vec/np.linalg.norm(vec)

			if np.dot(vec,oldVec)<0:
				factor = -value
			else:
				factor = value
			
			self.pts[i] = np.add(self.pts[i],vec*factor)
			oldVec = vec

		return self.pts


	#scaleCrv(origin,factor) scales the entire curve from an origin point by a specified factor in the [x,y,z] directions
	#origin = origin of scale
	#factor = list of three values defining how much the curve will be scaled in the x,y,z directions


	def scaleCrv(self,origin,factor):

		vecs = []

		for i in range(len(self.pts)):

			vecs.append([self.pts[i][0]-origin[0],self.pts[i][1]-origin[1],self.pts[i][2]-origin[2]])
			vecs[-1] = [vecs[-1][0]*factor[0], vecs[-1][1]*factor[1], vecs[-1][2]*factor[2]]

			self.pts[i] = np.add(origin,vecs[i])

		return self.pts


	#reverseCrv() reverses the direction of the curve (start point and end point)


	def reverseCrv(self):

		newPts = []

		for i in range(len(self.pts)):

			newPts.append(self.pts[len(self.pts)-1-i])

		self.pts = newPts




"""
interpSrf takes a series of contours (interpolated curves) and uses their points to form a grid that defines the UV coordinates of an interpolated srf
parameters:
- contours: interpolated curves
- resolution U: number of points between contours that define the surface grid
- resolution V: number of points along the contours that define the surface grid
- distDivide (optional, default = False): whether or not the surface grid is normalized based on distance 

"""

class interpSrf:


	def __init__(self,contours,resolutionU,resolutionV,distDivide = False):

		self.vCrv = contours
		self.uCrv = []
		self.resoU = resolutionU
		self.resoV = resolutionV
		self.divByDist = distDivide

		self.U = []
		self.V = []
		self.faces = []
		self.v = []
		self.vn = []
		self.vt = []
		self.utan = []
		self.vtan = []

		self.box_sets = []

		crosses = []
		columns = []

		# START: divides the original contours into a a set of resolution V points

		# step = self.vCrv[0].crvLen/self.resoV

		for i in range(len(self.vCrv)):

			if self.divByDist == True:

				columns.append(self.vCrv[i].genUniform(self.resoV))

			else:

				columns.append(self.vCrv[i].genPts(self.resoV-1))

		rowPts = []

		

		for i in range(len(columns[0])):

			cross = []

			for j in range(len(columns)):

				cross.append(columns[j][i])

			crosses.append(interpCrv(np.array(cross)))
			rowPts.append(cross)

		# END


		# START: divides the original contours into a a set of resolution U points


		for i in range(self.resoV):

			#step = crosses[i].crvLen/self.resoU

			row = crosses[i].genPts(self.resoU)

			#row = crosses[i].genUniform(step)

			self.U.append(row)

		# END

		self.genVCoord()
		self.boundingBox()


	# updates the V coordinates 


	def genVCoord(self):

		self.V = []

		for i in range(len(self.U[0])):

			row = []

			for j in range(len(self.U)):

				row.append(self.U[j][i])

			self.V.append(row)

	# calculates the bounding box around a interpolated surface. [lower corners 1-4 counterclockwise then upper corners 1-4 counterclockwise]

	def boundingBox(self):


		minX,minY,minZ = self.U[0][0][0],self.U[0][0][1],self.U[0][0][2]
		maxX,maxY,maxZ = self.U[0][0][0],self.U[0][0][1],self.U[0][0][2]

		for i in range(len(self.U)):
			for j in range(len(self.U[i])):

				pt = self.U[i][j]

				if pt[0]<minX: minX = pt[0]
				if pt[0]>maxX: maxX = pt[0]
				if pt[1]<minY: minY = pt[1]
				if pt[1]>maxY: maxY = pt[1]
				if pt[2]<minZ: minZ = pt[2]
				if pt[2]>maxZ: maxZ = pt[2]

		self.bounds = [[minX,minY,minZ],[maxX,minY,minZ],[maxX,maxY,minZ],[minX,maxY,minZ]]
		self.bounds.extend([[minX,minY,maxZ],[maxX,minY,maxZ],[maxX,maxY,maxZ],[minX,maxY,maxZ]])


	# correctAlign(thres,level,align): Aligns all the surface verticies within a specified tolerance to a defined plane in the X,Y or Z direction
	# thres = tolerance between vertex location and the level that you want to align
	# level = the X,Y,Z value at which vertices will be corrected
	# align =  X Y or Z (what plane you want to align to)

	def correctAlign(self,thres,level,align):

		correctedLocs = []

		for i in range(len(self.U)):
			for j in range(len(self.U[i])):

				if abs(self.U[i][j][align] - level)<thres:
					
					self.U[i][j][align] = self.bounds[0][align]

		self.genVCoord()
		self.boundingBox()

	def swapUV(self):

		self.U = self.V

		resoU = self.resoU
		self.resoU = self.resoV
		self.resoV = resoU

		self.genVCoord()
		self.boundingBox()


	def recordPts(self):

		f = open('gridPts.csv','w')

		for i in range(len(self.U)):
			for j in range(len(self.U[i])):

				f.write(str(self.U[i][j][0]))
				f.write(',')
				f.write(str(self.U[i][j][1]))
				f.write(',')
				f.write(str(self.U[i][j][2]))
				f.write('\n')

		f.close()

	# createFaces() generates the faces, vertex normals and vertex texture values that would be used to export a mesh 

	def createFaces(self):

		for i in range(len(self.U)):

			for j in range(len(self.U[i])):
				
				self.v.append(self.U[i][j])

				if i<len(self.U)-1:
					vecX  = np.subtract(self.U[i+1][j],self.U[i][j])
				else:
					vecX = np.subtract(self.U[i][j],self.U[i-1][j])

				if j<len(self.U[i])-1:
					vecY = np.subtract(self.U[i][j+1],self.U[i][j])
				else:
					vecY = np.subtract(self.U[i][j],self.U[i][j-1])

				
				vNorm = np.cross(vecX,vecY)

				self.vn.append(np.multiply(vNorm,1/np.linalg.norm(vNorm)))

				self.vt.append(np.array([i/(len(self.U)-1),j/(len(self.U[i])-1)]))

				self.utan.append(vecX/np.linalg.norm(vecX))

				self.vtan.append(vecY/np.linalg.norm(vecY))


		for i in range(len(self.U)-1):
			for j in range(len(self.U[i])):

				if j>0:

					face = []
					face.append(i*len(self.U[i])+j)
					face.append(i*len(self.U[i])+j+1)
					face.append((i+1)*len(self.U[i])+j+1)
					face.append((i+1)*len(self.U[i])+j)

					self.faces.append(face)

	# exportMesh(output) uses a standard .obj format to export the interpsrf as a quad mesh
	# output = the output .obj file

	def exportMesh(self,output):

		self.createFaces()

		out = open(output,'w')

		out.write('# Rhino')
		out.write('\n')
		out.write('\n')
		out.write('o object_1')
		out.write('\n')

		for i in range(len(self.v)):

			out.write('v ')
			out.write(str(self.v[i][0]))
			out.write(' ')
			out.write(str(self.v[i][1]))
			out.write(' ')
			out.write(str(self.v[i][2]))
			out.write('\n')

		for i in range(len(self.v)):

			out.write('vt ')
			out.write(str(self.vt[i][0]))
			out.write(' ')
			out.write(str(self.vt[i][1]))
			out.write('\n')

		for i in range(len(self.v)):

			out.write('vn ')
			out.write(str(self.vn[i][0]))
			out.write(' ')
			out.write(str(self.vn[i][1]))
			out.write(' ')
			out.write(str(self.vn[i][2]))
			out.write('\n')

		for i in range(len(self.faces)):

			out.write('f ')
			out.write(str(self.faces[i][0]))
			out.write(' ')
			out.write(str(self.faces[i][1]))
			out.write(' ')
			out.write(str(self.faces[i][2]))
			out.write(' ')
			out.write(str(self.faces[i][3]))
			out.write('\n')

		out.close()

		print("EXPORTED!!!")


	# offsetSrf(thickness) offsets the U,V grid based on the vertex normals and a specified thickness that can be (- or +)


	def offsetSrf(self,thickness,solid=True):

		self.createFaces()

		U,V = [],[]

		srfCrvs = []
		faces = []
		norms = []
		tans = []
		vertNum = len(self.v)

		for i in range(len(self.U)):
			
			row = []

			for j in range(len(self.U[i])):

				index = i*len(self.U[i])+j

				pt = np.add(self.U[i][j],list(self.vn[index]*thickness))

				if solid:

					self.v.append(pt)
					self.vn.append(list(np.multiply(self.vn[index],-1)))
					self.vt.append(self.vt[index])

					if i==0 and j<len(self.U[i])-1:

						faces.append([index+1,index+2,index+2+vertNum,index+vertNum+1])

					if i==len(self.U)-1 and j<len(self.U[i])-1:

						faces.append([index+1,index+2,index+2+vertNum,index+vertNum+1])

					if j==0 and i<len(self.U)-1:

						faces.append([index+1,index+len(self.U[i])+1,index+len(self.U[i])+vertNum+1,index+vertNum+1])

					if j==len(self.U[i])-1 and i<len(self.U)-1:

						faces.append([index+1,index+len(self.U[i])+1,index+len(self.U[i])+vertNum+1,index+vertNum+1])

				else:

					self.v[index] = pt
					self.vn[index] = list(np.multiply(self.vn[index],-1))
					self.vt[index] = self.vt[index]


				row.append(pt)

			U.append(row)
			srfCrvs.append(interpCrv(np.array(row)))

		if solid:

			faceNum = len(self.faces)

			for i in range(faceNum):

				newFace = []

				for j in range(len(self.faces[i])):

					newFace.append(self.faces[i][j]+vertNum)

				if newFace not in self.faces:

					self.faces.append(newFace)

			for i in range(len(faces)):

				if faces[i] not in self.faces:

					self.faces.append(faces[i])

		temp = []

		for i in range(len(self.faces)):

			if self.faces[i] not in temp:

				temp.append(self.faces[i])
		
		self.faces = temp


		return srfCrvs


	def offsetVerts(self,thickness, dim3 = True):

		midU = int(len(self.U)/2)

		baseVerts = len(self.U)*len(self.U[0])

		for i in range(len(self.U)):

			midV = int(len(self.U[i])/2)

			for j in range(len(self.U[i])):

				index = i*len(self.U[i])+j
				index02 = index + baseVerts

				if dim3:

					self.v[index] = list(np.add(self.v[index],np.multiply(self.vn[index],thickness)))
					self.v[index02] = list(np.add(self.v[index02],np.multiply(self.vn[index02],thickness)))

				if j<midV:

					vec = np.subtract(self.U[i][j+1],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(j-midV)/midV*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				if j>midV:

					vec = np.subtract(self.U[i][j-1],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(midV-j)/midV*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				if i<midU:

					vec = np.subtract(self.U[i+1][j],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(j-midU)/midU*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				if i>midU:

					vec = np.subtract(self.U[i-1][j],self.U[i][j])
					vec = vec/np.linalg.norm(vec)
					vec = vec*(j-midU)/midU*thickness
					self.v[index] = list(np.add(self.v[index],vec))

					if dim3:
						self.v[index02] = list(np.add(self.v[index02],vec))

				#self.U[i][j] = self.v[index]

		return self.v


	def createBoxes(self,atts,hAtts,dimU,dimV,dimW,height,limit=10,limitH=15,closed=True, merge_V_centers = [], merge_V_ranges = [], merge_V = [], merge_U_centers = [], merge_U_ranges = [], merge_U = []):

		u = []
		v = []
		w = []

		columns = []
		cols = []
		rows = []
		crosses = []

		layers = []

		for i in range(len(self.U)):

			columns.append(interpCrv(np.array(self.U[i])))
			cols.append(columns[-1].genPts(dimV))

		for i in range(len(cols[0])):
			cross = []
			for j in range(len(cols)):
				cross.append(cols[j][i])
			crosses.append(interpCrv(np.array(cross)))

		for i in range(len(crosses)):
			row = crosses[i].genPts(dimV)
			rows.append(row)


		for n in range(dimW):
			
			nxtRows = []

			for i in range(len(cols)):
				
				nxtRow = []
				prevN = np.array([0,0,0])

				for j in range(len(cols[i])):

					if j==len(cols[i])-1:
						vec01 = np.subtract(cols[i][j],cols[i][j-1])
					else:
						vec01 = np.subtract(cols[i][(j+1)%len(cols[i])],cols[i][j])

					if i==len(cols)-1:
						vec02 = np.subtract(cols[i][j],cols[i-1][j])
					else:
						vec02 = np.subtract(cols[(i+1)%len(cols)][j],cols[i][j])
					
					norm = np.cross(vec01,vec02)

					if i==len(cols)-1 or i==0:
						norm[2] = 0

					if j!=0 and np.dot(prevN,norm)<-.7:
						norm = -norm

					prevN = norm/np.linalg.norm(norm)
					norm = norm/np.linalg.norm(norm)

					if n>1:
						rise = 2.5*height
					else:
						rise = height

					nxtRow.append(np.add(cols[i][j],norm*n*rise/(dimW-1)))

				nxtRows.append(nxtRow)

			layers.append(nxtRows)


		cells = []
		upperCells = []
		finalCells = []
		cellAreas = []
		crvs = []
		num = 0

		if closed:
			closeAdjustment = 0
		else:
			closeAdjustment = -1
		

		for i in range(len(layers)):
			for j in range(len(layers[i])+closeAdjustment):
				for k in range(len(layers[i][j])-1):

					pt01 = layers[i][j][k]
					pt02 = layers[i][j][k+1]
					pt03 = layers[i][(j+1)%len(layers[i])][k+1]
					pt04 = layers[i][(j+1)%len(layers[i])][k]
					botPts = [pt01,pt02,pt03,pt04]

					for n in range(len(botPts)):
						if m.isnan(botPts[n][0]) or m.isnan(botPts[n][1]) or m.isnan(botPts[n][2]):
							botPts[n] = np.add(botPts[n-1],botPts[(n+1)%len(botPts)])
							botPts[n] = list(botPts[n]/2)
							print('holy shit bottom')

					pt05 = layers[(i+1)%len(layers)][j][k]
					pt06 = layers[(i+1)%len(layers)][j][k+1]
					pt07 = layers[(i+1)%len(layers)][(j+1)%len(layers[i])][k+1]
					pt08 = layers[(i+1)%len(layers)][(j+1)%len(layers[i])][k]
					topPts = [pt05,pt06,pt07,pt08]

					for n in range(len(topPts)):
						if m.isnan(topPts[n][0]) or m.isnan(topPts[n][1]) or m.isnan(topPts[n][2]):
							topPts[n] = np.add(topPts[n-1],topPts[(n+1)%len(topPts)])
							topPts[n] = list(topPts[n]/2)
							print('holy shit top')

					cellAreas.append(np.linalg.norm(np.subtract(pt01,pt02))*np.linalg.norm(np.subtract(pt04,pt01)))

					newBotPts = []
					newTopPts = []
					rises = []

					for p in range(len(botPts)):

						pt = botPts[p]
						factor = 0
						count = 0

						for q in range(len(hAtts)):

							attPt = hAtts[q].closestPt(pt)[1]
							vec = np.subtract(attPt,pt)

							if np.linalg.norm(vec)<limitH:
								factor = factor + 1-np.linalg.norm(vec)/limitH
								count = count + 1

						if count>0:
							factor = factor/count

						riseNorm = np.subtract(topPts[p],pt)

						rises.append(riseNorm*factor*2)
						newBotPts.append(np.add(rises[-1],pt))
						newTopPts.append(np.add(rises[-1],topPts[p]))


					botPts = newBotPts
					topPts = newTopPts


					if i>0:

						newTopPts = []
						valid = False
						edges = []
						rises = []

						for p in range(len(botPts)):

							factor = 0
							count = 0

							pt = botPts[p]

							for q in range(len(atts)):
								attPt = atts[q].closestPt(pt)[1]
								vec = np.subtract(attPt,pt)

								if np.linalg.norm(vec)<limit:
									factor = factor + np.linalg.norm(vec)/limit
									count = count + 1
									valid = True

							if count>0:
								factor = 1 - factor/count

							
							edges.append(np.subtract(botPts[(p+1)%len(botPts)],pt)*.5)
							rises.append(np.subtract(topPts[p],pt)*factor)
							newTopPts.append(np.add(rises[-1],pt))


						if valid:

							cell = []
							cell.extend(botPts)
							cell.extend(newTopPts)
							crvs.append(cell)
							#finalCells.append(cell)
							upperCells.append(cell)

					else:

						cell = []
						cell.extend(botPts)
						cell.extend(topPts)
						crvs.append(cell)
						finalCells.append(cell)

					print(str(int(num/((len(layers))*len(layers[i])*(len(layers[i][j])))*100*100)/100) + ' %')
					num = num + 1

		self.finalCells = finalCells
		self.cellMembers = crvs
		self.upperCells = upperCells
		self.cellU = dimU
		self.cellV = dimV


	def divideReso(self,domU,domV,thresX,thresY,loopLimit=4):

		for n in range(loopLimit):

			allSubCells = []
			maxCellX = 0
			maxCellY = 0

			for i in range(len(self.finalCells)):

				midB , midT , cntB , cntT= [] , [] , np.array([0,0,0]) , np.array([0,0,0])

				for j in range(int(len(self.finalCells[i])/2)):

					midB.append(np.add(self.finalCells[i][j],self.finalCells[i][(j+1)%4])/2)
					cntB = cntB + midB[-1]
					index = j+5
					if j == 3:
						index = 4
					midT.append(np.add(self.finalCells[i][j+4],self.finalCells[i][index])/2)
					cntT = cntT + midT[-1]

				distX = np.linalg.norm(np.subtract(self.finalCells[i][4],self.finalCells[i][7]))
				distY = np.linalg.norm(np.subtract(self.finalCells[i][4],self.finalCells[i][5]))

				cntT = cntT/4
				cntB = cntB/4

				divideX = False
				divideY = False


				if i%self.cellV<self.cellV*domV[1] and i%self.cellV>=self.cellV*domV[0] and n==0:

					if int(i/self.cellV)<self.cellU*domU[1] and int(i/self.cellV)>self.cellU*domU[0]:

						divideX = True
						divideY = True

				if distX > thresX: 
					divideX = True
				if distY > thresY: 
					divideY = True

				if distX > maxCellX: 
					maxCellX = distX
				if distY > maxCellY: 
					maxCellY = distY


				if divideX and divideY:

					subCells = []

					bot = [self.finalCells[i][0],midB[0],cntB,midB[3]]
					top = [self.finalCells[i][4],midT[0],cntT,midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.finalCells[i][1],midB[1],cntB]
					top = [midT[0],self.finalCells[i][5],midT[1],cntT]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [cntB,midB[1],self.finalCells[i][2],midB[2]]
					top = [cntT,midT[1],self.finalCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[3],cntB,midB[2],self.finalCells[i][3]]
					top = [midT[3],cntT,midT[2],self.finalCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX and divideY==False:

					subCells = []

					bot = [self.finalCells[i][0],self.finalCells[i][1],midB[1],midB[3]]
					top = [self.finalCells[i][4],self.finalCells[i][5],midT[1],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[1],self.finalCells[i][2],self.finalCells[i][3],midB[3]]
					top = [midT[1],self.finalCells[i][6],self.finalCells[i][7],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)


				if divideY and divideX==False:

					subCells = []

					bot = [self.finalCells[i][0],midB[0],midB[2],self.finalCells[i][3]]
					top = [self.finalCells[i][4],midT[0],midT[2],self.finalCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.finalCells[i][1],self.finalCells[i][2],midB[2]]
					top = [midT[0],self.finalCells[i][5],self.finalCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX==False and divideY==False:

					allSubCells.append(self.finalCells[i])

			self.finalCells = allSubCells

			allSubCells = []

			print(maxCellX)
			print(maxCellY)

			if maxCellX<thresX and maxCellY<thresY: 
				break

		return self.finalCells


	def divideResoUpper(self,thresX,thresY,loopLimit=4):

		for n in range(loopLimit):

			allSubCells = []
			maxCellX = 0
			maxCellY = 0

			for i in range(len(self.upperCells)):

				midB , midT , cntB , cntT= [] , [] , np.array([0,0,0]) , np.array([0,0,0])

				for j in range(int(len(self.upperCells[i])/2)):

					midB.append(np.add(self.upperCells[i][j],self.upperCells[i][(j+1)%4])/2)
					cntB = cntB + midB[-1]
					index = j+5
					if j == 3:
						index = 4
					midT.append(np.add(self.upperCells[i][j+4],self.upperCells[i][index])/2)
					cntT = cntT + midT[-1]

				distX = np.linalg.norm(np.subtract(self.upperCells[i][4],self.upperCells[i][7]))
				distY = np.linalg.norm(np.subtract(self.upperCells[i][4],self.upperCells[i][5]))

				cntT = cntT/4
				cntB = cntB/4

				divideX = False
				divideY = False

				if distX > thresX: 
					divideX = True
				if distY > thresY: 
					divideY = True

				if distX > maxCellX: 
					maxCellX = distX
				if distY > maxCellY: 
					maxCellY = distY


				if divideX and divideY:

					subCells = []

					bot = [self.upperCells[i][0],midB[0],cntB,midB[3]]
					top = [self.upperCells[i][4],midT[0],cntT,midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.upperCells[i][1],midB[1],cntB]
					top = [midT[0],self.upperCells[i][5],midT[1],cntT]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [cntB,midB[1],self.upperCells[i][2],midB[2]]
					top = [cntT,midT[1],self.upperCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[3],cntB,midB[2],self.upperCells[i][3]]
					top = [midT[3],cntT,midT[2],self.upperCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX and divideY==False:

					subCells = []

					bot = [self.upperCells[i][0],self.upperCells[i][1],midB[1],midB[3]]
					top = [self.upperCells[i][4],self.upperCells[i][5],midT[1],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[1],self.upperCells[i][2],self.upperCells[i][3],midB[3]]
					top = [midT[1],self.upperCells[i][6],self.upperCells[i][7],midT[3]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)


				if divideY and divideX==False:

					subCells = []

					bot = [self.upperCells[i][0],midB[0],midB[2],self.upperCells[i][3]]
					top = [self.upperCells[i][4],midT[0],midT[2],self.upperCells[i][7]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					bot = [midB[0],self.upperCells[i][1],self.upperCells[i][2],midB[2]]
					top = [midT[0],self.upperCells[i][5],self.upperCells[i][6],midT[2]]
					cell = bot
					cell.extend(top)
					subCells.append(cell)

					allSubCells.extend(subCells)

				if divideX==False and divideY==False:

					allSubCells.append(self.upperCells[i])

			self.upperCells = allSubCells

			allSubCells = []

			print(maxCellX)
			print(maxCellY)

			if maxCellX<thresX and maxCellY<thresY: 
				break

		return self.upperCells


	def exportJSON(self,dir):

		crvs = self.cellMembers
		finalCells = self.finalCells
		upperCells = self.upperCells

		f = open('testCrvs.csv','w')

		for i in range(len(crvs)):

			for j in range(len(crvs[i])):

				for k in range(len(crvs[i][j])):

					f.write(str(crvs[i][j][k]))
					if k<len(crvs[i][j])-1:
						f.write(',')

				if j<len(crvs[i])-1:
					f.write(' ')

			if i<len(crvs)-1:
				f.write('\n')

		f.close()

		f = open(dir,'w')

		f.write('[')

		for i in range(len(finalCells)):

			f.write('[')

			if finalCells[i] != None:

				for j in range(len(finalCells[i])):

					f.write('[')

					for k in range(len(finalCells[i][j])):

						f.write(str(finalCells[i][j][k]))
						if k<len(finalCells[i][j])-1:
							f.write(',')
						else:
							f.write(']')

					if j<len(finalCells[i])-1:
						f.write(',')
					else:
						f.write(']')

			if i<len(finalCells)-1:
				f.write(',')
			else:
				f.write(']')

		f.close()

		########### CREATE UPPER BOXES #############

		f = open(dir.split('.')[0] + '_h.json','w')
		f.write('[')

		for i in range(len(upperCells)):
			f.write('[')
			for j in range(len(upperCells[i])):
				f.write('[')
				for k in range(len(upperCells[i][j])):
					f.write(str(upperCells[i][j][k]))
					if k<len(upperCells[i][j])-1:
						f.write(',')
					else:
						f.write(']')

				if j<len(upperCells[i])-1:
					f.write(',')
				else:
					f.write(']')

			if i<len(upperCells)-1:
				f.write(',')
			else:
				f.write(']')
		f.close()

		return [finalCells,upperCells]



"""
interpPipe takes a curve (interpolated curve) and a starting radius to create a possibly tapered pipe with caps
- path: interpolated curve
- radius: starting radius of the pipe
"""

class interpPipe:

	def __init__ (self, path ,radius):

		self.pts = path.pts
		self.t = path.tans
		self.sections =  []
		self.radius = radius

	def genSections (self, taper = True, gradient = .5, start = .1):

		self.radii = []

		for i in range(len(self.pts)):

			if taper and i/len(self.pts)>start:

				f = 1 - .5*m.pow(i/len(self.pts)-start,gradient)
			
			else:
			
				f = 1

			if f>1: f=1
			
			self.radii.append(self.radius*f)

			self.sections.append(genCircleSection(self.pts[i],self.t[i],radius))

	def genPipe(self, reso = 8, cap = True):

		if len(self.sections) == 0:

			self.genSections()

		self.srf = interpSrf(self.sections,len(self.pts),reso)
		self.srf.createFaces()

		if cap:
			self.srf.v[self.U[-1]]





def importCrvsFromCSV(loc,reso):

	f = open(loc,'r')
	lines = f.read().split('\n')

	crvs = []

	for i in range(len(lines)):
		entries = lines[i].split(' ')
		crv = []
		for j in range(len(entries)):
			element = entries[j].split(',')
			pt = [float(element[0]),float(element[1]),float(element[2])]
			crv.append(pt)

		crvs.append(interpCrv(np.array(crv),reso))

	return crvs





def importCrvOBJ(loc):

	f = open(loc,'r')
	allLines = f.read().split('\n')

	crvs = []
	verts = []

	broken = False

	for i in range(len(allLines)):

		line = allLines[i].split(' ')

		if line[0] == 'v':

			if '\\' in line:

				pt = [float(line[1]),float(line[2])]
				pt.append(float(allLines[i+1]))

				broken = True

			else:

				pt = [float(line[1]),float(line[2]),float(line[3])]
				verts.append(pt)
				broken = False

		elif len(verts)>0 and broken == False:

			crvs.append(verts)
			verts = []

	return [np.array(crvs),np.array(verts)]


####
# importMeshOBJ (loc)
# imports .obj as a list of faces (quads and triangles) and a list of vertices (3D points)
# loc (str) = the path location and file name for the output .obj file
####

def importMeshOBJ(loc):

	f = open(loc,'r')
	allLines = f.read().split('\n')

	faces = []
	verts = []

	broken = False

	for i in range(len(allLines)):

		line = allLines[i].split(' ')

		if line[0] == 'v':

			pt = [float(line[1]),float(line[2]),float(line[3])]
			verts.append(pt)
			broken = False

		elif line[0] == 'f':

			face = []

			for j in range(len(line)):

				if j>0:

					f_traits = line[j].split('/')
					
					for k in range(len(f_traits)):

						if f_traits[k] != '' and k == 0: 
							face.append(int(f_traits[k]) - 1)

			if len(face) >= 3:
				faces.append(face)


	return [np.array(verts),faces]


####
# exportMeshOBJ (faces, verts, norms, loc)
# exports mesh as .obj file
# faces (tri or quad) = a list of triangle or quads represented by lists of 3 or 4 index values per face
# verts (numpy array 3D points) = an array of vertex points  
# loc (str) = the path location and file name for the output .obj file
####

def exportMeshOBJ(faces, verts, norms, loc):

	f = open(loc,'w')

	f.write('# Rhino')
	f.write('\n')
	f.write('\n')

	for i in range(verts.shape[0]):

		f.write('v')
		f.write(' ')
		f.write(str(verts[i, 0]))
		f.write(' ')
		f.write(str(verts[i, 1]))
		f.write(' ')
		f.write(str(verts[i, 2]))
		f.write('\n')

	for i in range(norms.shape[0]):

		f.write('vn')
		f.write(' ')
		f.write(str(norms[i, 0]))
		f.write(' ')
		f.write(str(norms[i, 1]))
		f.write(' ')
		f.write(str(norms[i, 2]))
		f.write('\n')

	for i in range(len(faces)):

		f.write('f')
		f.write(' ')

		for j in range(len(faces[i])):

			f.write(str(int(faces[i][j] + 1)))
			f.write('//')
			f.write(str(int(faces[i][j] + 1)))

			if j < len(faces[i])-1:
				f.write(' ')

		if i < len(faces) - 1:
			f.write('\n')

	f.close()

	return [np.array(verts),faces]



"""
BOXMORPH CLASS creates a mapping object for transforming 3D points from a unit box to a user defined box
Input:
crns (numpy array of (8,3) points) = points defining the corners of the morphing box
src (numpy array of 3D points) = the points placed in world coordinates defining the vertices to be mapped
faces (numpy array of 3 or 4 indicies) = a list defining the faces of the input and output meshes
"""

class boxMorph:

	def __init__(self, crns, src = np.array([]), faces = np.array([])):

		self.crns = crns

		# if there were no faces or vertices provided the box maps the corners of a normal 1x1x1 box

		if src.shape[0] == 0: 
			self.src = self.crns
		else:
			self.src = src

		# gets the general width and length of the box based on the starting corner
		self.W = np.linalg.norm(self.crns[1] - self.crns[0])
		self.L = np.linalg.norm(self.crns[3] - self.crns[0])

		# defines all the x_axis edges in the box
		self.x_ax = [self.crns[1] - self.crns[0]]
		self.x_ax.append(self.crns[2] - self.crns[3])
		self.x_ax.append(self.crns[5] - self.crns[4])
		self.x_ax.append(self.crns[6] - self.crns[7])

		self.lock = False

		self.x_ax = np.array(self.x_ax)
		self.v = self.mapV(self.src)  # defines the mapped vertices
		self.f = faces


	####
	# mapV (faces, verts, norms, loc)
	# exports mesh as .obj file
	# faces (tri or quad) = a list of triangle or quads represented by lists of 3 or 4 index values per face
	# verts (numpy array 3D points) = an array of vertex points  
	# loc (str) = the path location and file name for the output .obj file
	####

	def mapV(self, pts):

		mapPts = []
		minX = np.min(pts[:,0])
		minY = np.min(pts[:,1])
		minZ = np.min(pts[:,2])

		maxX = np.max(pts[:,0])
		maxY = np.max(pts[:,1])
		maxZ = np.max(pts[:,2])

		nrmPt = [0,0,0]

		# if the dimension between points is 0 in any direction ignore all the vertices (no 2D objects)

		if maxX-minX != 0 and maxY-minY != 0 and maxZ-minZ != 0:

			# SCROLLS THROUGH EACH POINT AND MAPS IT TO THE BOX BASED ON ITS PARAMETER LOCATION IN SPACE

			for i in range(pts.shape[0]):

				# normalize the vertex points

				nrmPt[0] = (pts[i,0] - minX)/(maxX-minX)
				nrmPt[1] = (pts[i,1] - minY)/(maxY-minY)
				nrmPt[2] = (pts[i,2] - minZ)/(maxZ-minZ)

				# find the point location on the bottom face of the target box

				st_botPt = nrmPt[0]*self.x_ax[0] + self.crns[0] # bottom first x axis evaluate location
				en_botPt = nrmPt[0]*self.x_ax[1] + self.crns[3] # bottom second x axis evaluate location

				# draw line between previous points and evaluate location on y axis using V coordinate
				bot_pt = st_botPt + (en_botPt - st_botPt)*nrmPt[1]

				# repeat process for the top face

				st_topPt = nrmPt[0]*self.x_ax[2] + self.crns[4]
				en_topPt = nrmPt[0]*self.x_ax[3] + self.crns[7]

				top_pt = st_topPt + (en_topPt - st_topPt)*nrmPt[1]

				# draw line between bottom point and top point then evaluate location on z axis for the W coordinate
				mapPts.append(bot_pt + (top_pt - bot_pt)*nrmPt[2])
				

			mapPts = np.array(mapPts)

		else:

			mapPts = np.array([])
			self.lock = True
			
		return mapPts



	####
	# getTransform (att, thres)
	# generates transformation matrix covering every vertex in the boxMorph
	# att (numpy array of 3D points) = points representing the attractor curve in question
	# thres (float) = a threshold distance for blending units 
	####

	def getTransform(self, att, thres):

		factors = []

		# generate blend values for every individual vertex

		for i in range(self.src.shape[0]):

			diff = att[:] - self.v[i]
			minimum = np.min(np.linalg.norm(att[:] - self.v[i], axis = 1))

			factor = 1 - minimum/thres
			if factor > 1: factor = 1
			if factor < 0: factor = 0

			factors.append(factor)

		self.transform = np.array(factors)

		return self.transform



	####
	# blendSrc (target, att, thres)
	# deals with blending different units together based off of attractors
	# target (numpy array of 3D points) = points representing the target unit that you are blending towards
	# att (numpy array of 3D points) = points representing the attractor curve in question
	# thres (float) = a threshold distance for blending units 
	####


	def blendSrc(self, target, att, thres):

		self.tar = target

		self.getTransform(att, thres)

		# checks if the number of vertices in each set for blending match
		if self.tar.shape[0] != self.src.shape[0]:
			print("BLENDING TARGETS DON'T MATCH")
			return -1

		# finds the difference between the target and source vertices in the unit box (before morphing)
		vecs = self.tar[:] - self.src[:]

		# scales each vector based on the transform matrix from above
		vecs[:,0] = self.transform[:]*vecs[:,0]
		vecs[:,1] = self.transform[:]*vecs[:,1]
		vecs[:,2] = self.transform[:]*vecs[:,2]

		self.src = self.src[:] + vecs
		self.v = self.mapV(self.src)

		return self.v


	####
	# divide (thres, divide)
	# divides boxes with widths or lengths greater than a specified threshold into 2 - 4 subboxes
	# thres (float) = the threshold above which the box will divide
	# divide (boolean) = if True the box will divide into four regardless of the threshold
	####


	def divide(self, thres = m.pow(10, 1000), divide = False):

		boxes = []
		boxcrns_0, boxcrns_1, boxcrns_2, boxcrns_3 = [], [], [], []

		midB_01 = (self.crns[1] - self.crns[0])/2 + self.crns[0]
		midB_02 = (self.crns[2] - self.crns[1])/2 + self.crns[1] 
		midB_03 = (self.crns[3] - self.crns[2])/2 + self.crns[2]
		midB_04 = (self.crns[0] - self.crns[3])/2 + self.crns[3]

		midT_01 = (self.crns[5] - self.crns[4])/2 + self.crns[4]
		midT_02 = (self.crns[6] - self.crns[5])/2 + self.crns[5]
		midT_03 = (self.crns[7] - self.crns[6])/2 + self.crns[6]
		midT_04 = (self.crns[4] - self.crns[7])/2 + self.crns[7]

		cnt_B = np.multiply((self.crns[0] + self.crns[1] + self.crns[2] + self.crns[3]), .25)
		cnt_T = np.multiply((self.crns[4] + self.crns[5] + self.crns[6] + self.crns[7]), .25)

		# if both width and length greater than threshold or divide on then split the box into four pieces

		if self.L > thres and self.W > thres or divide:

			boxcrns_0 = [self.crns[0], midB_01, cnt_B, midB_04, self.crns[4], midT_01, cnt_T, midT_04]
			boxcrns_1 = [midB_01, self.crns[1], midB_02, cnt_B, midT_01, self.crns[5], midT_02, cnt_T]
			boxcrns_2 = [cnt_B, midB_02, self.crns[2], midB_03, cnt_T, midT_02, self.crns[6], midT_03]
			boxcrns_3 = [midB_04, cnt_B, midB_03, self.crns[3], midT_04, cnt_T, midT_03, self.crns[7]]

			boxes = [boxcrns_0, boxcrns_1, boxcrns_2, boxcrns_3]

		# if length greater than threshold then split box along the Y axis (the length)

		if self.W < thres and self.L > thres and divide == False:

			boxcrns_0 = [self.crns[0], self.crns[1], midB_02, midB_04, self.crns[4], self.crns[5], midT_02, midT_04]
			boxcrns_1 = [midB_04, midB_02, self.crns[2], self.crns[3], midT_04, midT_02, self.crns[6], self.crns[7]]

			boxes = [boxcrns_0, boxcrns_1]

		# if length greater than threshold then split box along the X axis (the width)

		if self.W > thres and self.L < thres and divide == False:

			boxcrns_0 = [self.crns[0], midB_01, midB_03, self.crns[3], self.crns[4], midT_01, midT_03, self.crns[7]]
			boxcrns_1 = [midB_01, self.crns[1], self.crns[2], midB_03, midT_01, self.crns[5], self.crns[6], midT_03]

			boxes = [boxcrns_0, boxcrns_1]

		if self.W < thres and self.L < thres and divide == False:

			boxes = [self.crns]

		return boxes



	####
	# genDefault ()
	# generates a basic box unit of curves to preview the shape of the morphed box
	####

	def genDefault(self):

		members = []

		self.v = self.crns

		for i in range(4):

			members.append([self.v[i], self.v[(i+1)%4]])
			members.append([self.v[i + 4], self.v[(i+1)%4 + 4]])
			members.append([self.v[i], self.v[i+4]])

		return members




# def importMeshOBJ(loc):

# 	f = open(loc, 'r')
# 	allLines = f.read().split('\n')

# 	f = []
# 	v = []

# 	for i in range(len(allLines)):

# 		line = allLines[i].split(' ')

# 		if line[0] == 'v':