import sys
import numpy
import os
if __name__ == "__main__":
	sys.path[0] = '/'.join(sys.path[0].split('/')[0:-2]) #go up 2 level
from src.geometry.Geometry2D import *
from src.dynamics.dynamics2D import *
from src.functions.functions2D import *

class SuperLattice2D(SuperGeometry):
	def __init__(self,SuperGeometry):
		self.Cs2 = 1/3
		self.DnQn = 'D2Q9'
		self.distribution = numpy.array([4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36])
		self.cx = 			numpy.array([  0,   1,  0, -1,  0,    1,  -1,  -1,   1])
		self.cy = 			numpy.array([  0,   0,  1,  0, -1,    1,   1,  -1,  -1])
		self.opposite = 	numpy.array([  0,   3,  4,  1,  2,    7,   8,   5,   6]) #should be moved to boundary condition part
		self.dynamics = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1]])
		self.omega = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1]]) #omega = 1/tau, which determines the kinematic viscosity
		self.initialization = 0
		self.rhoMap = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1]])
		self.UMap = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1],2])# velocity field
		self.materialCoordsDic  = SuperGeometry.materialCoordDic #coordinates for each material number
		self.materialMap = SuperGeometry.materialMap


	def defineDynamics(self, SuperGeometry, materialNum, dynamics):
		for coord in self.materialCoordsDic[materialNum]:
			self.dynamics[coord[0],coord[1]] = dynamics.index
			if hasattr(dynamics, 'omega'):
				self.omega[coord[0],coord[1]] = dynamics.omega
			else:
				self.omega[coord[0],coord[1]] = 0

			# indexlist:
			#		index = 0: noDynamics
			# 		index = 1: bulkDynamics
			#		index = 2: bounceBack - wall
			#		index = 3: bounceback - velocity/pressure

	def initialize_(self):
		#define distribution function
		self.f = numpy.zeros([self.rhoMap.shape[0],self.rhoMap.shape[1],9])
		#define surrounding dynamics - for bounceback
		self.surrundingDynamics = numpy.zeros(self.f.shape)
		for i in numpy.arange(9):
			self.f[:,:,i] = numpy.multiply(self.rhoMap,self.distribution[i])
			self.surrundingDynamics[:,:,i] = numpy.roll(numpy.roll(self.dynamics,-self.cx[i],0),-self.cy[i],1)
		self.feq = numpy.array(self.f)

		#define number of bulk lattices
		self.numBulkLattices = self.dynamics[self.dynamics == 1].sum()
		# increase speed with cost of mem
		self.omega_9 = numpy.zeros([self.omega.shape[0],self.omega.shape[1],9]) 
		self.cx_9 = numpy.zeros(self.f.shape)
		self.cy_9 = numpy.zeros(self.f.shape)
		for i in numpy.arange(9):
			self.cx_9[:,:,i] = self.cx[i]
			self.cy_9[:,:,i] = self.cy[i]
			self.omega_9[:,:,i] = self.omega
		
		# self.dynamicsCoordsDic = {}
		# size_ = [1,self.dynamics.shape[0]*self.dynamics.shape[1]]
		# for dynamics in set(numpy.reshape(self.dynamics,size_)[0]):
		# 	self.dynamicsCoordsDic[dynamics] = self.getDynamicsCoords(dynamics)

		self.filterWall = numpy.zeros(self.f.shape)
		self.filterBulk = numpy.ones(self.f.shape)
		for i in self.getDynamicsCoords(2): #bounceback
			self.filterWall[i[0],i[1],:] = 1
			self.filterBulk[i[0],i[1],:] = 0

		self.initialization = 1



	def defineRhoU(self,SuperGeometry,materialNum,rho,u):
		for coord in self.materialCoordsDic[materialNum]:
			self.rhoMap[coord[0],coord[1]] = rho
			self.UMap[coord[0],coord[1],0] = u[0]
			self.UMap[coord[0],coord[1],1] = u[1]

	def defineU_BC(self,SuperGeometry,materialNum,u,BCmethod = 'ZH'):
		self.vloc = 0

		#only support symmetrical velocity profile
		for coord in self.materialCoordsDic[materialNum]:
			if u.size == 2:
				self.UMap[coord[0],coord[1],0] = u[0]
				self.UMap[coord[0],coord[1],1] = u[1]
				self.givenUFindRho(coord[0],coord[1],u,BCmethod)
			else:
				self.UMap[coord[0],coord[1],0] = u[self.vloc,0]
				self.UMap[coord[0],coord[1],1] = u[self.vloc,1]
				self.givenUFindRho(coord[0],coord[1],u[self.vloc,:],BCmethod)
				self.vloc = self.vloc + 1


	def defineRho_BC(self,SuperGeometry,materialNum,rho,BCmethod = 'ZH'):
		for coord in self.materialCoordsDic[materialNum]:
			self.rhoMap[coord[0],coord[1]] = rho
			self.givenRhoFindU(coord[0],coord[1],rho,BCmethod)


	def givenUFindRho(self,i,j,u,BCmethod): #need modidfication
		#find rho
		#arrow = self.findArrow(i,j)

		if BCmethod == 'ZH':
			for p in numpy.arange(1,5):
				if self.dynamics[(i+self.cx[p])%self.dynamics.shape[0],(j+self.cy[p])%self.dynamics.shape[1]] == 1:

					# self.rhoMap[i,j] = 1.5*self.rhoMap[i+self.cx[p],j+self.cy[p]] - 0.5*self.rhoMap[i+2*self.cx[p],j+2*self.cy[p]]
					# # the following codes have sign error !!!!!!
					if p == 1 : #left/west
						self.rhoMap[i,j] = 1/(1-u[0])*(self.f[i,j,0]+self.f[i,j,2]+self.f[i,j,4]+2*(self.f[i,j,3]+self.f[i,j,6]+self.f[i,j,7]))
						self.f[i,j,1] = self.f[i,j,3]+2/3*self.rhoMap[i,j]*u[0]
						self.f[i,j,5] = self.f[i,j,7]-1/2*(self.f[i,j,2]-self.f[i,j,4])+1/6*self.rhoMap[i,j]*u[0]+1/2*self.rhoMap[i,j]*u[1]
						self.f[i,j,8] = self.f[i,j,6]+1/2*(self.f[i,j,2]-self.f[i,j,4])+1/6*self.rhoMap[i,j]*u[0]-1/2*self.rhoMap[i,j]*u[1]
					elif p == 2: #bottom/south
						self.rhoMap[i,j] = 1/(1-u[1])*(self.f[i,j,0]+self.f[i,j,1]+self.f[i,j,3]+2*(self.f[i,j,4]+self.f[i,j,7]+self.f[i,j,8]))
						self.f[i,j,2] = self.f[i,j,4]+2/3*self.rhoMap[i,j]*u[1]
						self.f[i,j,5] = self.f[i,j,7]-1/2*(self.f[i,j,1]-self.f[i,j,3])+1/2*self.rhoMap[i,j]*u[0]+1/6*self.rhoMap[i,j]*u[1]
						self.f[i,j,6] = self.f[i,j,8]+1/2*(self.f[i,j,1]-self.f[i,j,3])-1/2*self.rhoMap[i,j]*u[0]+1/6*self.rhoMap[i,j]*u[1]

					elif p == 3: #right/east
						#print(self.f[i,j,:])
						self.rhoMap[i,j] = 1/(1+u[0])*(self.f[i,j,0]+self.f[i,j,2]+self.f[i,j,4]+2*(self.f[i,j,1]+self.f[i,j,5]+self.f[i,j,8]))
						self.f[i,j,3] = self.f[i,j,1]-2/3*self.rhoMap[i,j]*u[0]
						self.f[i,j,7] = self.f[i,j,5]+1/2*(self.f[i,j,2]-self.f[i,j,4])-1/6*self.rhoMap[i,j]*u[0]-1/2*self.rhoMap[i,j]*u[1]
						self.f[i,j,6] = self.f[i,j,8]-1/2*(self.f[i,j,2]-self.f[i,j,4])-1/6*self.rhoMap[i,j]*u[0]+1/2*self.rhoMap[i,j]*u[1]
						#print(self.rhoMap[i,j])
					elif p == 4: #top/north
						self.rhoMap[i,j] = 1/(1+u[1])*(self.f[i,j,0]+self.f[i,j,1]+self.f[i,j,3]+2*(self.f[i,j,2]+self.f[i,j,5]+self.f[i,j,6]))
						self.f[i,j,4] = self.f[i,j,2]-2/3*self.rhoMap[i,j]*u[1]
						self.f[i,j,7] = self.f[i,j,5]+1/2*(self.f[i,j,1]-self.f[i,j,3])-1/2*self.rhoMap[i,j]*u[0]-1/6*self.rhoMap[i,j]*u[1]
						self.f[i,j,8] = self.f[i,j,6]-1/2*(self.f[i,j,1]-self.f[i,j,3])+1/2*self.rhoMap[i,j]*u[0]-1/6*self.rhoMap[i,j]*u[1]
					# self.rhoMap[i,j] = self.f[i,j].sum()
					# print(self.rhoMap[i,j])


	def givenRhoFindU(self,i,j,rho,BCmethod): #need modidfication
		#find U
		#arrow = findArrow(i,j)
		if BCmethod == 'ZH':
			for p in numpy.arange(1,5):
				if self.dynamics[(i+self.cx[p])%self.dynamics.shape[0],(j+self.cy[p])%self.dynamics.shape[1]] == 1:
					# self.UMap[i,j] = 1.5*self.UMap[i+self.cx[p],j+self.cy[p]]-0.5*self.UMap[i+2*self.cx[p],j+2*self.cy[p]]
					if p == 1:
						self.UMap[i,j,0] = 1 - (self.f[i,j,0]+self.f[i,j,2]+self.f[i,j,4]+2*(self.f[i,j,3]+self.f[i,j,6]+self.f[i,j,7]))/rho
						self.UMap[i,j,1] = 0
						self.f[i,j,1] = self.f[i,j,3]+2/3*rho*self.UMap[i,j,0]
						self.f[i,j,5] = self.f[i,j,7]-1/2*(self.f[i,j,2]-self.f[i,j,4])+1/6*rho*self.UMap[i,j,0]
						self.f[i,j,8] = self.f[i,j,6]+1/2*(self.f[i,j,2]-self.f[i,j,4])+1/6*rho*self.UMap[i,j,0]
					elif p == 2:
						self.UMap[i,j,0] = 0
						self.UMap[i,j,1] = 1 - (self.f[i,j,0]+self.f[i,j,1]+self.f[i,j,3]+2*(self.f[i,j,4]+self.f[i,j,7]+self.f[i,j,8]))/rho
						self.f[i,j,2] = self.f[i,j,4]+2/3*rho*self.UMap[i,j,1]
						self.f[i,j,5] = self.f[i,j,7]-1/2*(self.f[i,j,1]-self.f[i,j,3])+1/6*rho*self.UMap[i,j,1]
						self.f[i,j,6] = self.f[i,j,8]+1/2*(self.f[i,j,1]-self.f[i,j,3])+1/6*rho*self.UMap[i,j,1]
					elif p == 3:
						#print('good')
						self.UMap[i,j,0] = -1 + (self.f[i,j,0]+self.f[i,j,2]+self.f[i,j,4]+2*(self.f[i,j,1]+self.f[i,j,5]+self.f[i,j,8]))/rho
						self.UMap[i,j,1] = 0
						self.f[i,j,3] = self.f[i,j,1]-2/3*rho*self.UMap[i,j,0]
						self.f[i,j,7] = self.f[i,j,5]+1/2*(self.f[i,j,2]-self.f[i,j,4])-1/6*rho*self.UMap[i,j,0]
						self.f[i,j,6] = self.f[i,j,8]-1/2*(self.f[i,j,2]-self.f[i,j,4])-1/6*rho*self.UMap[i,j,0]
					elif p == 4:
						self.UMap[i,j,0] = 0
						self.UMap[i,j,1] = -1 + (self.f[i,j,0]+self.f[i,j,1]+self.f[i,j,3]+2*(self.f[i,j,2]+self.f[i,j,5]+self.f[i,j,6]))/rho
						self.f[i,j,4] = self.f[i,j,2]-2/3*rho*self.UMap[i,j,1]
						self.f[i,j,7] = self.f[i,j,5]+1/2*(self.f[i,j,1]-self.f[i,j,3])-1/6*rho*self.UMap[i,j,1]
						self.f[i,j,8] = self.f[i,j,6]-1/2*(self.f[i,j,1]-self.f[i,j,3])-1/6*rho*self.UMap[i,j,1]
	
	def openBC(self,SuperGeometry,materialNum):
		for coord in self.materialCoordsDic[materialNum]:
			for p in numpy.arange(1,5):
				if self.dynamics[(coord[0]+self.cx[p])%self.dynamics.shape[0],(coord[1]+self.cy[p])%self.dynamics.shape[1]] == 1:
					if p == 1:
						#print(2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],1] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],1])
						self.f[coord[0],coord[1],1] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],1] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],1]
						self.f[coord[0],coord[1],5] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],5] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],5]
						self.f[coord[0],coord[1],8] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],8] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],8]
					elif p == 2:
						self.f[coord[0],coord[1],2] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],2] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],2]
						self.f[coord[0],coord[1],5] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],5] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],5]
						self.f[coord[0],coord[1],6] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],6] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],6]
					elif p == 3:
						self.f[coord[0],coord[1],3] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],3] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],3]
						self.f[coord[0],coord[1],7] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],7] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],7]
						self.f[coord[0],coord[1],6] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],6] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],6]
					elif p == 4:
						self.f[coord[0],coord[1],4] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],4] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],4]
						self.f[coord[0],coord[1],7] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],7] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],7]
						self.f[coord[0],coord[1],8] = 2 * self.f[coord[0]+self.cx[p],coord[1]+self.cy[p],8] - self.f[coord[0]+2*self.cx[p],coord[1]+2*self.cy[p],8]



	def findArrow(self,i,j):
		arrow = numpy.zeros(8)
		for p in numpy.arange(1,9):
			if self.dynamics[(i+self.cx[p])%self.dynamics.shape[0],(j+self.cy[p])%self.dynamics.shape[1]] == 1:
				#print((i+self.cx[p])%self.dynamics.shape[0],(j+self.cy[p])%self.dynamics.shape[1])
				arrow[p-1] = 1

		#print('==========')
		return(arrow)


	def collideAndStream(self):
		pass


	def collide(self):

		#initialize distribution function
		if self.initialization == 0:
			self.initialize_()

		self.updateRhoU() #13
		self.BGKcollide() #25

	def BGKcollide(self):

		self.t1 = self.UMap[:,:,0]*self.UMap[:,:,0] + self.UMap[:,:,1]*self.UMap[:,:,1]
		if not hasattr(self,'t2'):
			self.t2 = numpy.zeros(self.f.shape)

		self.t2 = self.UMap[:,:,0]*numpy.swapaxes([[self.cx]],0,2)+self.UMap[:,:,1]*numpy.swapaxes([[self.cy]],0,2)
		self.feq = self.rhoMap*numpy.swapaxes([[self.distribution]],0,2)*(1+3*self.t2+4.5*self.t2*self.t2-1.5*self.t1)
		self.feq = numpy.swapaxes(numpy.swapaxes(self.feq,0,1),1,2)

		# for k in numpy.arange(9):
		# 	self.t2[:,:,k] = self.UMap[:,:,0]*self.cx[k]+self.UMap[:,:,1]*self.cy[k]
		# 	self.feq[:,:,k] = self.rhoMap*self.distribution[k]*(1+3*self.t2[:,:,k]+4.5*self.t2[:,:,k]*self.t2[:,:,k]-1.5*self.t1)
		self.f = self.f + self.omega_9*(self.feq - self.f)



	def stream(self):
		#bounceback 2 # 75% computation time
		# for i in numpy.arange(self.dynamics.shape[0]):
		# 	for j in numpy.arange(self.dynamics.shape[1]):
		# 		if self.dynamics[i,j] == 2:	
		# 			for k in numpy.arange(1,9):
		# 				#if self.surrundingDynamics[i,j,k] in [1,3,4]: #half way bounceback: modify f on near bulk fluids
		# 				self.f[i,j,k] = self.f[(i+self.cx[k])%self.f.shape[0],(j+self.cy[k])%self.f.shape[1],self.opposite[k]] 

		

		#pre-stream for bounceBack
		self.tmp_f = self.f.copy()
		for k in numpy.arange(1,9):
			self.tmp_f[:,:,self.opposite[k]] = numpy.roll(numpy.roll(self.f[:,:,k],self.cx[k],0),self.cy[k],1)
		self.f = self.tmp_f*self.filterWall + self.f*self.filterBulk

		#stream 
		for k in numpy.arange(1,9):
			self.f[:,:,k] = numpy.roll(numpy.roll(self.f[:,:,k],self.cx[k],0),self.cy[k],1)

		#calculate rhoMap given f after stream
		



	def updateRhoU(self):
		self.rhoMap = self.f.sum(2)
		self.UMap[:,:,0] = (self.f*self.cx_9).sum(2)/self.rhoMap
		self.UMap[:,:,1] = (self.f*self.cy_9).sum(2)/self.rhoMap

	def getAverageRho(self):
		if not hasattr(self,'numBulkLattices'):
			self.initialize_()

		self.averageRho = 0
		for i in numpy.arange(self.dynamics.shape[0]):
			for j in numpy.arange(self.dynamics.shape[1]):
				if self.dynamics[i,j] == 1:
					self.averageRho = self.averageRho + self.rhoMap[i,j]
		return self.averageRho/self.numBulkLattices


	def communicate(self):
		pass

	def executeCoupling(self):
		pass

	def getRhoMap(self):
		if 2 in self.dynamics:
			for coord in self.getDynamicsCoords(2):
				self.rhoMap[coord[0],coord[1]] = 1
		return(self.rhoMap)

	def getUxMap(self):
		if 2 in self.dynamics:
			for coord in self.getDynamicsCoords(2):
				self.UMap[coord[0],coord[1],0] = 0
		return(self.UMap[:,:,0])

	def getUyMap(self):
		if 2 in self.dynamics:
			for coord in self.getDynamicsCoords(2):
				self.UMap[coord[0],coord[1],1] = 0
		return(self.UMap[:,:,1])

	def getSpeedMap(self):
		return(numpy.power(numpy.power(self.UMap[:,:,0],2)+numpy.power(self.UMap[:,:,1],2),1/2))

	 #return N x 2 matrix 
	def getDynamicsCoords(self,dynamicsIndex):
		dynamicsCoords = numpy.zeros([0,2])
		for x in numpy.arange(self.dynamics.shape[0]):
			for y in numpy.arange(self.dynamics.shape[1]):
				if self.dynamics[x,y] == dynamicsIndex:
					dynamicsCoords = numpy.append(dynamicsCoords,[[x,y]],0)
		return numpy.int_(dynamicsCoords)


if __name__ == "__main__":
	#parameters
	numpy.set_printoptions(3)
	nx = 30
	ny = 10
	center_x = 3
	center_y = 3
	radius = 2
	omega = 1
	#define geometry
	topPlate = Indicator.cuboid(0,ny,nx,ny) #x1,y1,x2,y2
	circle = Indicator.circle(center_x,center_y,radius)
	rightPlate = Indicator.cuboid(nx,0,nx,ny)
	leftPlate = Indicator.cuboid(0,0,0,ny)

	cGeometry = CuboidGeometry2D(0,0,nx,ny)
	cGeometry.setPeriodicity()
	superG = SuperGeometry(cGeometry)
	superG.rename(0,5,topPlate)
	superG.rename(0,1,circle)
	superG.rename(0,2)
	superG.rename(2,3,rightPlate)
	superG.rename(2,4,leftPlate)
	#print(superG.materialMap)
	#superG.print()
	print('================================================')
	#lattice
	rho = 1
	u = [0,0.]
	sLattice = SuperLattice2D(superG)
	sLattice.defineRhoU(superG,1,rho,u)
	sLattice.defineRhoU(superG,2,rho,u)
	sLattice.defineRhoU(superG,3,rho,u)
	sLattice.defineRhoU(superG,4,rho,u)
	sLattice.defineRhoU(superG,5,rho,u)
	#print(sLattice.getRhoMap())
	#print(sLattice.getRhoMap().sum())
	#print(sLattice.getUxMap())
	#print(sLattice.getUyMap())
	#print(sLattice.getSpeedMap())
	bulk1 = BGKdynamics(omega)
	

	sLattice.defineDynamics(superG,1,bulk1)# SuperGeometry, materialNum, dynamics
	sLattice.defineDynamics(superG,2,bulk1)
	sLattice.defineDynamics(superG,4,BBvelocity(omega))
	#sLattice.defineDynamics(superG,4,BBpressure(omega))
	sLattice.defineDynamics(superG,3,BBpressure(omega))
	sLattice.defineDynamics(superG,5,BBwall())

	#print(sLattice.dynamics)
	#print(sLattice.omega)
	#print(sLattice.dynamics)

	#print(sLattice.getUyMap())



	outputDirectory = 'data'
	if not os.path.exists(outputDirectory):
		os.makedirs(outputDirectory)

	
	print('initial average rho: {}'.format(sLattice.getAverageRho()))
	maxVelocity = numpy.array([0.1,0])
	#poV = Poiseuille2D(superG,3,maxVelocity,0.5).getVelocityField() #SuperGeometry,materialNum,maxVelocity,distance2Wall
	poV = numpy.array(maxVelocity)

	sLattice.defineU_BC(superG,4,poV)
	sLattice.defineRho_BC(superG,3,1)
	#sLattice.openBC(superG,4)
	#print('initial Ux:\n{}\n==============================='.format(sLattice.getUxMap()))
	numpy.set_printoptions(3)


	print(sLattice.getUxMap())
	print(sLattice.rhoMap)
	for iT in numpy.arange(1000 ):
		# if iT%1000 == 0:
		# 	print('{}/1000'.format(iT))
		# 	print(sLattice.getUxMap())

		sLattice.collide()
		sLattice.stream()
		sLattice.defineU_BC(superG,4,poV)
		sLattice.defineRho_BC(superG,3,1)
		#sLattice.defineRho_BC(superG,3,1)
		#sLattice.openBC(superG,4)
		# if iT%500==0:
		# 	numpy.savetxt('{}/VelocityProfile_{}'.format(outputDirectory,iT),sLattice.getUxMap())
		# 	print('{}/10000'.format(iT))

		#print(sLattice.getUxMap())
		#print(sLattice.getAverageRho())
	#print(poV.getVelocityField())

	#print(sLattice.dynamics)
	#print(sLattice.surrundingDynamics[:,:,1])

	#print(sLattice.f[:,:,1])
	print('===============================final Ux:\n{}\n'.format(sLattice.getUxMap()))
	print('===============================final Uy:\n{}\n'.format(sLattice.getUyMap()))
	#print(sLattice.getUyMap())
	#print(sLattice.getRhoMap())
	#print(sLattice.getRhoMap().sum())
	print('final average rho: {}'.format(sLattice.getAverageRho()))
	print(sLattice.getRhoMap())

	#print(sLattice.dynamics)
	#numpy.savetxt('tmp_file_Ux.txt',sLattice.getUxMap())
	#print(sLattice.dynamics)

