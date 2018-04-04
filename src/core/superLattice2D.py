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
		#self.omega = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1]]) #omega = 1/tau, which determines the kinematic viscosity
		self.initialization = 0
		self.rhoMap = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1]])
		self.UMap = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1],2])# velocity field
		self.materialCoordsDic  = SuperGeometry.materialCoordsDic #coordinates for each material number
		self.materialMap = SuperGeometry.materialMap

	def defineDynamics(self, SuperGeometry, materialNum, dynamics):
		for coord in self.materialCoordsDic[materialNum]:
			self.dynamics[coord[0],coord[1]] = dynamics.index
			if hasattr(dynamics, 'omega'):
				self.omega = dynamics.omega
			else:
				self.omega = 0

			# indexlist:
			#		index = 0: noDynamics
			# 		index = 1: bulkDynamics
			#		index = 2: bounceBack - wall
			#		index = 3: bounceback - velocity/pressure

	def initialize_(self):
		#define distribution function
		self.f = numpy.zeros([self.rhoMap.shape[0],self.rhoMap.shape[1],9])
		self.tmp_f = self.f.copy()
		#define surrounding dynamics - for bounceback
		self.surrundingDynamics = numpy.zeros(self.f.shape)
		for i in numpy.arange(9):
			self.f[:,:,i] = numpy.multiply(self.rhoMap,self.distribution[i])
			self.surrundingDynamics[:,:,i] = numpy.roll(numpy.roll(self.dynamics,-self.cx[i],0),-self.cy[i],1)
		self.feq = numpy.array(self.f)

		#define number of bulk lattices
		self.numBulkLattices = self.dynamics[self.dynamics == 1].sum()
		# increase speed with cost of mem
		self.cx_9 = numpy.zeros(self.f.shape)
		self.cy_9 = numpy.zeros(self.f.shape)
		for i in numpy.arange(9):
			self.cx_9[:,:,i] = self.cx[i]
			self.cy_9[:,:,i] = self.cy[i]

		
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
					if p == 1 : #left/west
						#print('good')
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
						self.f[i,j,7] = self.f[i,j,5]+1/2*(self.f[i,j,2]-self.f[i,j,4])-1/6*rho*self.UMap[i,j,0] # - 1/2*rho*self.UMap[i,j,1]
						self.f[i,j,6] = self.f[i,j,8]-1/2*(self.f[i,j,2]-self.f[i,j,4])-1/6*rho*self.UMap[i,j,0] # + 1/2*rho*self.UMap[i,j,1]
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


		self.BGKcollide() #25

	def BGKcollide(self):

		self.t1 = self.UMap[:,:,0]*self.UMap[:,:,0] + self.UMap[:,:,1]*self.UMap[:,:,1]
		self.t2 = self.UMap[:,:,0]*numpy.swapaxes([[self.cx]],0,2)+self.UMap[:,:,1]*numpy.swapaxes([[self.cy]],0,2) # (9,x,y)
		self.feq = self.rhoMap*numpy.swapaxes([[self.distribution]],0,2)*(1+3*self.t2+4.5*self.t2*self.t2-1.5*self.t1) # (9,x,y)
		self.feq = numpy.swapaxes(numpy.swapaxes(self.feq,0,1),1,2) # (x,y,9)

		# for k in numpy.arange(9):
		# 	self.t2[:,:,k] = self.UMap[:,:,0]*self.cx[k]+self.UMap[:,:,1]*self.cy[k]
		# 	self.feq[:,:,k] = self.rhoMap*self.distribution[k]*(1+3*self.t2[:,:,k]+4.5*self.t2[:,:,k]*self.t2[:,:,k]-1.5*self.t1)
		self.f = self.f + self.omega*(self.feq - self.f)

	def stream(self):
		#bounceback 2 # 75% computation time
		# for i in numpy.arange(self.dynamics.shape[0]):
		# 	for j in numpy.arange(self.dynamics.shape[1]):
		# 		if self.dynamics[i,j] == 2:	
		# 			for k in numpy.arange(1,9):
		# 				#if self.surrundingDynamics[i,j,k] in [1,3,4]: #half way bounceback: modify f on near bulk fluids
		# 				self.f[i,j,k] = self.f[(i+self.cx[k])%self.f.shape[0],(j+self.cy[k])%self.f.shape[1],self.opposite[k]] 

		
		#pre-stream for bounceBack
		
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

	def addLatticeCoupling(self,SuperGeometry,G,SuperLattice2D):
		self.G = G
		self.G_omega = G/self.omega
		SuperLattice2D.G_omega = G/SuperLattice2D.omega

	def prepareCoupling(self):
		#calculate moment
		self.rhoMap = self.f.sum(2)
		self.momentX = (self.cx * self.f).sum(2) # (x,y)
		self.momentY = (self.cy * self.f).sum(2) # (x,y)

		self.tmp_rho = self.rhoMap * numpy.swapaxes([[self.distribution]],0,2) # (9,x,y)

		self.rhoContribX =  numpy.zeros(self.rhoMap.shape)
		self.rhoContribY =  numpy.zeros(self.rhoMap.shape)
		for k in numpy.arange(1,9):  # SC force component excluding " -Psi*G "


			self.rhoContribX = self.rhoContribX - numpy.roll(numpy.roll(self.tmp_rho[k,:,:],self.cx[k],0),self.cy[k],1) * self.cx[k]
			self.rhoContribY = self.rhoContribY - numpy.roll(numpy.roll(self.tmp_rho[k,:,:],self.cx[k],0),self.cy[k],1) * self.cy[k]



			# # or: 
			# self.rhoContribX = self.rhoContribX + numpy.roll(numpy.roll(self.tmp_rho[k,:,:],self.cx[self.opposite[k]],0),self.cy[self.opposite[k]],1) * self.cx[k]
			# self.rhoContribY = self.rhoContribY + numpy.roll(numpy.roll(self.tmp_rho[k,:,:],self.cx[self.opposite[k]],0),self.cy[self.opposite[k]],1) * self.cy[k]




	def executeCoupling(self,SuperLattice2D):

		self.totalRho_omega = self.rhoMap*self.omega + SuperLattice2D.rhoMap*SuperLattice2D.omega

		self.commonUx = (self.momentX*self.omega + SuperLattice2D.momentX*SuperLattice2D.omega)/ self.totalRho_omega
		self.commonUy = (self.momentY*self.omega + SuperLattice2D.momentY*SuperLattice2D.omega)/ self.totalRho_omega

		self.UMap[:,:,0] = self.commonUx - self.G_omega*SuperLattice2D.rhoContribX
		self.UMap[:,:,1] = self.commonUy - self.G_omega*SuperLattice2D.rhoContribY
		SuperLattice2D.UMap[:,:,0] = self.commonUx - SuperLattice2D.G_omega*self.rhoContribX
		SuperLattice2D.UMap[:,:,1] = self.commonUy - SuperLattice2D.G_omega*self.rhoContribY


if __name__ == "__main__":

	#parameters
	numpy.set_printoptions(3)
	nx = 200
	ny = 200
	center_x = 5
	center_y = 5
	radius = 2

	omega1 = 0.55
	omega2 = 0.55
	#define geometry
	circle = Indicator.circle(center_x,center_y,radius)
	centerBox = Indicator.cuboid(80,80,120,120)

	cGeometry = CuboidGeometry2D(0,0,nx,ny)
	cGeometry.setPeriodicity()
	superG = SuperGeometry(cGeometry)
	superG.rename(0,1)
	superG.rename(1,2,centerBox)

	print('================================================')
	#lattice
	rho = 1
	rhoZero = 0
	u = [0,0.]
	sLattice1 = SuperLattice2D(superG)
	sLattice2 = SuperLattice2D(superG)

	sLattice1.defineRhoU(superG,1,rho,u)
	sLattice1.defineRhoU(superG,2,rhoZero,u)

	sLattice2.defineRhoU(superG,1,rhoZero,u)
	sLattice2.defineRhoU(superG,2,rho,u)


	bulk1 = ShanChenBGKdynamics(omega1)
	bulk2 = ShanChenBGKdynamics(omega2)

	sLattice1.defineDynamics(superG,1,bulk1)# SuperGeometry, materialNum, dynamics
	sLattice1.defineDynamics(superG,2,bulk1)
	sLattice2.defineDynamics(superG,1,bulk2)
	sLattice2.defineDynamics(superG,2,bulk2)


	outputDirectory = 'data'
	if not os.path.exists(outputDirectory):
		os.makedirs(outputDirectory)

	
	print('initial average rho1: {0:.5f}'.format(sLattice1.getAverageRho())) #initialize
	print('initial average rho2: {0:.5f}'.format(sLattice2.getAverageRho()))
	# maxVelocity = numpy.array([0.01,0])
	#poV = Poiseuille2D(superG,3,maxVelocity,0.5).getVelocityField() #SuperGeometry,materialNum,maxVelocity,distance2Wall
	# poV = numpy.array(maxVelocity)

	numpy.set_printoptions(3)



	G = 3
	sLattice1.addLatticeCoupling(superG,G,sLattice2)

	# MAIN LOOP
	for iT in numpy.arange( 5000 ):

		# sLattice.openBC(superG,3)
		# sLattice1.updateRhoU() #13 #needed for defineU_BC / defineRho_BC
		# sLattice2.updateRhoU() #13


		# sLattice.defineU_BC(superG,4,poV)
		# sLattice.defineRho_BC(superG,3,1)

		if iT%100==0:
			numpy.savetxt('{}/rho1_{}'.format(outputDirectory,iT),sLattice1.getRhoMap())
			print('{}/5000'.format(iT))


		sLattice1.prepareCoupling()
		sLattice2.prepareCoupling()
		sLattice1.executeCoupling(sLattice2)

		sLattice1.collide() #update rhoU, collide
		sLattice2.collide() #update rhoU, collide
		sLattice1.stream()
		sLattice2.stream()	
		
	# print('===============================final Ux:\n{}\n'.format(sLattice1.getUxMap()))

	print('final average rho1: {0:.5f}'.format(sLattice1.getAverageRho()))
	print('final average rho2: {0:.5f}'.format(sLattice2.getAverageRho()))

	# for i in numpy.arange(sLattice1.rhoMap.shape[0]):
	# 	for j in numpy.arange(sLattice1.rhoMap.shape[1]):
	# 		print('%10.4f' %sLattice1.rhoMap[i,j],end='')
	# 	print('')

	# print('=========================================================================================')
	# for i in numpy.arange(sLattice1.rhoMap.shape[0]):
	# 	for j in numpy.arange(sLattice1.rhoMap.shape[1]):
	# 		print('%10.4f' %sLattice2.rhoMap[i,j],end='')
	# 	print('')

