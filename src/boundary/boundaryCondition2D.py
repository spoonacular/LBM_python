import numpy`
import sys
if __name__ == "__main__":
	sys.path[0] = '/'.join(sys.path[0].split('/')[0:-2]) #go up 2 level
from src.geometry.Geometry2D import *
from src.core.superLattice2D import *



# class addBoundaryCondtion():

# 	def __init__(self):
# 		pass

# 	def addVelocityBoundary(self,SuperGeometry,materialNum,omega):
# 		return(self.toBeFind = 'rho')

# 	def addPressureBoundary(self,SuperGeometry,materialNum,omega):
# 		return(self.toBeFind = 'u')


class OnLatticeBoundaryCondtion(SuperLattice2D):
	def __init__(self,SuperLattice2D):
		self.paras = numpy.array([[0,2,4,3,6,7],[0,1,3,4,7,8],[0,2,4,1,5,8],[0,1,3,2,5,6]]) #find Rho/U
		self.tmp_f = SuperLattice2D.f
		self.tmp_dynamics = SuperLattice2D.dynamics

	def addVelocityBoundary(self):
		#find rho
		for i in numpy.arange(SuperLattice2D.dynamics.shape[0]):
			for j in numpy.arange(SuperLattice2D.dynamics.shape[1]):
				if SuperLattice2D.dynamics[i,j] == 2:	
					self.arrow = findArrow(i,j)
					SuperLattice2D.rhoMap[i,j] = 1/(1-SuperLattice2D.UMap[i,j,0])*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][0]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][1]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][2]]+ 2*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][3]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][4]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][5]]
							))+																																1/(1-SuperLattice2D.UMap[i,j,1])*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][0]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][1]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][2]]+ 2*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][3]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][4]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][5]]
							))

	def addPressureBoundry(self):
		for i in numpy.arange(SuperLattice2D.dynamics.shape[0]):
			for j in numpy.arange(SuperLattice2D.dynamics.shape[1]):
				if SuperLattice2D.dynamics[i,j] == 2:	
					self.arrow = findArrow(i,j)			
					SuperLattice2D.UMap[i,j] = [
						-1+1/SuperLattice2D.rhoMap[i,j]*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][0]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][1]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][2]]+ 2*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][3]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][4]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[0]][5]]
							)),
						-1+1/SuperLattice2D.rhoMap[i,j]*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][0]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][1]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][2]]+ 2*(
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][3]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][4]]+
							SuperLattice2D.f[i,j,self.paras[self.arrow[1]][5]]
							))
					]


	def findArrow(self,i,j):
		self.arrow = numpy.zeros(2)
		if self.dynamics[i+1,j] == 2:
			self.arrow[0] = 0
		elif self.dynamics[i-1,j] == 2:
			self.arrow[0] = 2
		if self.dynamics[i,j+1] == 2:
			self.arrow[1] = 1
		elif self.dynamics[i,j-1] == 2:
			self.arrwo[1] = 3

		return(self.arrow)



class BounceBackMethod(): #bounce back approach - wall at midway
	def __init__(self,u):
		pass

	def findRho(self):
		pass

	def findU(self):
		pass


class ZouHeMethod(): #wet node approach - wall at boundary
	def __init__(self):
		pass

	def findRho(self):
		pass


	def findU(self):
		pass



