import numpy
import sys
if __name__ == "__main__":
	sys.path[0] = '/'.join(sys.path[0].split('/')[0:-2]) #go up 2 level
from src.geometry.Geometry2D import *


# class OnLatticeBoundaryCondition(SuperLattice2D):

# 	def __init__(self,SuperLattice2D):
# 		pass


# 	def velocityBoundary(self,SuperGeometry,materialNum,omega):
# 		pass

# 	def pressureBoundary(self,SuperGeometry,materialNum,omega):
# 		pass


class Poiseuille2D():
	def __init__(self,SuperGeometry,materialNum,maxVelocity,distance2Wall):
		self.velocityField = numpy.zeros([numpy.maximum(SuperGeometry.getMaterialCoords(materialNum,SuperGeometry).shape[0],SuperGeometry.getMaterialCoords(materialNum,SuperGeometry).shape[1]),2])
		self.distance2Wall = distance2Wall
		self.L = self.velocityField.shape[0]-1+self.distance2Wall*2
		self.K = 4*maxVelocity/self.L**2
		for i in numpy.arange(self.velocityField.shape[0]):
			self.x = i + self.distance2Wall
			self.velocityField[i] = self.K*self.x*(self.L-self.x)

	def getVelocityField(self):
		return(self.velocityField)

