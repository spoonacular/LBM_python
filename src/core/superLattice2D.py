import sys
import numpy
if __name__ == "__main__":
	sys.path[0] = '/'.join(sys.path[0].split('/')[0:-2]) #go up 2 level
from src.geometry.Geometry2D import *



class SuperLattice2D(SuperGeometry):
	def __init__(self,SuperGeometry):
		self.DnQn = 'D2Q9'
		self.opposite = 	numpy.array([  1,   4,  5,  2,  3,    8,   9,   6,   7]) #should be moved to boundary condition part
		self.df = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1],9])

	def defineDynamics(self, SuperGeometry, materialNum, dynamics):
		self.omega = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1]])
		for coord in self.getMaterialCoords(materialNum,SuperGeometry):
			self.omega[coord[0],coord[1]] = dynamics.omega


	def defineRhoU(self,SuperGeometry,materialNum,rho,u):
		self.rhoMap = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1]])
		self.UMap = numpy.zeros([SuperGeometry.materialMap.shape[0],SuperGeometry.materialMap.shape[1],2])

		for coord in self.getMaterialCoords(materialNum,SuperGeometry):
			self.rhoMap[coord[0],coord[1]] = rho
			self.UMap[coord[0],coord[1],0] = u[0]
			self.UMap[coord[0],coord[1],1] = u[1]


	def collideAndStream(self):
		pass

	def collide(self):

		pass

	def stream(self):
		pass


	def communicate(self):
		pass

	def executeCoupling(self):
		pass

	def getRhoMap(self):
		return(self.rhoMap)

	def getUxMap(self):
		return(self.UMap[:,:,0])

	def getUyMap(self):
		return(self.UMap[:,:,1])

	def getSpeedMap(self):
		return(numpy.power(numpy.power(self.UMap[:,:,0],2)+numpy.power(self.UMap[:,:,1],2),1/2))



	@staticmethod #return N x 2 matrix 
	def getMaterialCoords(materialNum,SuperGeometry):
		materialCoords = numpy.zeros([0,2])
		for x in range(SuperGeometry.materialMap.shape[0]):
			for y in range(SuperGeometry.materialMap.shape[1]):
				if SuperGeometry.materialMap[x][y] == materialNum:
					materialCoords = numpy.append(materialCoords,[[x,y]],0)
		return numpy.int_(materialCoords)



if __name__ == "__main__":
	#parameters
	nx = 10
	ny = 5
	center_x = 3
	center_y = 3
	radius = 2
	#define geometry
	topPlate = Indicator.cuboid(0,ny,nx,ny) #x1,y1,x2,y2
	circle = Indicator.circle(center_x,center_y,radius)

	cGeometry = CuboidGeometry2D(0,0,nx,ny)
	cGeometry.setPeriodicity()
	superG = SuperGeometry(cGeometry)
	superG.rename(0,5,topPlate)
	superG.rename(0,1,circle)
	#print(superG.materialMap)
	superG.print()
	print('================================================')
	#lattice
	rho = 1
	u = [2,3]
	sLattice = SuperLattice2D(superG)
	sLattice.defineRhoU(superG,1,rho,u)

	#print(sLattice.getRhoMap())
	#print(sLattice.getUxMap())
	#print(sLattice.getUyMap())
	#print(sLattice.getSpeedMap())

	print(superG.periodicity)

	