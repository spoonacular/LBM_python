import numpy

class Poiseuille2D():
	def __init__(self,SuperGeometry,materialNum,maxVelocity,distance2Wall):
		self.velocityField = numpy.zeros([numpy.maximum(SuperGeometry.getMaterialCoords(materialNum).shape[0],SuperGeometry.getMaterialCoords(materialNum).shape[1]),2])
		self.distance2Wall = distance2Wall
		self.L = self.velocityField.shape[0]-1+self.distance2Wall*2
		self.K = 4*maxVelocity/self.L**2
		for i in numpy.arange(self.velocityField.shape[0]):
			self.x = i + self.distance2Wall
			self.velocityField[i] = self.K*self.x*(self.L-self.x)

	def getVelocityField(self):
		return(self.velocityField)