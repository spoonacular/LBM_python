import numpy

class Indicator():
	def __init__(self,_shape=0,p1=0,p2=0,p3=0,p4=0):
		self._shape = _shape
		if _shape:
			if _shape == 'cuboid':
				self.vertices = [p1,p2,p3,p4] #x1 y1 x2 y2
			elif _shape == 'circle':
				self.center = [p1,p2]
				self.radius = p3

	@classmethod
	def cuboid(cls,x1,y1,x2,y2):
		return cls('cuboid',x1,y1,x2,y2)

	@classmethod
	def circle(cls,x,y,radius):
		return cls('circle',x,y,radius)

class CuboidGeometry2D():
	dimension = 2
	periodicity = [False, False]
	def __init__(self,x1,x2,y1,y2):
		self.boundary = [x1, x2, y1+1, y2+1]

	def setPeriodicity(self,tf_x = False, tf_y = False):
		periodicity = [tf_x, tf_y]

class SuperGeometry(CuboidGeometry2D):
	def __init__(self,CuboidGeometry2D):
		self.boundary = CuboidGeometry2D.boundary
		self.materialMap = numpy.zeros([self.boundary[2]-self.boundary[0],self.boundary[3]-self.boundary[1]])

	def rename(self,fromM,toM,*args):
		if len(args) == 0:
			numpy.place(self.materialMap,self.materialMap==fromM,toM)
		elif len(args) == 1:
			if args[0]._shape == 'cuboid':
				toBeReplaced = self.materialMap[args[0].vertices[0]:args[0].vertices[2]+1,args[0].vertices[1]:args[0].vertices[3]+1]
				numpy.place(toBeReplaced,toBeReplaced==fromM,toM)
			elif args[0]._shape == 'circle':
				center_x0 = args[0].center[0]
				center_y0 = args[0].center[1]
				for coord_x in numpy.arange(self.materialMap.shape[0]):
					for coord_y in numpy.arange(self.materialMap.shape[1]):
						distance = numpy.sqrt(numpy.power((numpy.float64(coord_x)-center_x0),2)+numpy.power((numpy.float64(coord_y)-center_y0),2))
						if self.materialMap[coord_x][coord_y] == fromM:
							if distance <= args[0].radius:
								self.materialMap[coord_x][coord_y] = toM
		else:
			print('rename: unidentified number of input')

	@staticmethod #return N x 2 matrix 
	def getMaterialCoords(materialNum,SuperGeometry):
		materialCoords = numpy.zeros([0,2])
		for x in numpy.arange(SuperGeometry.materialMap.shape[0]):
			for y in numpy.arange(SuperGeometry.materialMap.shape[1]):
				if SuperGeometry.materialMap[x][y] == materialNum:
					materialCoords = numpy.append(materialCoords,[[x,y]],0)
		return numpy.int_(materialCoords)
		
	def print(self):
		non_duplicate = list(map(numpy.unique, self.materialMap.reshape(1,self.materialMap.size)))

		for i in non_duplicate[0]:
			count_material = numpy.isin(self.materialMap, i)
			print('Number of material {}: {}'.format(int(i),sum(sum(count_material))))


if __name__ == "__main__":

	nx = 10
	ny = 5
	center_x = 3
	center_y = 3
	radius = 2

	topPlate = Indicator.cuboid(0,ny,nx,ny) #x1,y1,x2,y2
	circle = Indicator.circle(center_x,center_y,radius)

	cGeometry = CuboidGeometry2D(0,0,nx,ny)
	cGeometry.setPeriodicity()
	superG = SuperGeometry(cGeometry)
	superG.rename(0,5,topPlate)
	superG.rename(0,1,circle)
	print(superG.materialMap)
	superG.print()




