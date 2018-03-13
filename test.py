from LBM_python import *

#input parameter 
nx = 10
ny = 5
center_x = 3
center_y = 3
radius = 2

#geometry definition
topPlate = Indicator.cuboid(0,ny,nx,ny) #x1,y1,x2,y2
circle = Indicator.circle(center_x,center_y,radius)

cGeometry = CuboidGeometry2D(0,0,nx,ny)
cGeometry.setPeriodicity()
superG = SuperGeometry(cGeometry)
superG.rename(0,5,topPlate)
superG.rename(0,1,circle)
#print(superG.materialMap)
superG.print()

#define lattice
sLattice = SuperLattice2D(superG)
print(superG.materialMap)

print('=================')
class test_A():
	tt = 111
	def __init__(self):
		self.x =1
		print('nice')
	def bark(self):
		print('bark')


class parent_A(test_A):
	def __init__(self):
		self.y = 10
		print('end')



pig = parent_A()
pig.bark()
print(pig.tt)