from LBM_python import *
import time                                                

def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print ('%r (%r, %r) %2.2f sec' % (method.__name__, args, kw, te-ts))
              
        return result

    return timed


class Foo(object):

    @timeit
    def foo(self, a=2, b=3):
        pass


numpy.set_printoptions(3)
nx = 10
ny = 5
center_x = 3
center_y = 3
radius = 2
omega = 1
#define geometry
topPlate = Indicator.cuboid(0,ny,nx,ny) #x1,y1,x2,y2
circle = Indicator.circle(center_x,center_y,radius)
rightPlate = Indicator.cuboid(nx,0,nx,ny)

cGeometry = CuboidGeometry2D(0,0,nx,ny)
cGeometry.setPeriodicity()
superG = SuperGeometry(cGeometry)
superG.rename(0,5,topPlate)
superG.rename(0,1,circle)
superG.rename(0,2)
superG.rename(2,3,rightPlate)
#print(superG.materialMap)
superG.print()
print('================================================')
#lattice
rho = 1
u = [0,0.]
sLattice = SuperLattice2D(superG)
sLattice.defineRhoU(superG,1,1.0,u)
#print(sLattice.getRhoMap())
sLattice.defineRhoU(superG,5,1.0,u)
sLattice.defineRhoU(superG,3,1,u)
sLattice.defineRhoU(superG,2,rho,u)
#print(sLattice.getRhoMap())
#print(sLattice.getRhoMap().sum())
#print(sLattice.getUxMap())
#print(sLattice.getUyMap())
#print(sLattice.getSpeedMap())
bulk1 = BGKdynamics(omega)
bounceb = bounceBack()
sLattice.defineDynamics(superG,1,bulk1)# SuperGeometry, materialNum, dynamics
sLattice.defineDynamics(superG,2,bulk1)
sLattice.defineDynamics(superG,5,bounceb)
#print(sLattice.dynamics)
#print(sLattice.omega)
#print(sLattice.dynamics)
#print(sLattice.getUxMap())
#print(sLattice.getUyMap())
#print(sLattice.getAverageRho())
#print(sLattice.getRhoMap())
maxVelocity = numpy.array([-0.1,0])
poV = Poiseuille2D(superG,3,maxVelocity,0.5) #SuperGeometry,materialNum,maxVelocity,distance2Wall

@timeit
def f1():
	pass

	#print(sLattice.getAverageRho())
@timeit
def collide(iter):
	for i in numpy.arange(iter):
		sLattice.collide()
		
		
################################################################################################
@timeit
def stream(iter):
	sLattice.stream()
@timeit
def define(iter):
	sLattice.defineU(superG,3,poV.getVelocityField())

f1()
collide(10)
stream(10000)
define(10000)
Foo().foo()





