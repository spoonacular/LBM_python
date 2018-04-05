import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from LBM_python import *



#parameters
numpy.set_printoptions(3)
nx = 200
ny = 200
center_x = 5
center_y = 5
radius = 2
omega1 = 1
omega2 = 1
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


G = 3
sLattice1.addLatticeCoupling(superG,G,sLattice2)

# MAIN LOOP
fig = plt.figure()
ims = []
for iT in numpy.arange( 1000 ):


	if iT%50==0:
		im = plt.imshow(sLattice1.rhoMap, animated=True)
		ims.append([im])
		plt.draw()
		plt.pause(0.00001)
		# plt.draw()
		# plt.pause(0.01)
		# numpy.savetxt('{}/rho1_{}'.format(outputDirectory,iT),sLattice1.getRhoMap())
		print('{}/1000'.format(iT))

	sLattice1.prepareCoupling()
	sLattice2.prepareCoupling()
	sLattice1.executeCoupling(sLattice2)

	sLattice1.collide() 
	sLattice2.collide() 
	sLattice1.stream()
	sLattice2.stream()	

ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=1000)
plt.show()

# print('===============================final Ux:\n{}\n'.format(sLattice1.getUxMap()))

print('final average rho1: {0:.5f}'.format(sLattice1.getAverageRho()))
print('final average rho2: {0:.5f}'.format(sLattice2.getAverageRho()))

plt.show()