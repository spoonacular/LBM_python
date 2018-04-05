from LBM_python import *
import matplotlib.pyplot as plt
#parameters
nx = 200
ny = 100
center_x = 10
center_y = 5
radius = 2
omega1 = 1
omega2 = 1
#define geometry
circle = Indicator.circle(center_x,center_y,radius)
botBox = Indicator.cuboid(80,1,120,40)
bottomPlate = Indicator.cuboid(0,0,nx,0)

cGeometry = CuboidGeometry2D(0,0,nx,ny)
cGeometry.setPeriodicity()
superG = SuperGeometry(cGeometry)
superG.rename(0,1)
superG.rename(1,3,bottomPlate)
superG.rename(1,2,botBox)

print('================================================')
#lattice
rho = 1
rhoZero = 0
u = [0,0.]
sLattice1 = SuperLattice2D(superG)
sLattice2 = SuperLattice2D(superG)

sLattice1.defineRhoU(superG,1,rho,u)
sLattice1.defineRhoU(superG,2,rhoZero,u)
#sLattice1.defineRhoU(superG,3,0.5,u)


sLattice2.defineRhoU(superG,1,rhoZero,u)
sLattice2.defineRhoU(superG,2,rho,u)
#sLattice2.defineRhoU(superG,3,0.5,u)

bulk1 = ShanChenBGKdynamics(omega1)
bulk2 = ShanChenBGKdynamics(omega2)


bounceback1 = BBwall(0.6)
bounceback2 = BBwall(0.4)
sLattice1.defineDynamics(superG,1,bulk1)# SuperGeometry, materialNum, dynamics
sLattice1.defineDynamics(superG,2,bulk1)
sLattice1.defineDynamics(superG,3,bounceback1)

sLattice2.defineDynamics(superG,1,bulk2)
sLattice2.defineDynamics(superG,2,bulk2)
sLattice2.defineDynamics(superG,3,bounceback2)


outputDirectory = 'data'
if not os.path.exists(outputDirectory):
	os.makedirs(outputDirectory)


# maxVelocity = numpy.array([0.01,0])
# poV = Poiseuille2D(superG,3,maxVelocity,0.5).getVelocityField() #SuperGeometry,materialNum,maxVelocity,distance2Wall
# poV = numpy.array(maxVelocity)

G = 3

fig = plt.figure()
# MAIN LOOP


sLattice1.addLatticeCoupling(superG,G,sLattice2)

print('initial average rho1: {0:.5f}'.format(sLattice1.getAverageRho())) #initialize
print('initial average rho2: {0:.5f}'.format(sLattice2.getAverageRho()))

#ok now

for iT in numpy.arange( 2000 ):

	# sLattice.openBC(superG,3)
	# sLattice1.updateRhoU() #13 #needed for defineU_BC / defineRho_BC
	# sLattice2.updateRhoU() #13

	# sLattice.defineU_BC(superG,4,poV)
	# sLattice.defineRho_BC(superG,3,1)

	if iT%50==0:
		im = plt.imshow(sLattice1.rhoMap, animated=True)
		plt.draw()
		if iT == 0:
			plt.pause(5)
		else:
			plt.pause(0.00001)
		# numpy.savetxt('{}/rho1_{}'.format(outputDirectory,iT),sLattice1.getRhoMap())
		print('{}/2000'.format(iT))



	sLattice1.prepareCoupling()
	sLattice2.prepareCoupling()

	sLattice1.executeCoupling(sLattice2)

	sLattice1.collide() 
	sLattice2.collide() 
	sLattice1.stream()
	sLattice2.stream()	
	
# print('===============================final Ux:\n{}\n'.format(sLattice1.getUxMap()))

print('final average rho1: {0:.5f}'.format(sLattice1.getAverageRho()))
print('final average rho2: {0:.5f}'.format(sLattice2.getAverageRho()))

# for i in numpy.arange(sLattice1.rhoMap.shape[0]):
# 	for j in numpy.arange(sLattice1.rhoMap.shape[1]):
# 		print('%10.4f' %sLattice1.rhoMap[i,j],end='')
# 	print('')

print('=========================================================================================')
# for i in numpy.arange(sLattice1.rhoMap.shape[0]):
# 	for j in numpy.arange(sLattice1.rhoMap.shape[1]):
# 		print('%10.4f' %sLattice2.rhoMap[i,j],end='')
# 	print('')