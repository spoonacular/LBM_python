import matplotlib.pyplot as plt
from LBM_python import *



#not really a rayleigh taylor :(



#parameters
nx = 400
ny = 200

omega1 = 1
omega2 = 1
#define geometry

bothalf = Indicator.cuboid(int(nx/3),0,int(nx/2),int(ny/2))
topPlate = Indicator.cuboid(0,ny,nx,ny)
bottomPlate = Indicator.cuboid(0,0,nx,0)

cGeometry = CuboidGeometry2D(0,0,nx,ny)
cGeometry.setPeriodicity()
superG = SuperGeometry(cGeometry)
superG.rename(0,1)
superG.rename(1,3,bottomPlate)
superG.rename(1,4,topPlate)
superG.rename(1,2,bothalf)
#444444444444444444444444
#111111111111111111111111
#111111111111111111111111
#222222222222222222222222
#222222222222222222222222
#333333333333333333333333
print('================================================')
#lattice
LatticeTop = SuperLattice2D(superG)
LatticeBottom = SuperLattice2D(superG)

rho = 1
rhoZero = 0
u = [0,0.]
f = [0,-30/(ny*ny)]

LatticeTop.defineRhoU(superG,1,rho,u)
LatticeTop.defineRhoU(superG,2,rhoZero,u)
LatticeTop.defineExternalField(superG,1,f)
LatticeTop.defineExternalField(superG,2,f)	

LatticeBottom.defineRhoU(superG,1,rhoZero,u)
LatticeBottom.defineRhoU(superG,2,rho,u)


bulk1 = ShanChenBGKdynamics(omega1)
bulk2 = ShanChenBGKdynamics(omega2)


bounceback1 = BBwall(1)
bounceback2 = BBwall(0)
LatticeTop.defineDynamics(superG,1,bulk1)# SuperGeometry, materialNum, dynamics
LatticeTop.defineDynamics(superG,2,bulk1)
LatticeTop.defineDynamics(superG,3,bounceback1)
LatticeTop.defineDynamics(superG,4,bounceback2)

LatticeBottom.defineDynamics(superG,1,bulk2)
LatticeBottom.defineDynamics(superG,2,bulk2)
LatticeBottom.defineDynamics(superG,3,bounceback2)
LatticeBottom.defineDynamics(superG,4,bounceback1)

outputDirectory = 'data'
if not os.path.exists(outputDirectory):
	os.makedirs(outputDirectory)
outputDirectory = 'figures'
if not os.path.exists(outputDirectory):
	os.makedirs(outputDirectory)


G = 3

fig = plt.figure()
# MAIN LOOP


LatticeTop.addLatticeCoupling(superG,G,LatticeBottom)

print('initial average rho1: {0:.5f}'.format(LatticeTop.getAverageRho())) #initialize
print('initial average rho2: {0:.5f}'.format(LatticeBottom.getAverageRho()))

for iT in numpy.arange( 20000 ):

	# sLattice.openBC(superG,3)
	# LatticeTop.updateRhoU() #13 #needed for defineU_BC / defineRho_BC
	# LatticeBottom.updateRhoU() #13

	# sLattice.defineU_BC(superG,4,poV)
	# sLattice.defineRho_BC(superG,3,1)

	# if iT == 2:
	# 	print(LatticeBottom.externalField[:,:,1])
	# 	# print(LatticeTop.materialCoordsDic)

	if iT%100==0:
		im = plt.imshow(LatticeTop.rhoMap, animated=True)
		#plt.draw()

		# if iT == 0:
		# 	plt.pause(1)
		# else:
		# 	plt.pause(0.00001)
		# numpy.savetxt('{}/rho1_{}'.format(outputDirectory,iT),LatticeTop.getRhoMap())
		fig.savefig('figures/density_{0:.0f}.png'.format(iT))
		print('{}/20000'.format(iT))


	LatticeTop.prepareCoupling()
	LatticeBottom.prepareCoupling()

	LatticeTop.executeCoupling(LatticeBottom)

	LatticeTop.collide() 
	LatticeBottom.collide() 
	LatticeTop.stream()
	LatticeBottom.stream()	
	
# print('===============================final Ux:\n{}\n'.format(LatticeTop.getUxMap()))

print('final average rho1: {0:.5f}'.format(LatticeTop.getAverageRho()))
print('final average rho2: {0:.5f}'.format(LatticeBottom.getAverageRho()))



# for i in numpy.arange(LatticeTop.rhoMap.shape[0]):
# 	for j in numpy.arange(LatticeTop.rhoMap.shape[1]):
# 		print('%10.4f' %LatticeTop.rhoMap[i,j],end='')
# 	print('')

print('=========================================================================================')
# for i in numpy.arange(LatticeTop.rhoMap.shape[0]):
# 	for j in numpy.arange(LatticeTop.rhoMap.shape[1]):
# 		print('%10.4f' %LatticeBottom.UMap[i,j,0],end='')
# 	print('')
