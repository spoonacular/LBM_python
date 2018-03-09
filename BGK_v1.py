import numpy
import matplotlib 
# test : parabolic profile in Poiseuille flow using lBM

# Parameters for the simulation setup, all in lattice units
nx = 200
ny = 50
maxIter = 10000
writeDataIter = 100
showDataIter = 50
omega = 1
rhoWater = 1 	#fluid 1
rhoOil = 0.8 	#fluid 2
inject_V = 0.05
rhoWaterBB == 0.5

def prepareGeometry(superGeometry):
	print('prepareGeometry')
	topPlate = Indicator.cuboid([0,ny],[nx,ny])
	botPlate = Indicator.cuboid([0,0],[nx,0])
	leftInlet = Indicator.cuboid([0,1],[0,199])
	rightOutlet = Indicator.cuboid([200,1],[200,199])
	middleCircle = Indicator.circle([25,25],10)

	bouncebacks = topPlate + botPlate + middleCircle

	superGeometry.rename(0,1)
	superGeometry.rename(1,3,leftInlet)
	superGeometry.rename(1,4,rightOutlet)
	superGeometry.rename(1,5,bouncebacks)

	superGeometry.clean() #assume material 0 is getNoDynamics
	superGeometry.print()

	print('Prepare Geometry ... OK')

def prepareLattice( sLattice1, sLattice2, bulkDynamics1,bulkDynamics2, onLatticeBoundaryCondition1, onLatticeBoundaryCondition2, bounceBackRho1, bounceBackRho2, superGeometry):
	print('prepareLattice')

		#no dynamics
	sLattice1.defineDynamics( superGeometry,0,getNoDynamics)
	sLattice2.defineDynamics( superGeometry,0,getNoDynamics)

		#bulkDynamics
	sLattice1.defineDynamics( superGeometry,1,bulkDynamics1)
	#sLattice1.defineDynamics( superGeometry,2,bulkDynamics1)
	sLattice1.defineDynamics( superGeometry,3,bulkDynamics1)
	sLattice1.defineDynamics( superGeometry,4,bulkDynamics1)
	sLattice2.defineDynamics( superGeometry,1,bulkDynamics1)
	#sLattice2.defineDynamics( superGeometry,2,bulkDynamics1)
	sLattice2.defineDynamics( superGeometry,3,bulkDynamics1)
	sLattice2.defineDynamics( superGeometry,4,bulkDynamics1)

		#boundaries
	sLattice1.defineDynamics( superGeometry,5,bounceBackRho1)
	sLattice2.defineDynamics( superGeometry,5,bounceBackRho2)

	bc_1.addVelocityBoundary( superGeometry, 3, omega)
	bc_1.addPressureBoundary( superGeometry, 4, omega)
	bc_2.addPressureBoundary( superGeometry, 4, omega)

	print('Prepare Lattice ... OK')


def setBoundaryValues( sLattice1, sLattice2, iT, superGeometry):
	if iT == 0:
		rho1 = rhoWater
		rho2 = rhoOil
		vInitial = inject_V
		zeroV = 0





