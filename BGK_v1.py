import numpy
import matplotlib 
# test : parabolic profile in Poiseuille flow using lBM
# Parameters for the simulation setup, all in lattice units
nx = 200
ny = 50
maxIter = 10000

dataIter = 100
statIter = 50

omega = 1
rhoWater = 1 	#fluid 1
rhoOil = 0.8 	#fluid 2
inject_V = 0.05
rho1_BB == 0.5

def prepareGeometry(superGeometry):
	print('prepareGeometry')
	topPlate = Indicator.cuboid([0,ny],[nx,ny])
	botPlate = Indicator.cuboid([0,0],[nx,0])
	leftInlet = Indicator.cuboid([0,1],[0,199])
	rightOutlet = Indicator.cuboid([200,1],[200,199])
	middleCircle = Indicator.circle([25,25],10)

	superGeometry.rename(0,1)
	superGeometry.rename(1,3,leftInlet)
	superGeometry.rename(1,4,rightOutlet)
	superGeometry.rename(1,5,topPlate)
	superGeometry.rename(1,5,botPlate)
	superGeometry.rename(1,5,middleCircle)

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

	onLatticeBoundaryCondition1.addVelocityBoundary( superGeometry, 3, omega)
	onLatticeBoundaryCondition1.addPressureBoundary( superGeometry, 4, omega)
	onLatticeBoundaryCondition2.addPressureBoundary( superGeometry, 4, omega)

	print('Prepare Lattice ... OK')


def setBoundaryValues( sLattice1, sLattice2, iT, superGeometry):
	if iT == 0:
		rho1 = rhoWater
		rho2 = rhoOil
		vInitial = [inject_V,0]
		zeroV = [0,0]

		sLattice1.defineRhoU( superGeometry, 1, rho1, zeroV)
		#sLattice1.defineRhoU( superGeometry, 2, zero, zeroV)
		sLattice1.defineRhoU( superGeometry, 3, rho1, vInitial)
		sLattice1.defineRhoU( superGeometry, 4, rho1, zeroV)

		#sLattice2.defineRhoU( superGeometry, 1, zero, zeroV)
		#sLattice2.defineRhoU( superGeometry, 2, zero, zeroV)
		#sLattice2.defineRhoU( superGeometry, 3, zero, zeroV)
		#sLattice2.defineRhoU( superGeometry, 4, zero, zeroV)


def getResults(sLattice1,sLattice2,iT,superGeometry,timer):
	print('getResults')
	vtmWriter = VTMwriter2D('fluid_1')
	gifWriter = GIFwriter2D('fluid_1')
	if (iT == 0):
		vtmWriter.geometry(sLattice1,superGeometry)
		vtmWriter.write( geometry )
	if (iT%statIter == 0)
		timer.update( iT )
		timer.printStep()

		print('average rho of fluid 1 = {}'.format(sLattice1.getStatistics().getAverageRho()))
		print('average rho of fluid 2 = {}'.format(sLattice2.getStatistics().getAverageRho()))

	if (it%dataIter == 0)
		velocity = SuperLatticeVelocity(sLattice1)
		density = SuperLatticeDensity(sLattice2)
		vtmWriter.write(velocity,iT)
		vtmWriter.write(density,iT)
		gifWriter.write(density,iT)

	printf('writing data ... OK')


def main():
	outputDir = './CPv1/'
	omega1 = 1
	omega2 = 1
	G = 3

	cGeometry = CuboidGeometry2D(0,0,nx,ny)
	cGeometry.setPeriodicity( False, False) #perodic in x and y direction

	superGeometry = superGeometry(cGeometry)
	prepareGeometry( superGeometry )

	sLattice1 = SuperLattice2D( superGeometry )
	sLattice2 = SuperLattice2D( superGeometry )

	bulk1 = BGKdynamics(omega1)
	bulk2 = BGKdynamics(omega2)

	bc1 = onLatticeBoundaryCondition2D(sLattice1)
	bc2 = onLatticeBoundaryCondition2D(sLattice2)

	bb1 = BounceBack( rho1_BB )
	bb2 = BounceBack( 1-rho1_BB )
	

	rho0 = [rho1, rho2]
	interactionPotential = PsiEqualsRho()
	coupling = ShanChen93(G,rho0, interactionPotential)

	prepareLattice(sLattice1,sLattice2,bulk1,bulk2,bc1,bc2,bb1,bb2,superGeometry)


	iT = 0
	timer = Timer(maxIter,superGeometry.getStatistics().getNvoxel())
	timer.start()

	while iT < maxIter:
		setBoundaryValues(sLattice1,sLattice2,iT,superGeometry)

		sLattice1.collideAndStream()
		sLattice2.collideAndStream()

		sLattice1.communicate()
		sLattice2.communicate()

		sLattice1.executeCoupling()

		getResults(sLattice1,sLattice2,iT,superGeometry,timer)

	timer.stop()
	timer.printSummary()



if __name__ == "__main__":
	main()























