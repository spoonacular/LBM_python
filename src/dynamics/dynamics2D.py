import numpy

dynamicsIndex = {
	0:'noDynamics',
	1:'BGKdynamics',
	2:'bounceBack (half-way)',
	3:'velocityBoundary',
	4:'pressureBoundary',
	5:'movingWall'
}

class noDynamics():
	def __init__(self):
		self.index = 0

class BGKdynamics():
	def __init__(self,omega):
		self.index = 1
		self.omega = omega #counterclockwise

class BBwall():
	def __init__(self, bounceBackRho = 1):
		self.index = 2
		self.bounceBackRho = bounceBackRho;

class BBvelocity():
	def __init__(self,omega):
		self.index = 3
		self.omega = 1

class BBpressure():
	def __init__(self,omega):
		self.index = 4
		self.omega = 1
class movingWall():
	def __init__(self):
		self.index = 5


class ShanChenBGKdynamics():
	def __init__(self,omega):
		self.index = 1
		self.omega = omega




class ShanChenGenerator():
	def __init__(self,G,interactionPotential):
		self.G = G














	# def calEqulibriumFunction(self,df):

	# 	for i in range(nx):
	# 		for i in range(ny):
	# 			t1 = u(i,j)*u(i,j)+v(i,j)*v(i,j); #u^2
	# 			for k in range(9):
	# 				t2 = u(i,j)*c_xy(1,k)+v(i,j)*c_xy(2,k); %ck * u
	# 				feq(i,j,k) = rho(i,j)*w(k)*(1+3*t2+4.5*t2*t2-1.5*t1);
	# 	            f_new(i,j,k) = omega*feq(i,j,k) + (1-omega)*f(i,j,k);