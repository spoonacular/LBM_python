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

class bounceBack():
	def __init__(self):
		self.index = 2

class velocityBoundary():
	def __init__(self):
		self.index = 3

class pressureBoundary():
	def __init__(self):
		self.index = 4

class movingWall():
	def __init__(self):
		self.index = 5

	# def calEqulibriumFunction(self,df):

	# 	for i in range(nx):
	# 		for i in range(ny):
	# 			t1 = u(i,j)*u(i,j)+v(i,j)*v(i,j); #u^2
	# 			for k in range(9):
	# 				t2 = u(i,j)*c_xy(1,k)+v(i,j)*c_xy(2,k); %ck * u
	# 				feq(i,j,k) = rho(i,j)*w(k)*(1+3*t2+4.5*t2*t2-1.5*t1);
	# 	            f_new(i,j,k) = omega*feq(i,j,k) + (1-omega)*f(i,j,k);