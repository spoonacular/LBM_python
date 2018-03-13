class BGKdynamics():
	def __init__(self,omega):
		self.name = 'BGK'
		self.omega = omega
		self.distribution = numpy.array([4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]) #counterclockwise
		self.cx = 			numpy.array([  0,   1,  0, -1,  0,    1,  -1,  -1,   1])
		self.cy = 			numpy.array([  0,   0,  1,  0, -1,    1,   1,  -1,  -1])

	def calEqulibriumFunction():

		for i in range(nx):
			for i in range(ny):
				t1 = u(i,j)*u(i,j)+v(i,j)*v(i,j); #u^2
				for k in range(9):
					t2 = u(i,j)*c_xy(1,k)+v(i,j)*c_xy(2,k); %ck * u
					feq(i,j,k) = rho(i,j)*w(k)*(1+3*t2+4.5*t2*t2-1.5*t1);
		            f_new(i,j,k) = omega*feq(i,j,k) + (1-omega)*f(i,j,k);