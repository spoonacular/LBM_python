import numpy
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


@timeit
def method1():
	x = numpy.arange(1000000)
	x = x*x

@timeit
def method2():
	y = numpy.arange(1000000)
	for i in numpy.arange(1000000):
		y[i] = i*i



method1()
method2()
Foo().foo()



#after roll

# def getMaterialCoords():
# 	materialCoords = numpy.zeros([0,2])
# 	for x in range(5):
# 		for y in range(10):
# 			print(x)
# 			print(y)
# 			print(materialCoords.shape)
# 			print(numpy.array([x,y]))
# 			materialCoord = numpy.append(materialCoords,numpy.array([x,y]),0)
# 	return materialCoords


# materialCoords = getMaterialCoords()
# print(materialCoords)