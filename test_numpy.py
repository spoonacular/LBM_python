import numpy

y = numpy.array([1,2,3])
x = numpy.zeros([10,5])+1
z = x*numpy.swapaxes([[y]],0,2)

t = numpy.array([[1,2,3],[4,5,6],[7,8,9]])



print(t)
print(numpy.roll(t,-1,1))


w = 1

print('%.2f'.format(w))