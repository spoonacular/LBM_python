import sys
import numpy
import matplotlib

def return_specific(x1):
	for x in range(1,x1):
		if x%2==0:
			yield(x)


print(list(return_specific(10)))

#random changes