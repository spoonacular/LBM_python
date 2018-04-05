import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

x = np.array([[[1,2,3],[4,5,6],[7,8,9]],[[3,2,1],[6,5,4],[7,8,9]]])
y = np.array([[1,0,1],[1,1,1],[0,1,0]])

#x = np.swapaxes(x,0,2)


print(x.shape)
print(y.shape)

z = x*y

print(x.max())