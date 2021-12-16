from scipy.signal import convolve2d, fftconvolve
import numpy as np
from numpy.fft import fft2, ifft2

np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})

f = open("input/input.txt", "r")
ah, aw, bh, bw = [int(x) for x in next(f).split()]
i = 0
matrix = []
kernel = []
line = f.readline()
while i < ah:
    matrix.append(list(map(float, line.split())))
    line = f.readline()
    i += 1
i = 0
while i < bh:
    kernel.append(list(map(float, line.split())))
    line = f.readline()
    i += 1
matrix = np.matrix(matrix)
kernel = np.matrix(kernel)

conv = fftconvolve(matrix, kernel, mode='full')
np.savetxt('outputs/python_matrix.txt', conv, fmt='%.2f')
# print(conv)
