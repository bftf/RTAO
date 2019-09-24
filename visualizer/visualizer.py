
import numpy as np 
import matplotlib.pyplot as plt

plt.switch_backend('TkAgg')

SizeX = 1022
SizeY = 1022

counter = 0;
array = []
cur_line = []

# assume that a framebuffer was generated!
with open('../bin/frameBuffer.txt') as fp:
  for line in fp:
    cur_line.append(int(line.strip('\n')))
    counter += 1

    if counter % SizeX == 0:
      array.append(cur_line)
      cur_line = []
      counter = 0

# print(array)

print('min: ', np.amin(array))
print('max: ', np.amax(array))

# transpose?
# array = [*zip(*array)]

plt.imshow(array, cmap='gray', vmin=0, vmax=32)
plt.show()
