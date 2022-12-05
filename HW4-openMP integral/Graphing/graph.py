# importing the required module
from cProfile import label
from ctypes import sizeof
from tabnanny import filename_only
import matplotlib.pyplot as plt
import os

paralle = []
linear = []

for filename in os.listdir(r'C:\\Users\\masoud-pc\Documents\\Parallel-Programming\\HW4-openMP integral\\Graphing\\Inputs\\'):
    with open(os.path.join(r'C:\\Users\\masoud-pc\Documents\\Parallel-Programming\\HW4-openMP integral\\Graphing\\Inputs\\', filename), 'r') as f:
        res = f.read()
        paralle.append([res.splitlines()[0], res.splitlines()[1]])

paralle.sort(key=lambda x: float(x[0]))

speedup = []
core = []

for i in range(1):
    for j in paralle:
        core.append(float(j[0]))
        speedup.append(float(j[1]))

    # x axis values
    x = core
    # corresponding y axis values
    y = speedup
    # plotting the points
    plt.plot(x, y, label='N = 10^6', marker='o', markersize=5)
    plt.legend()
    speedup = []
    core = []

plt.xticks(x, [1, 2, 3, 4, 8, 12, 16])
plt.locator_params(axis='x', nbins=14)

plt.ylim(0, 16)
plt.locator_params(axis='y', nbins=14)

plt.axline((1, 1), (2, 2), linestyle='dashed', color="0", label="S = P")
plt.legend()

# naming the x axis
plt.xlabel('number of threads')
# naming the y axis
plt.ylabel('Speedup')

# giving a title to my graph
plt.title('Speedup graph for solving integral')

# function to show the plot
plt.show()
