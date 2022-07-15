import sys
path = "/Users/martinkristiansen/Desktop/Simula_2022"
sys.path.insert(0, path+"/Graphs")
from RRT import RRT, connect_RRT
import matplotlib.pyplot as plt

import numpy as np

def getEquidistantPoints(p1, p2, parts):
    new_points = zip(np.linspace(p1[0], p2[0], parts+1)[1: -1], np.linspace(p1[1], p2[1], parts+1)[1: -1])
    new_points = np.asarray(list(new_points))
    return new_points 

def increase_resolution(nodes, lines, num_points):
    x = nodes[:, 0]
    y = nodes[:, 1]
    connections = []
    for i, j in lines:
        new_points = getEquidistantPoints(nodes[i, :], nodes[j, :], num_points+1)
        x = np.insert(x, j, new_points[:, 0])
        y = np.insert(y, j, new_points[:, 1])
    return x,y

starting_point = np.array([0, 0])
points = RRT(2, 1, starting_point)
nodes, lines = connect_RRT(points, 1)
x = nodes[:, 0]
y = nodes[:, 1]
#plt.scatter(points[:, 0], points[:, 1], s = 90, facecolors='none', edgecolors='b')
#plt.scatter(x, y, c = 'r', s = 1)
#plt.show()

partition = [10, 10]

M, xedges, yedges = np.histogram2d(x, y, partition)     #Make 2D image
M = np.where(M == 0 , M, 1)

M_new = np.zeros(partition) #just to start
i = 1
while ((M_new == M).all() == False):
    print(i)
    M = M_new
    x,y = increase_resolution(nodes, lines, i)
    M_new = np.histogram2d(x, y, partition)[0]
    M_new = np.where(M_new == 0 , M_new, 1)
    i+=1

plt.imshow(M_new, origin = "lower")
plt.show()
