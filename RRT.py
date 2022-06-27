import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.spatial as ss
#np.random.seed(1)

def intersect(a1, a2, b1, b2):
    """
    Returns the point of intersection of the lines passing through a2,a1 and b2,b1.
    a1: [x, y] a point on the first line
    a2: [x, y] another point on the first line
    b1: [x, y] a point on the second line
    b2: [x, y] another point on the second line
    """

    s = np.vstack([a1,a2,b1,b2])
    h = np.hstack((s, np.ones((4, 1))))
    l1 = np.cross(h[0], h[1])
    l2 = np.cross(h[2], h[3])
    x, y, z = np.cross(l1, l2)
    print(z)
    if z == 0:
        return False #no intersection
    return True


def RRT(num_points, connectivity, starting_point):
    """ Generates random points using the a version of the Rapid Random Trees algorithmself.
    low connectivity means high spread and vice versa.
    """
    N = num_points
    a = np.multiply(np.ones((N, 2)), starting_point)
    points = np.zeros((N, 2))
    points[0] = starting_point
    for i in range(1, N):
        b = (connectivity*num_points)*np.random.uniform(low = 0, high = 1, size = (1, 2))
        w = np.argmin(np.sum((a-b)**2, 1))
        n = a[w, :]
        k = np.arctan2(b[0][1]-n[1], b[0][0]-n[0])
        c = np.array([n[0]+np.cos(k), n[1]+np.sin(k)])
        a[i, :] = c
        points[i] = (c[0]+n[0])/2, (c[1]+n[1])/2,
    return points




def connect_RRT(points, max_norm):
    connections = [] #store connections [(point1, point1), ..]
    N = len(points[:, 0])
    for i in range(1, N): #skip starting point
        iter = 0
        for j in range(N):
            if i != j:
                if np.linalg.norm(points[i, :]-points[j, :]) < max_norm:
                    iter += 1
                    if iter < 2: #restrict connections in clusters and dont connect both ways
                        connections.append((j, i))

    #exclude nodes without connection
    for i in range(len(points)):
        if (i in np.array(connections)) == False:
            points[i, :] = None

    #TODO
    #exclude nodes with intersect connections

    return points, connections


#testing
starting_point = np.array([0, 0])
points = RRT(30, 0.3, starting_point)
nodes, lines = connect_RRT(points, 1)
plt.scatter(nodes[:, 0], nodes[:, 1])
for i, j in lines:
    plt.plot([nodes[i, 0], nodes[j, 0]], [nodes[i, 1], nodes[j, 1]], 'b-', lw = 0.1)
plt.show()



def create_animation(points):
    fig = plt.figure("RRT")
    ax = fig.add_subplot(111, aspect='equal',
                            xlim=(0, np.max(points[:, 0])), ylim=(0, np.max(points[:, 1])))
    p, = plt.plot(points[0, 0], points[0, 1], 'ko')
    line, = plt.plot([points[0, 0], points[0, 1]], 'k--')

    def animate(i):
        p.set_data(points[i, 0], points[i, 1])
        line.set_data(points[0:i, 0], points[0:i, 1])
        return p, line

    anim = animation.FuncAnimation(fig, animate,  frames=25, \
                                           interval=25, blit=True)
    anim.save('rrt.mp4', fps = 1)
