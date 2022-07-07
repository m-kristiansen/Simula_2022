import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#np.random.seed(5)

def intersect(p1, p2, p3, p4):
    """
    return True if lines intersect
    """

    x1, y1 = p1 # a point on the first line
    x2, y2 = p2 # another point on the first line
    x3, y3 = p3 # a point on the second line
    x4, y4 = p4 # another point on the second line

    #check paralell:
    if (x1-x2)*(y3-y4) == (y1-y2)*(x3-x4):
        return False

    px = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
    py = ((x1*y2-y1*x2)*(y3-y4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))

    #scalar value desciding point of intersection on line 1
    t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)) #scalar value desciding point of intersection on line 1

    #scalar value desciding point of intersection on line 2
    u = ((x1-x3)*(y1-y2)-(y1-y3)*(x1-x2))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))

    _ = [t, u]

    if all(0<= i <= 1 for i in _):
        #print('points = {}'.format(((x1, y1),(x2, y2), (x3, y3), (x4, y4))))
        return True

    return False


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

    #exclude nodes with intersect connections
    for i in range(len(connections)):
        for j in range(len(connections)):
            if i<j:
                if connections[i][0] not in connections[j] and connections[i][1] not in connections[j]:
                    if intersect(points[connections[i][0], :], points[connections[i][1], :], \
                        points[connections[j][0], :], points[connections[j][1], :]):
                        print(connections[i],connections[j])#print intersect connection

    return points, connections


def plot_graph(nodes, lines, view = True):
    fig, ax = plt.subplots()
    ax.scatter(nodes[:, 0], nodes[:, 1], s = 0.1)
    n = range(len(points[:, 0]))
    for i, txt in enumerate(n):
        ax.annotate(txt, (nodes[i, 0], nodes[i, 1]), fontsize = 6)
    for i, j in lines:
        ax.plot([nodes[i, 0], nodes[j, 0]], [nodes[i, 1], nodes[j, 1]], 'b-', lw = 0.1)

    if view:
        plt.show()



starting_point = np.array([0, 0])
points = RRT(20, 0.1, starting_point)
nodes, lines = connect_RRT(points, 1)
plot_graph(nodes, lines)



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
