from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def plotme2d(wss_init, wss_fin):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel("Criterion 1")
    ax1.set_ylabel("Criterion 2")
    ax1.set_title("Initial clusterings and Pareto clusterings")

    px, py = [], []
    for i,j in wss_init:
        px.append(i)
        py.append(j)

    fx, fy = [], []
    for i,j in wss_fin:
        fx.append(i)
        fy.append(j)

    plt.plot(px, py, 'g^', fx, fy, 'bs')
    plt.show()


def plotme3d(wss_init, wss_fin):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    lx, ly, lz = [], [], []
    for i,j,z in wss_init:
        lx.append(i)
        ly.append(j)
        lz.append(z)

    gx, gy, gz = [], [], []
    for i,j,z in wss_fin:
        gx.append(i)
        gy.append(j)
        gz.append(z)

    ax.scatter(lx, ly, lz, c='b', marker='^')
    ax.scatter(gx, gy, gz, c='r', marker='o')

    ax.set_xlabel('X Criterion 1')
    ax.set_ylabel('Y Criterion 2')
    ax.set_zlabel('Z Criterion 3')
    ax.set_title("Initial clusterings and Pareto clusterings")

    plt.show()








