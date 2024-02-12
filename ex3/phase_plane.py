import numpy as np
from plotly import graph_objects as go
import plotly.figure_factory as ff
import matplotlib.pyplot as plt

def dx(x, y):
    return 3 * x - (2 * x * y) + y

def dy(x, y):
    return x - (y ** 2) + 2

def dx_nullcline(x):
    return (3 * x) / (2 * x - 1)

def dy_nullcline(x):
    return np.sqrt(x + 2)


if __name__ == '__main__':
    x = np.linspace(-4, 4, 25)
    y = np.linspace(-3, 3, 12)

    X, Y = np.meshgrid(x, y)

    Z = dx(X, Y)
    W = dy(X, Y)
    fig, ax = plt.subplots()
    ax.quiver(X, Y, Z, W)

    # add the nullclines
    x_null_x1 = np.linspace(-3, 0.5, 50)
    y_1 = dx_nullcline(x_null_x1)
    x_null_x2 = np.linspace(0.51, 3, 50)
    y_2 = dx_nullcline(x_null_x2)
    # concatenate the two parts, putting none in the middle, so that the line is not connected
    x_null_x = np.concatenate((x_null_x1, [0.5], x_null_x2))
    y_null_y = np.concatenate((y_1, [None], y_2))

    # plot with matplotlib
    ax.plot(x_null_x, y_null_y, label='dx nullcline', color='red')
    ax.plot(x, dy_nullcline(x), label='dy first nullcline', color='royalblue')
    ax.plot(x, -dy_nullcline(x), label='dy second nullcline', color='darkviolet')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Phase Plane of the given system of differential equations')
    ax.legend()
    # set limits for the plot
    ax.set_xlim(-4, 4)
    ax.set_ylim(-3, 3)
    # draw X and Y axes
    ax.axhline(0, color='black', lw=1)
    ax.axvline(0, color='black', lw=1)
    # Add the equilibrium points
    ax.plot(2, 2, 'ro')
    ax.plot(-1, 1, 'ro')
    ax.plot(0.25, -1.5, 'ro')
    # increase the size of the plot, and save it
    fig.set_size_inches(10, 7)
    plt.show()
    fig.savefig('phase_plane.png')



