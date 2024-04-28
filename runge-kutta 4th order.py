

#importing libaries
import numpy as np
from math import sin, cos, pi

#given Lotka-Volterra equations
def equations(t, y):
    alpha = 1.1
    beta = 0.4
    delta = 0.1
    gamma = 0.1
    epsilon = 0.01
    zeta = 0.05
    x = y[0]
    y = y[1]
    x_t = (x * (alpha - (beta * y))) + (x * gamma * sin(2 * np.pi * t))
    y_t = (-y * (delta - (epsilon * x))) + (y * zeta * cos(2 * np.pi * t))
    return np.array([x_t, y_t])

#solves using Runge-Kutta 4th ordee method
def runge_kutta_method(f, t0, tn, y, h):

    #num is the number of iterations, h is step size
    num = int((tn - t0) / h)

    #initializes the list of t and y values
    t_values = [t0 + i * h for i in range(num + 1)]
    y_values = [None] * (num + 1)

    #initializes y and t
    y_values[0] = y
    t = t0

    #calculates slope at each point and updates y and t
    for i in range(num):
        k1 = f(t, y)
        k2 = f(t + h/2, y + (k1 * h)/2)
        k3 = f(t + h/2, y + (k2 * h)/2)
        k4 = f(t + h, y + k3 * h)
        y = y + (h * (k1 + 2*k2 + 2*k3 + k4)) / 6
        t = t + h
        y_values[i + 1] = y
    return t_values, y_values

#main function
def main():

    #initial given values for prey (x(0)) and predator (y(0))
    y0 = np.array([40, 9])

    #chosen step size
    h = 0.1

    #prints the result values
    t_values, y_values = runge_kutta_method(equations, 0, 5, y0, h)
    for i in range(len(t_values)):
        t = t_values[i]
        x = y_values[i][0]
        y = y_values[i][1]
        print(f"At t = {t:.1f}, prey(x(t)) = {x:.4f}, predator(y(t)) = {y:.4f}")

if __name__ == "__main__":
    main()
