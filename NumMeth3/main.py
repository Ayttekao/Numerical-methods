from sympy import symbols, Matrix
from copy import copy
from typing import List
import numpy as np
from prettytable import PrettyTable

x1, x2, x3 = symbols("x1,x2,x3")

A = Matrix([
    [2.7, 3.3, 1.3],
    [3.5, -1.7, 2.8],
    [4.1, 5.8, -1.7]])

b = Matrix([
    [2.1],
    [1.7],
    [0.8]])

X = Matrix([
    [x1],
    [x2],
    [x3]])

matrix = np.array([
    [2.7, 3.3, 1.3],
    [3.5, -1.7, 2.8],
    [4.1, 5.8, -1.7]])

coefficients = np.array([2.1, 1.7, 0.8], dtype=float)

correct_res = np.linalg.solve(matrix, coefficients)


def gause_method(m: np.ndarray, c: np.ndarray):
    m = np.array(m)
    c = np.array(c)

    for i in range(m.shape[0]):
        for j in range(i + 1, m.shape[0]):
            c[j] -= c[i] * m[j][i] / m[i][i]
            m[j] -= m[j][i] / m[i][i] * m[i]

    for i in range(m.shape[0]):
        c[i] /= m[i][i]
        m[i] /= m[i][i]

    x = np.zeros(m.shape[0])

    for i in range(m.shape[0] - 1, -1, -1):
        x[i] = c[i]
        for j in range(m.shape[0] - 1, i, -1):
            x[i] -= m[i][j] * x[j]
    return x


def steepest_descent_method(m, c, epsilon):
    g: np.array = np.dot(m.transpose(), c)
    B: np.array = np.dot(m.transpose(), m)
    n = 0
    x: List[np.array] = [np.zeros(B.shape[0])]
    r: List = [np.dot(B, x[0]) - g]
    # print(epsilon)
    while n == 0 or not (
            np.sqrt(np.dot(r[-1], r[-1])) < epsilon and np.sqrt(np.dot(x[-1] - x[-2], x[-1] - x[-2])) < epsilon):
        tau = np.dot(r[-1], r[-1]) / np.dot(np.dot(B, r[-1]), r[-1])
        x.append(x[-1] - tau * r[-1])
        r.append(np.dot(B, x[-1]) - g)
        n += 1

    return x[-1], n


def conjugate_gradient_method(A, c, epsilon):
    A = np.array(copy(A))
    g = np.dot(A.transpose(), c)
    A = np.dot(A.transpose(), A)

    x: List[np.array] = [np.zeros(A.shape[0])]
    r: List[np.array] = [g - np.dot(A, x[0])]
    z: List[np.array] = [r[0]]
    n = 0
    while np.sqrt(np.dot(r[-1], r[-1])) / np.sqrt(np.dot(g, g)) > epsilon:
        alpha = np.dot(r[-1], r[-1]) / np.dot(np.dot(A, z[-1]), z[-1])
        x.append(x[-1] + alpha * z[-1])
        r.append(r[-1] - alpha * np.dot(A, z[-1]))
        beta = np.dot(r[-1], r[-1]) / np.dot(r[-2], r[-2])
        z.append(r[-1] + beta * z[-1])
        n += 1

    return x[-1], n


gaussTable = PrettyTable()
gaussTable.field_names = ["x1", "x2", "x3"]
result = gause_method(matrix, coefficients)
gaussTable.add_row([result[0], result[1], result[2]])
gaussTable.float_format = '.7'
print("\t\t\t\t\tGauss Method")
print(gaussTable, "\n")
#print('{0}, {1}, {2}, {3}\n'.format(result[0], result[1], result[2], result[3]))

steepestDescentTable = PrettyTable()
steepestDescentTable.field_names = ["Eps", "x1", "x2", "x3", "Iterations"]
steepestDescentTable.float_format = '.7'

for i in range(2, 8):
    result = steepest_descent_method(matrix, coefficients, pow(10, -i))
    steepestDescentTable.add_row([pow(10, -i), result[0][0], result[0][1], result[0][2], result[1]])
    #print('{0}, {1}, {2}, {3}'.format(result[0][0], result[0][1], result[0][2], result[0][3]))

print("\t\t\t\t\t\tSteepest Descent Method")
print(steepestDescentTable, "\n")


conjugateGradientDescentTable = PrettyTable()
conjugateGradientDescentTable.field_names = ["Eps", "x1", "x2", "x3", "Iterations"]
conjugateGradientDescentTable.float_format = '.7'

for i in range(2, 8):
    result = conjugate_gradient_method(matrix, coefficients, pow(10, -i))
    conjugateGradientDescentTable.add_row([pow(10, -i), result[0][0], result[0][1], result[0][2], result[1]])
    #print(conjugate_gradient_method(matrix, coefficients, pow(10, -i)))

print("\t\t\t\t\t\tConjugate Gradient Method")
print(conjugateGradientDescentTable, "\n")
