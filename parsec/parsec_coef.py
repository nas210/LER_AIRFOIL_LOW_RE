# Airfoil coefficients for baseline and for pins

from math import tan, sqrt
import numpy as np

# baseline
def pcoef(r_le, x_te, y_te, x_c, y_c, dydx_c, dy2dx2_c, te_ang, surface):

    # init all coef

    coef = np.zeros(6)

    # checks if coef is for suction or pressure surface
    if surface.startswith('p'):
        coef[0]= -sqrt(2*r_le)

    else:
        coef[0]= sqrt(2*r_le)

    #Au = B

    # unknowns, A
    A = np.array([
        [x_te**(3/2), x_te**(5/2), x_te**(7/2), x_te**(9/2), x_te**(11/2)],
        [x_c**(3/2), x_c**(5/2), x_c**(7/2), x_c**(9/2), x_c**(11/2)],
        [(3/2), (5/2), (7/2), (9/2), (11/2),],
        [(3/2)*x_c**0.5, (5/2)*x_c**(3/2), (7/2)*x_c**(5/2), (9/2)*x_c**(7/2), (11/2)*x_c**(9/2)],
        [(3/4)*x_c**(-1/2), (15/4)*x_c**(1/2), (35/4)*x_c**(3/2), (63/4)*x_c**(5/2), (99/4)*x_c**(7/2)]
        ])

    # knowns, B
    B = np.array([
        y_te - coef[0]*x_te**(1/2),
        y_c - coef[0]*x_c**(1/2),
        tan(np.radians(te_ang)) - (1/2)*coef[0],
        dydx_c - (1/2)*coef[0]*x_c**(-1/2),
        dy2dx2_c + (1/4)*coef[0]*x_c**(-3/2)])

    u = np.linalg.solve(A,B)

    coef[1:6] = u

    return coef

# pins

def ppcoef(r_le, x_te, y_te, x_c, y_c, dydx_c, dy2dx2_c, te_ang, surface, pins = []):

    # checks if coef is for suction or pressure surface

    if surface.startswith('p'):
        a1= -sqrt(2*r_le)
    else:
        a1= sqrt(2*r_le)

    n = 5 + len(pins)

    # exponents list 
    
    exponents = []
    for i in range(n):
        exponents.append(1.5+i)

    A = np.zeros((n,n))
    B = np.zeros(n)

    # y(x_te) = y_te

    for j in range(n):
        p = exponents[j]
        A[0,j] = x_te**(p)
    B[0] = y_te - a1*x_te**(1/2)

    # y(x_c) = y_c

    for j in range(n):
        p = exponents[j]
        A[1,j] = x_c**(p)
    B[1] = y_c - a1*x_c**(1/2)

    # y'(1) = tan(te_ang) - (1/2)a1*x_te^(-1/2)

    for j in range(n):
        p = exponents[j]
        A[2,j] = p*x_te**(p-1)
    B[2] = tan(np.radians(te_ang)) - (1/2)*a1*x_te**(-1/2)

    # dydx_c

    for j in range(n):
        p = exponents[j]
        A[3,j] = p*x_c**(p-1)
    B[3] = dydx_c - (1/2)*a1*x_c**(-1/2)

    # dy2dx2_c

    for j in range(n):
        p = exponents[j]
        A[4,j] = p*(p-1)*x_c**(p-2)
    B[4] = dy2dx2_c + (1/4)*a1*x_c**(-3/2)

    # rows past 5 for pins
    for i in range(len(pins)):
        x_pin = pins[i][0]
        y_pin = pins[i][1]
        row = 5 + i 
        for j in range(n):
            p = exponents[j]
            A[row,j] = x_pin**p
        B[row] = y_pin - a1*x_pin**(1/2)

    u = np.linalg.solve(A,B)

    coef = [a1]
    for co in u:
        coef.append(co)
    return coef

def evaluate(coord, coef):
    exponents = []
    result = 0

    for i in range(len(coef)):
        exponents.append(i+0.5)

    for i in range(len(coef)):
        p = exponents[i]
        result += coef[i]*coord**p

    return result
