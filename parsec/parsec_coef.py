from math import tan, sqrt
import numpy as np

def pcoef(r_le, x_te, y_te, x_c, y_c, dydx_c, dy2dx2_c, te_ang, surface):

    # checks if coef is for suction or pressure surface
    if surface.startswith('p'):
        coef[0]= -sqrt(2*r_le)
    else:
        coef[0]= sqrt(2*r_le)

    #Au = B

    # init all coef
    coef = np.zeros(6)

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
        [y_te - coef[0]*x_te**(1/2)],
        [y_c - coef[0]*x_c**(1/2)],
        [tan(np.radians(te_ang)) - (1/2)*coef[0]],
        [dydx_c - (1/2)*coef[0]*x_c**(-1/2)],
        [dy2dx2_c + (1/4)*x_c**(-3/2)]
    ])

    u = np.linalg.solve(A,B)

    coef[1:6] = u[0:5:0]

    return coef