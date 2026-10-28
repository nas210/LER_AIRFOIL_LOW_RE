import numpy as np
import parsec_coef as pc

xte = 1
yte = 0
xle = 0
yle = 0

airfoils = {}

# x_te, y_te, x_c, y_c, dydx_c, dy2dx2_c, te_ang
params_su = dict(x_te=xte, y_te=yte, x_c=0.29663, y_c= 0.06002, dydx_c=0, dy2dx2_c=-0.4515, te_ang=-8.27)
params_pre = dict(x_te=xte, y_te=yte, x_c=0.29663, y_c=-0.06002, dydx_c=0, dy2dx2_c= 0.4515, te_ang= 8.27)

r_le = [0.004, 0.008, 0.0158, 0.024]

# pin locations and y-values

coef_base = pc.pcoef(r_le[2], **params_su, surface= 's')

x_p1 = 0.5
x_p2 = 0.7
x_p3 = 0.9


y_p1 = pc.evaluate(x_p1, coef_base)
y_p2 = pc.evaluate(x_p2, coef_base)
y_p3 = pc.evaluate(x_p3, coef_base)


pins_su = [(x_p1, y_p1), (x_p2, y_p2), (x_p3, y_p3)]
pins_pre = [(x_p1, -y_p1), (x_p2, -y_p2), (x_p3, -y_p3)]

# radii sweep
for r in r_le:
    
    c_pre = pc.ppcoef(r_le = r, **params_pre, surface='p', pins = pins_pre)
    c_su = pc.ppcoef(r_le = r, **params_su, surface='s', pins = pins_su)

    n_points = 200
    theta = np.linspace(0, np.pi, n_points)
    x = 0.5 * (1.0 - np.cos(theta))

    yy_pre = pc.evaluate(x,c_pre)
    yy_su = pc.evaluate(x,c_su)

    pre = np.column_stack((x, yy_pre))
    su = np.column_stack((x, yy_su))

    airfoils[r] = {
        'pre' : pre,
        'su' : su,
        'contour' : np.vstack((pre[::-1], su[1:]))
    }
