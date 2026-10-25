import numpy as np
import parsec_coef as pc

xte = 1
yte = 0
xle = 0
yle = 0 

airfoils = {}

#  x_te, y_te, x_c, y_c, dydx_c, dy2dx2_c, te_ang
params_su = [xte, yte, 0.29663, 0.06002, 0, -0.4515, -8.27]
params_pre = [xte, yte, 0.29663, -0.06002, 0, 0.4515, 8.27]

r_le = [0.004, 0.008, 0.0158, 0.024]


for r in r_le:
    c_pre = pc.pcoef(r_le = r, **params_pre, surface='p')
    c_su = pc.pcoef(r_le = r, **params_su, surface='s')

    xx_pre = np.linspace(xte,xle,101)
    yy_pre = (c_pre[0]*xx_pre**(1/2) + 
          c_pre[1]*xx_pre**(3/2) + 
          c_pre[2]*xx_pre**(5/2) + 
          c_pre[3]*xx_pre**(7/2) + 
          c_pre[4]*xx_pre**(9/2) + 
          c_pre[5]*xx_pre**(11/2)
         ) 

    xx_su = np.linspace(xle,xte,101)
    yy_su = (c_su[0]*xx_su**(1/2) + 
          c_su[1]*xx_su**(3/2) + 
          c_su[2]*xx_su**(5/2) + 
          c_su[3]*xx_su**(7/2) + 
          c_su[4]*xx_su**(9/2) + 
          c_su[5]*xx_su**(11/2)
         )

    pre = np.column_stack(xx_pre, yy_pre)
    su = np.column_stack(xx_su, yy_su)

    airfoils[r] = {
        
        'pre' : pre,
        'su' : su,
        'contour' : np.vstack((pre, su[1:]))
        
        }