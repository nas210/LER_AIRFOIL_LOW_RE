import os                                        
import numpy as np
import matplotlib.pyplot as plt                  
import parsec_coef as pc
 
xte = 1
yte = 0
xle = 0
yle = 0
 
airfoils_wo_p = {}
 
# x_te, y_te, x_c, y_c, dydx_c, dy2dx2_c, te_ang
params_su = dict(x_te=xte, y_te=yte, x_c=0.29663, y_c= 0.06002, dydx_c=0, dy2dx2_c=-0.4515, te_ang=-8.27)
params_pre = dict(x_te=xte, y_te=yte, x_c=0.29663, y_c=-0.06002, dydx_c=0, dy2dx2_c= 0.4515, te_ang= 8.27)
 
r_le = [0.004, 0.006, 0.008, 0.01, 0.0158, 0.02, 0.024, 0.028]
 
# pin locations and y-values
 
coef_base = pc.pcoef(r_le[4], **params_su, surface= 's')
 
x_p1 = 0.5
x_p2 = 0.7
x_p3 = 0.9

n_points = 200
theta = np.linspace(0, np.pi, n_points)
x = 0.5 * (1.0 - np.cos(theta))
 
# output folder
out_dir = 'output'
os.makedirs(out_dir, exist_ok=True)
 
# radii sweep
for r in r_le:
    
    c_pre = pc.pcoef(r_le = r, **params_pre, surface='p')
    c_su = pc.pcoef(r_le = r, **params_su, surface='s')
 
    yy_pre = pc.evaluate(x,c_pre)
    yy_su = pc.evaluate(x,c_su)
 
    pre = np.column_stack((x, yy_pre))
    su = np.column_stack((x, yy_su))
 
    airfoils_wo_p[r] = {
        'c_pre': c_pre,
        'c_su': c_su,
        'pre' : pre,
        'su' : su,
        'contour' : np.vstack((pre[::-1], su[1:]))
    }

    x_dense = np.linspace(0, 1, 10001)
    y_dense_su = pc.evaluate(x_dense, c_su)
    y_dense_pre = pc.evaluate(x_dense, c_pre)
    i = np.argmax(y_dense_su)
    print(f"r_le={r}: max y={y_dense_su[i]:.5f} at x={x_dense[i]:.4f}, "
          f"y(TE)={y_dense_su[-1]:.1e}, "
          f"mirror err={np.max(np.abs(y_dense_su + y_dense_pre)):.1e}")
 

# log
with open(os.path.join(out_dir, "parsec_log.csv"), "w") as f:
    header = ["r_le", "surface"] + [f"a{k+1}" for k in range(9)]
    f.write(",".join(header) + "\n")
    for r in r_le:
        for surf, key in (("suction", "c_su"), ("pressure", "c_pre")):
            c = airfoils_wo_p[r][key]
            f.write(f"{r},{surf}," + ",".join(f"{v:.12e}" for v in c) + "\n")
    f.write("\n")
    f.write(f"pins x: {x_p1},{x_p2},{x_p3}\n")
    f.write(f"baseline r_le: {r_le[4]}\n")
 
# error in y wo pins
base_su = airfoils_wo_p[r_le[4]]['su']

fig = plt.figure(figsize=(13, 9))
yerr_wop_ax = fig.add_subplot(111)  

for r in r_le:
    diff = airfoils_wo_p[r]['su'][:, 1] - base_su[:, 1]
    yerr_wop_ax.plot(x, diff * 1e3, label=f"r_le = {r}")

yerr_wop_ax.set_title("Error In y-Coordinate Along x Without Pins")
yerr_wop_ax.set_xlabel("x/c")
yerr_wop_ax.set_ylabel("difference (x 1e-3 chord)")
yerr_wop_ax.axhline(0, color='k', linewidth=0.5)
yerr_wop_ax.grid(alpha=0.3)
yerr_wop_ax.legend(ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.1))

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "error_wo_pins.png"), dpi=130)
