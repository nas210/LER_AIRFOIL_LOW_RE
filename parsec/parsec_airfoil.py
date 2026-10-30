import os                                        
import numpy as np
import matplotlib.pyplot as plt                  
import parsec_coef as pc
 
xte = 1
yte = 0
xle = 0
yle = 0
 
airfoils_p = {}
 
# x_te, y_te, x_c, y_c, dydx_c, dy2dx2_c, te_ang
params_su = dict(x_te=xte, y_te=yte, x_c=0.29663, y_c= 0.06002, dydx_c=0, dy2dx2_c=-0.4515, te_ang=-8.27)
params_pre = dict(x_te=xte, y_te=yte, x_c=0.29663, y_c=-0.06002, dydx_c=0, dy2dx2_c= 0.4515, te_ang= 8.27)
 
r_le = [0.004, 0.006, 0.008, 0.01, 0.0158, 0.02, 0.024, 0.028]
 
# pin locations and y-values
 
coef_base = pc.pcoef(r_le[4], **params_su, surface= 's')
 
x_p1 = 0.5
x_p2 = 0.7
x_p3 = 0.9
 
 
y_p1 = pc.evaluate(x_p1, coef_base)
y_p2 = pc.evaluate(x_p2, coef_base)
y_p3 = pc.evaluate(x_p3, coef_base)
 
 
pins_su = [(x_p1, y_p1), (x_p2, y_p2), (x_p3, y_p3)]
pins_pre = [(x_p1, -y_p1), (x_p2, -y_p2), (x_p3, -y_p3)]
 
 
n_points = 200
theta = np.linspace(0, np.pi, n_points)
x = 0.5 * (1.0 - np.cos(theta))
 
# output folder
out_dir = 'output'
cfd_dir = os.path.join(out_dir, 'cfd')
dat_dir = os.path.join(out_dir, '.dat format')
os.makedirs(cfd_dir, exist_ok=True)
os.makedirs(dat_dir, exist_ok=True)

# radii sweep
for r in r_le:
    
    c_pre = pc.ppcoef(r_le = r, **params_pre, surface='p', pins = pins_pre)
    c_su = pc.ppcoef(r_le = r, **params_su, surface='s', pins = pins_su)
 
    yy_pre = pc.evaluate(x,c_pre)
    yy_su = pc.evaluate(x,c_su)
 
    pre = np.column_stack((x, yy_pre))
    su = np.column_stack((x, yy_su))
 
    airfoils_p[r] = {
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
 
    # export
    contour = airfoils_p[r]['contour']
 
    # plain x y file
    fname = os.path.join(dat_dir, f"airfoil_rle_{r:.4f}.dat")
    np.savetxt(fname, contour, fmt="%.8f")
 
    rows = []
    for k in range(len(pre)):
        p = pre[::-1][k]
        rows.append([1, p[0], p[1]])
    for k in range(len(su)):
        p = su[k]
        rows.append([2, p[0], p[1]])

    # cad text file -> output/cfd/
    fname_cad = os.path.join(cfd_dir, f"airfoil_rle_{r:.4f}_cad.txt")

    with open(fname_cad, "w") as f:
        f.write("polyline=true\n") 
        
        for surf_id, x_pt, y_pt in rows:
            f.write(f"{surf_id} {x_pt:.8f} {y_pt:.8f}\n")
    
# log
with open(os.path.join(out_dir, "parsec_log.csv"), "w") as f:
    header = ["r_le", "surface"] + [f"a{k+1}" for k in range(9)]
    f.write(",".join(header) + "\n")
    for r in r_le:
        for surf, key in (("suction", "c_su"), ("pressure", "c_pre")):
            c = airfoils_p[r][key]
            f.write(f"{r},{surf}," + ",".join(f"{v:.12e}" for v in c) + "\n")
    f.write("\n")
    f.write(f"pins x: {x_p1},{x_p2},{x_p3}\n")
    f.write(f"pins y (suction): {y_p1:.12f},{y_p2:.12f},{y_p3:.12f}\n")
    f.write(f"baseline r_le: {r_le[4]}\n")
 
# plot
base_su = airfoils_p[r_le[4]]['su']
 
fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.3])
ax0 = fig.add_subplot(gs[0, :])       # full contours (top row)
ax1 = fig.add_subplot(gs[1, 0])       # nose zoom
ax2 = fig.add_subplot(gs[1, 1])       # difference from baseline
 
for r in r_le:
    cont = airfoils_p[r]['contour']
    ax0.plot(cont[:, 0], cont[:, 1], label=f"r_le = {r}")
    ax1.plot(cont[:, 0], cont[:, 1])
    diff = airfoils_p[r]['su'][:, 1] - base_su[:, 1]
    ax2.plot(x, diff * 1e3)
 
ax0.set_title("airfoils_p")
ax0.set_aspect('equal')
ax0.grid(alpha=0.3)
ax0.legend(ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.1))
 
ax1.set_title("Nose Zoom")
ax1.set_xlim(0, 0.12)
ax1.set_ylim(-0.055, 0.055)
ax1.set_aspect('equal')
ax1.grid(alpha=0.3)
 
ax2.set_title("Error In y-Coordinate Along x")
ax2.set_xlabel("x/c")
ax2.set_ylabel("difference (x 1e-3 chord)")
ax2.axhline(0, color='k', linewidth=0.5)
ax2.grid(alpha=0.3)
 
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "parsec_contours.png"), dpi=130)
plt.show()
