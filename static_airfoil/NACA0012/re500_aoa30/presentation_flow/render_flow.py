"""Reproducible 2-D Q/vorticity presentation; run in the flowkit conda env."""
import argparse
import json
from pathlib import Path
import numpy as np
import xarray as xr
from scipy.spatial import Delaunay, cKDTree
from scipy.ndimage import binary_dilation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import Polygon
from matplotlib.animation import FFMpegWriter
from fluidfoam import readmesh, readvector
from flowkit.dataprocessing import Dataset, _rotate

ap = argparse.ArgumentParser()
ap.add_argument('--case', type=Path, default=Path('/home/raji/Research/Thesis/static_airfoil/re500_aoa30'))
ap.add_argument('--out', type=Path, default=Path(__file__).parent)
ap.add_argument('--preview', action='store_true')
args = ap.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
times = sorted((float(p.name), p.name) for p in args.case.iterdir() if p.is_dir() and p.name.replace('.', '', 1).isdigit() and float(p.name)>200)
times = times[-450:]
x,y,_ = readmesh(str(args.case), verbose=False)
dummy = xr.Dataset({'U': (('time','cell','comp'), np.zeros((1,len(x),3)))}, coords={'x':('cell',x), 'y':('cell',y), 'time':[0], 'cell':np.arange(len(x))})
region = Dataset(dummy).rotate(30).crop(-1,8,-2.8,2.8).data
ids = region.cell.values
points = np.column_stack((region.x.values,region.y.values))
tri = Delaunay(points)
h = 0.0125
gx = np.arange(-.75,5.75+h/2,h)
gy = np.arange(-1.9,1.9+h/2,h)
xx,yy = np.meshgrid(gx,gy)
query = np.column_stack((xx.ravel(), yy.ravel()))
simplex = tri.find_simplex(query)
assert (simplex>=0).all()
b = np.einsum('ijk,ik->ij', tri.transform[simplex,:2], query-tri.transform[simplex,2])
weights = np.column_stack((b,1-b.sum(axis=1)))
vertices = tri.simplices[simplex]
bx,by,_ = readmesh(str(args.case), boundary='airfoil', verbose=False)
order = np.argsort(np.arctan2(by, bx-.5))
body = np.column_stack(_rotate(bx[order],by[order],30))
solid = MPath(body).contains_points(query).reshape(xx.shape)
# Remove interpolation stencils that intersect the solid, plus one derivative cell.
crossing = np.zeros(len(query),dtype=bool)
for f in (.25,.5,.75):
    for a,c in ((0,1),(1,2),(2,0)):
        crossing |= MPath(body).contains_points(f*points[vertices[:,a]]+(1-f)*points[vertices[:,c]])
invalid = binary_dilation(solid | crossing.reshape(xx.shape),iterations=1)

# Native weighted least-squares velocity gradients; exclude solid-crossing edges.
dist, neighbors = cKDTree(points).query(points,k=13)
neighbors,dist = neighbors[:,1:],dist[:,1:]
delta = points[neighbors]-points[:,None,:]
lsw = 1/np.maximum(dist,1e-12)**2
for f in (.25,.5,.75):
    crosses = MPath(body).contains_points((points[:,None,:]+f*delta).reshape(-1,2)).reshape(dist.shape)
    lsw[crosses]=0
normal = np.einsum('nki,nk,nkj->nij',delta,lsw,delta)
operator = np.einsum('nij,nkj,nk->nik',np.linalg.pinv(normal),delta,lsw)

def fields(name):
    uv = readvector(str(args.case),name,'U',verbose=False)[:2,ids].T
    u,v = _rotate(uv[:,0],uv[:,1],30)
    velocity = np.column_stack((u,v))
    grad = np.einsum('nik,nkj->nij',operator,velocity[neighbors]-velocity[:,None,:])
    ux,uy,vx,vy = grad[:,0,0],grad[:,1,0],grad[:,0,1],grad[:,1,1]
    native_w = vx-uy
    native_q = -.5*(ux**2+vy**2+2*uy*vx)
    omega = np.einsum('ij,ij->i',weights,native_w[vertices]).reshape(xx.shape)
    q = np.einsum('ij,ij->i',weights,native_q[vertices]).reshape(xx.shape)
    omega[invalid] = np.nan
    q[invalid] = np.nan
    return omega,q

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':15,'axes.labelcolor':'#263343','text.color':'#263343','axes.edgecolor':'#9ba5b0'})
fig = plt.figure(figsize=(16,9),dpi=120,facecolor='white')
ax = fig.add_axes([.075,.19,.85,.66])
fig.text(.075,.925,'Vortex shedding from a stationary airfoil',fontsize=25,weight='normal')
fig.text(.075,.874,r'$Re_c=500\qquad\alpha=30^\circ$',fontsize=18)
stamp = fig.text(.925,.88,'',ha='right',fontsize=18)
cmap = plt.get_cmap('RdBu_r').copy()
cmap.set_bad('white')
im = ax.imshow(np.zeros(xx.shape),extent=(gx[0]-h/2,gx[-1]+h/2,gy[0]-h/2,gy[-1]+h/2),origin='lower',cmap=cmap,vmin=-8,vmax=8,interpolation='nearest')
ax.add_patch(Polygon(body,facecolor='#344354',edgecolor='#172331',lw=.8,zorder=5))
ax.set(xlim=(-.75,5.75),ylim=(-1.9,1.9),xlabel=r'$x_s/c$',ylabel=r'$y_s/c$',aspect='equal')
ax.spines[['top','right']].set_visible(False)
ax.tick_params(length=4,color='#9ba5b0')
ax.annotate('',xy=(.1,1.35),xytext=(-.5,1.35),arrowprops={'arrowstyle':'->','color':'#344354','lw':1.5})
ax.text(-.2,1.52,r'$U_\infty$',ha='center',fontsize=16)
cax = fig.add_axes([.32,.105,.40,.019])
cb = fig.colorbar(im,cax=cax,orientation='horizontal',extend='both',ticks=[-8,-4,0,4,8])
cb.set_label(r'Spanwise vorticity, $\omega_z c/U_\infty$',fontsize=15,labelpad=7)
cb.outline.set_visible(False)
fig.text(.075,.016,r'Vortex regions: $Q c^2/U_\infty^2 > 0.10$   |   Freestream-aligned coordinates',fontsize=12,color='#586576')
contours = None
def update(t,name):
    global contours
    w,q = fields(name)
    im.set_data(np.ma.masked_where((q<=.1)|~np.isfinite(q),w))
    if contours is not None:
        contours.remove()
    contours = ax.contour(xx,yy,q,levels=[.1],colors=['#66717e'],linewidths=.35,alpha=.65)
    stamp.set_text(r'$tU_\infty/c = %.2f$'%t)
    return w,q

sample = times[len(times)//2]
w,q = update(*sample)
fig.savefig(args.out/'re500_aoa30_q_vorticity_preview.png',dpi=120)
print('Preview:',sample,'Q percentiles:',np.nanpercentile(q,[50,90,99]),'vorticity:',np.nanpercentile(w,[1,99]),flush=True)
if not args.preview:
    writer = FFMpegWriter(fps=25,codec='libx264',extra_args=['-crf','18','-pix_fmt','yuv420p','-movflags','+faststart'],metadata={'title':'Re 500, alpha 30: Q-positive regions colored by spanwise vorticity'})
    with writer.saving(fig,str(args.out/'re500_aoa30_q_vorticity.mp4'),dpi=120):
        for i,(t,name) in enumerate(times):
            update(t,name)
            writer.grab_frame()
            if i%50==0: print(f'{i+1}/{len(times)} frames',flush=True)
    (args.out/'metadata.json').write_text(json.dumps({'case':str(args.case),'frames':len(times),'fps':25,'resolution':[1920,1080],'start':times[0][0],'end':times[-1][0],'Q_threshold':.1,'vorticity_limits':[-8,8],'grid_spacing':h,'method':'Flowkit rotation/crop; 12-neighbor inverse-distance weighted least-squares native velocity gradients; linear interpolation of Q and vorticity; body-crossing stencils excluded; Q=-0.5 tr(grad(U)^2).','times':[t for t,_ in times]},indent=2))
    print('Finished',flush=True)
