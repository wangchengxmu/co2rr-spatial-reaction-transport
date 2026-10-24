"""Editable, idealized geometry for the existing transport calculation.

The copper network is conceptual, not a reconstructed micrograph. Its shape
does not enter the calculation. Dimensions and boundary labels come from the
current model; transverse distances are deliberately not drawn to scale.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MPath
from matplotlib.colors import LinearSegmentedColormap

HERE=Path(__file__).resolve().parent
MM=1/25.4

def arrow(ax,start,end,color='#246C8B',width=1.15,scale=7,style='-|>'):
    a=FancyArrowPatch(start,end,arrowstyle=style,mutation_scale=scale,
                      linewidth=width,color=color,shrinkA=0,shrinkB=0,zorder=8)
    ax.add_patch(a)
    return a

def molecule(ax,x,y):
    # Species icons are illustrative CO2, not quantities or atomistic evidence.
    ax.plot([x-1.45,x+1.45],[y,y],color='#7B6F63',lw=.8,zorder=4)
    for dx,fc,r in [(-1.45,'#C66743',.68),(0,'#545C63',.60),(1.45,'#C66743',.68)]:
        ax.add_patch(Circle((x+dx,y),r,facecolor=fc,edgecolor='white',lw=.35,zorder=5))


def porous_metal(ax):
    """Author-drawn coalesced grains; no measured morphology or pore statistics.

    Smooth grain envelopes and their surface normals provide vector shading.
    The scalar field is a drawing device only and never enters the model.
    """
    x=np.linspace(34.7,74.6,300);y=np.linspace(10.4,29.6,170)
    X,Y=np.meshgrid(x,y);field=np.zeros_like(X)
    grains=[(37,12,3.4,3.6),(44,12.5,4.0,3.0),(52,11.6,4.7,3.5),
            (61,12.4,4.1,3.1),(70,11.8,4.6,3.8),
            (36.2,20.5,3.2,4.2),(44.8,20,3.8,3.4),(54,20,4.1,3.0),
            (63.2,19.4,3.5,3.6),(73,20.4,3.0,4.5),
            (38,28,4.2,3.5),(47.5,27.8,3.7,3.9),(56.3,28.5,4.5,3.6),
            (65,28,4.1,3.2),(73.3,28,3.4,3.5)]
    for cx,cy,rx,ry in grains:
        field+=np.exp(-1.3*((X-cx)**2/rx**2+(Y-cy)**2/ry**2))
    field+=.05*np.sin(1.7*X+.4*Y)*np.cos(1.35*Y)
    # Treat these envelopes as pore cavities in a continuous metal matrix.
    # A narrow, meandering throat joins the central cavities to the liquid.
    throat=20+.7*np.sin((X-35)/4)
    field+=.75*np.exp(-((Y-throat)/.72)**2)
    solid=field<.58
    fy,fx=np.gradient(field,y,x)
    shading=.55-.24*(-.65*fx+.75*fy)/np.sqrt(fx*fx+fy*fy+.10)
    shading+=.08*np.sin(.48*X+.27*Y)
    copper=LinearSegmentedColormap.from_list('copper_surface',
          ['#6F493B','#9D7060','#BF9682','#D5B8A2','#E8D4BE'])
    ax.contourf(X,Y,np.ma.masked_where(~solid,shading),levels=np.linspace(.15,.95,22),
                cmap=copper,vmin=.15,vmax=.95,zorder=2,antialiased=False)
    ax.contour(X,Y,field,levels=[.58],colors=['#815D4D'],linewidths=.48,zorder=3)

def draw_geometry(ax):
    ax.set_xlim(0,168);ax.set_ylim(0,43);ax.set_aspect('equal');ax.axis('off')
    gas_end=34;cat_end=75;outer=157;bottom=10;top=30
    ax.text(1,42,'a',weight='bold',fontsize=11,va='top')
    ax.text(7,42,'Model geometry and transport',weight='bold',fontsize=9.5,va='top')
    ax.text(166,42,'Liquid channel gap: 2 mm',fontsize=8,color='#56646C',ha='right',va='top')
    ax.add_patch(Rectangle((4,bottom),gas_end-4,top-bottom,facecolor='#FFF4E8',edgecolor='none',zorder=0))
    ax.add_patch(Rectangle((gas_end,bottom),outer-gas_end,top-bottom,facecolor='#E3F2F8',edgecolor='none',zorder=0))
    ax.add_patch(Rectangle((outer,bottom),10,top-bottom,facecolor='#CAE5F0',edgecolor='none',zorder=0))
    porous_metal(ax)
    # A conceptual meniscus makes the gas/liquid boundary visible without adding
    # an extra phase or an experimentally determined contact-angle claim.
    verts=[(gas_end,30),(gas_end-1,28),(gas_end+1,26),(gas_end,24),
           (gas_end-1,22),(gas_end+1,20),(gas_end,18),
           (gas_end-1,16),(gas_end+1,13),(gas_end,10)]
    codes=[MPath.MOVETO]+[MPath.CURVE4]*9
    ax.add_patch(PathPatch(MPath(verts,codes),facecolor='none',edgecolor='#2980A1',lw=1.25,zorder=6))
    ax.plot([cat_end,cat_end],[bottom,top],color='#567988',lw=.9,ls=(0,(3,2)),zorder=6)
    ax.plot([outer,outer],[bottom,top],color='#567988',lw=.9,ls=(0,(3,2)),zorder=6)
    for x,y in [(12,25),(22,16),(26,27),(12,13)]:molecule(ax,x,y)
    ax.text(18.5,20,'CO₂ gas',fontsize=9,ha='center',va='center',color='#864526',
            bbox=dict(facecolor='#FFF4E8',edgecolor='none',pad=1))
    ax.text(54,33.5,'Porous Cu',fontsize=9,weight='bold',ha='center',color='#855038')
    ax.text(113,33.5,'Flowing catholyte',fontsize=9,weight='bold',ha='center',color='#1D627F')
    # Keep chemistry labels in open liquid, not on top of solid particles.
    ax.text(104,26.5,'CO₂  ⇌  HCO₃⁻  ⇌  CO₃²⁻',fontsize=8.5,ha='center',color='#24566A')
    ax.text(103,13,'Diffusion + ion migration',fontsize=8,ha='center',color='#24566A')
    arrow(ax,(80,18.5),(118,18.5),style='<->',width=1.0)
    # Flow arrows grow away from the no-slip catalyst wall. Their direction is
    # parallel to the electrode, not through the catalyst.
    for x,height in [(126,5),(135,9),(144,14)]:arrow(ax,(x,12),(x,12+height),width=1.35,scale=8)
    ax.text(139,29,'Flow, x',fontsize=8,ha='center',color='#1D627F')
    ax.text(163,20,'Bulk',rotation=90,fontsize=8.5,ha='center',va='center',color='#24566A')
    # Pore liquid callout is outside the copper drawing, with a precise leader.
    ax.annotate('Pore liquid',xy=(61,12.4),xytext=(55,7),fontsize=8,ha='center',va='center',
                arrowprops=dict(arrowstyle='-',lw=.7,color='#336A80'),color='#24566A')
    ax.annotate('Gas–liquid contact',xy=(34,28),xytext=(20,35.5),ha='center',fontsize=8,
                arrowprops=dict(arrowstyle='-',lw=.7,color='#3D626E'))
    # Dimensions anchor the schematic to the coordinate change in every plot.
    arrow(ax,(34,4.5),(75,4.5),color='#56646C',style='<->',width=.8,scale=6)
    ax.text(54.5,2,'200 nm catalyst',fontsize=8,ha='center',va='center')
    arrow(ax,(75,6),(157,6),color='#56646C',style='<->',width=.8,scale=6)
    ax.text(116,4.5,'500 µm resolved liquid',fontsize=8,ha='center',va='top')
    ax.text(33,7,'z = 0',fontsize=8,ha='right',va='center')
    ax.text(76,8,'y = 0',fontsize=8,ha='left',va='center')

def export_panel():
    plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'svg.fonttype':'none'})
    fig=plt.figure(figsize=(168*MM,43*MM));ax=fig.add_axes([0,0,1,1])
    draw_geometry(ax)
    for ext in ['png','svg']:fig.savefig(HERE/f'panel_a_model_geometry.{ext}',dpi=400,facecolor='white')
    plt.close(fig)

if __name__=='__main__':export_panel()
