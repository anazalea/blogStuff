# -*- coding: utf-8 -*-
"""
Created on Mon May 30 13:07:08 2016

@author: anneya

                      .-.
                 .--.(   ).--.
      <-.  .-.-.(.->          )_  .--.
       `-`(     )-'             `)    )
         (o  o  )                `)`-'
        (      )                ,)
        ( ()  )                 )
         `---"\    ,    ,    ,/`
               `--' `--' `--'
                |  |   |   |
                |  |   |   |
                '  |   '   |

"""
from __future__ import print_function, division
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly as py
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
plt.rc('font',family='serif')

# read snapshot data into pandas dataframe 
df = pd.read_csv('RefL0050N0752_Subhalo.csv',header=19)

# include only galaxies more massive than 10**8.5 Msun
massive=df[np.log10(df['MassType_Star'])>8.5]

# generate values to use to assign a random, unique color to each group
groups = massive['GroupID'].unique()
np.random.shuffle(groups)
groupCol = []
for groupid in massive['GroupID']:
    groupCol.append(list(groups).index(groupid))
groupCol = np.array(groupCol)

# calculate Delaunay mesh
points = zip(massive['CentreOfMass_x'],massive['CentreOfMass_y'],massive['CentreOfmass_z'])
dela = Delaunay(points)

# create plotly trace including all delaunay vertices connecting galaxies with sep<1.5Mpc
X = np.array(points)
nebInd,nebIndptr = dela.vertex_neighbor_vertices 
Xdel = []
for k,pt in enumerate(X):
    for n in nebIndptr[nebInd[k]:nebInd[k+1]]:
        d = np.sqrt(np.sum(np.power(X[k]-X[n],2)))
        if d<1.25:
            Xdel.append(X[k])
            Xdel.append(X[n])
            Xdel.append([None,None,None])
Xdel = np.array(Xdel)
traceDel=go.Scatter3d(x=Xdel[:,0],y=Xdel[:,1],z=Xdel[:,2],
                      mode='lines',
                      line=go.Line(color='black', width=3),
                      hoverinfo='none')  
 

# create plotly 3d scatter plot of galaxy positions, marker sizes will reflect
# the mass of each galaxy,                  
trace = go.Scatter3d(
        x = massive['CentreOfMass_x'],
        y = massive['CentreOfMass_y'],
        z = massive['CentreOfmass_z'],
        mode = 'markers',
        marker = dict(
        cmin = 0,
        cmax = len(groups),
        color = groupCol,
        colorscale = 'Rainbow',
        size=(((massive['MassType_Star'])*(3/4.)*(1./np.pi))**(1./3.))/50.,
        line=dict(width=0),
        opacity = 0.75
        ))

# generate plot
layout = go.Layout(
    showlegend=False,
    autosize=False,
    width=900,
    height=900,
    title='EAGLE RefL0050N0752 z=0.37 Galaxy Groups',
    scene = go.Scene(
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0.1, y=2.5, z=0.5)
    ),
    aspectmode='cube',
    xaxis=dict(
    title = "x [Mpc]",
    range = [0,50]
    ),
    yaxis=dict(
    title = "y [Mpc]",
    range=[0,50]
    ),
    zaxis=dict(
    title = "z [Mpc]",
    range=[0,50]
    )),
)
data = [trace,traceDel]
fig = go.Figure(data=data, layout=layout)
py.offline.plot(fig, filename='eagle.html')

