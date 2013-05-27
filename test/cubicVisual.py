from visual import sphere,color
from itertools import product

L = 1
R = 0.25

xvals = range(-L, L+1)
yvals = range(-L, L+1)
zvals = range(-L, L+1)

colorfn = lambda *args: [color.yellow, color.green][sum(args)%2]

for pos in product(xvals, yvals, zvals):
    sphere(pos=pos, radius=R, color=colorfn(*pos))
