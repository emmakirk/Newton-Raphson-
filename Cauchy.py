import matplotlib.pyplot as plt
import numpy as np
import math
import argparse
import sys

#-------------------------------------------------
xi = [1.77, -.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44, 3.29, 3.71, -2.40, 4.53, -.07, -1.05, -13.87, -2.53, -1.75, .27, 43.21]
xList = np.arange(-100,100,1)

def cauchyLogLik(xi, theta):
	yList = []
	for theta in xList:
		y = 0
		for x in xi:
			y += (-(math.log(math.pi))-math.log((1+(x-theta)**2)))
		yList.append(y)
	return yList



yList = cauchyLogLik(xi,xList)

start = [-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38]
#--------first derivative----------------------------------------
def d1(xi, theta):
	val = 0
	for x in xi:
		val += (2*(x-theta))/(1+(x-theta)**2)
	return val
#--------second derivative---------------------------------------
def d2(xi,theta):
	val = 0
	for x in xi:
		n = 2*(theta-x-1)*(theta-x+1)
		d = (1+(theta-x)**2)**2
		val += n/d
	return val
#--------newton method-------------------------------------------
def dx(xi,start):
	return start - d1(xi,start)/d2(xi,start)

def newton(xi,start):
	i = 0
	update = -1
	prev = start
	deriv = d2(xi,start)
	while(abs(deriv)>.00001):
		i+=1
		update = dx(xi,prev)
		if abs(update-prev) <= .00001:
			break
		if i > 100:
			break
		deriv = d2(xi,prev)
		prev = update
	return (start, update, i)

for x in start:
	print(newton(xi,x))

#--------figure def----------------------------------------------
figureHeight=4
figureWidth=4
panel1 = plt.figure(figsize=(figureWidth,figureHeight))

#main pannel
panelWidth=2.5
panelHeight=2
relativePanelWidth=panelWidth/figureWidth
relativePanelHeight=panelHeight/figureHeight


#--------plot----------------------------------------------------
#scatter plot
panelx = 0.25
panely = 0.15
panel1=plt.axes([panelx,panely,relativePanelWidth,relativePanelHeight])
panel1.plot(xList, yList)
panel1.set_title('Log Likelihood')



plt.savefig('cauchy.png',dpi=600)
