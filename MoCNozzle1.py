import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
from matplotlib.widgets import Slider, Button, TextBox

# there's like zero dcumentation in this thing
# I'll try to add some at a later time (i.e. not midnight)
# values seem to match up pretty well with JDA example
# (program inits with JDA example values)

## buncha boilerplate functions to handle flow properties, root finding, etc

def bisect(low,high,tol,maxn, func, debug = False):
    if debug:
        print("Bisection Method")
        print(' i\t\ta\t\t\tb\t\t\tp\t\t(b-a)/2\t\t\tf(p)')
    ncnt = 0
    guess = (low + high)/2
    while ncnt < maxn and abs(high-low)/2 >= tol:
        if func(guess) * func(low) > 0:# if guess and low have same sign
            low = guess # then guess is better than low, so move low to guess
        else:
            high = guess # otherwise guess was better than high, etc.
        guess = (low + high)/2
        ncnt += 1
        if debug:
            print('{0:2d}\t{1: 3f}\t{2: 3f}\t{3: 3f}\t{4: 3f}\t\t{5: 3f}'\
            .format(ncnt,low,high,guess,(high-low)/2, func(guess)))
    return guess

def secant(p0,p1,tol,maxn,func, debug = False):
    # Solve nonlinear function using Secant Method
    # ok I've changed this to be an actual secant method
    # instead of a masquerading false position method
    # sorry about that last time!
    if debug:
        print("Secant Method")
        print(' i\t\tp0\t\t\tp1\t\t\tq0\t\t\tq1\t\t\tp\t\tp0-p1')
    ncnt = 0
    p = p0
    plast = 2*(p + tol) # this value will get overwritten, just need to ensure we don't prematurely exit the while loop
    while ncnt < maxn and abs(p - plast) >= tol:
        q0 = func(p0) # implement secant method blah blah
        q1 = func(p1)
        plast = p
        p = p1 - q1 * (p1 - p0) / (q1 - q0)
        p0 = p1
        p1 = p
        ncnt += 1
        if debug:
            print('{0:2d}\t{1:>9.6f}\t{2:>9.6f}\t{3:>9.6f}\t{4:> 9.6f}\t{5:>9.6f}\t{6:>9.6f}'\
            .format(ncnt,p0,p1,q0,q1,p,(p0-p1)))
    return p

def combiBiSecant(low, high, crossover, tol, maxn, func, debug = False):
    # Solve a nonlinear function roughly by using bisection, and then
    # refine with faster convergance with secant method
    # it's just the bisect and secant functions copy/pasted
    # after each other
    if debug:
        print("Combined Bisection-Secant Method")
    tempAns = bisect(low, high, crossover, maxn, func, debug)
    return secant(tempAns - crossover/2, tempAns + crossover/2, tol, maxn, func, debug)

# bisect(0, 10, 0.01, 20, lambda x : np.cos(x) + x - 3, True)
# secant(0, 10, 0.000001, 20, lambda x : np.cos(x) + x - 3, True)
# combiBiSecant(0, 10, 0.1, 0.00001, 20, lambda x : np.cos(x) + x - 3, True)

def goldenSection(low, high, tol, maxn, func, debug = False):
    if debug:
        print("Golden Section Method")
        print(' i\t\tlow\t\t\tmidlow\t\t\tmidhigh\t\thigh\t\t\tf(p)')
    phi = (1+np.sqrt(5))/2
    ncnt = 0
    lowVal = func(low)
    highVal = func(high)
    midHigh = low + (high - low)/phi
    midLow = high - (high - low)/phi
    midLowVal = func(midLow)
    midHighVal = func(midHigh)
    while(abs(high-low)/2 >= tol and ncnt < maxn):
        if debug:
            print('{0:2d}\t{1: 3f}\t{2: 3f}\t{3: 3f}\t{4: 3f}\t\t{5: 3f}'\
            .format(ncnt,low,midLow,midHigh,high, min(midLowVal, midHighVal)))
        ncnt += 1
        if(midHighVal > midLowVal):
            high = midHigh
            highVal = midHighVal
            midHigh = midLow
            midHighVal = midLowVal
            midLow = high - (high - low)/phi
            midLowVal = func(midLow)
        else:
            low = midLow
            lowVal = midLowVal
            midLow = midHigh
            midLowVal = midHighVal
            midHigh = low + (high - low)/phi
            midHighVal = func(midHigh)

    return min(midLow, midHigh)

def findThetaMax(Mach, gamma):
    return goldenSection(np.arcsin(1/Mach), np.radians(90), np.radians(0.01), 40, lambda x : -findTheta(Mach, x, gamma), False)

def findTheta(Mach, beta, gamma):
    if(Mach == np.inf):
        return np.atan((2/np.tan(beta))*(pow(np.sin(beta), 2)) / ((gamma + np.cos(2*beta))))

    return np.atan((2/np.tan(beta))*(pow(Mach*np.sin(beta), 2) - 1) / (Mach*Mach*(gamma + np.cos(2*beta)) + 2))

def findMach(theta, beta, gamma):
    return  np.sqrt((-2/np.tan(beta)-2*np.tan(theta))/((gamma+np.cos(2*beta))*np.tan(theta) - 2/np.tan(beta)*pow(np.sin(beta), 2)))

def findBeta(Mach, beta, gamma, weak, thetaMax):
    if(thetaMax == 0):
        thetaMax = findThetaMax(Mach, gamma)
    low = np.arcsin(1/Mach)
    high = thetaMax
    if weak == 1:
        low = high
        high = np.radians(90)
    ans = combiBiSecant(low, high, 0.1, 0.00001, 20, lambda x : findTheta(Mach, x, gamma) - beta, True)
    return ans

def findBetaExact(Mach, beta, gamma, weak):
    m = Mach
    g = gamma
    weak = -weak + 1
    lamb = np.sqrt(pow(m*m-1,2) - 3*(1+(g-1)/2*m*m)*(1+(g+1)/2*m*m)*pow(np.tan(theta),2))
    chi = (pow(m*m-1,3)-9*(1+(g-1)/2*m*m)*(1+(g-1)/2*m*m+(g+1)/4*m*m*m*m)*pow(np.tan(theta),2))/pow(lamb,3)
    beta = np.atan((m*m-1+2*lamb*np.cos((4*np.pi*weak+np.acos(chi))/3))/(3*(1+(g-1)/2*m*m)*np.tan(theta)))
    return beta

def normp2p1(Mach, gamma):
    return 1 + 2 * gamma / (gamma + 1) * (Mach * Mach - 1)

def normr2r1(Mach, gamma):
    return (gamma + 1) * Mach * Mach / (2 + (gamma - 1) * Mach * Mach)

def normT2T1(Mach, gamma):
    return normp2p1(Mach, gamma) / normr2r1(Mach, gamma)

def normp02p01(Mach, gamma):
    return normp2p1(Mach, gamma) * pow(normT2T1(Mach, gamma), ((-gamma / (gamma - 1))))

def norms2s1(Mach, gamma):
    return -R * np.log(normp02p01(Mach, gamma))

def normp02p1(Mach, gamma):
    return normp02p01(Mach, gamma) * pow((1 + (gamma - 1) / 2 * Mach*Mach), (gamma / (gamma - 1)))

def normalShock(Mach, gamma):
    return np.sqrt((1+(gamma-1)/2*Mach*Mach)/(gamma*Mach*Mach - (gamma-1)/2))

def m2GivenM1Beta(Mach, beta, gamma):
    return normalShock(Mach*np.sin(beta), gamma)/np.sin(beta-findTheta(Mach, beta, gamma))

# print(m2GivenM1Beta(3, np.radians(30), gamma))

def BetaFM1M2(M1, M2, gamma):
    return combiBiSecant(np.asin(1/M1), np.radians(90), 0.01, 0.0001, 20, lambda x : m2GivenM1Beta(M1, x, gamma) - M2, False)

def findNu(Mach, gamma):
    g = gamma
    M = Mach
    findnu1 = np.sqrt((g+1)/(g-1)) * np.atan( np.sqrt((g-1)/(g+1) * (pow(M, 2) - 1))) - np.atan( np.sqrt(pow(M,2) - 1))
    return findnu1

# print(np.degrees(findNu(2.5, 1.4)))

def findPMM(nuM, gamma):
    return combiBiSecant(1, 100, 3, 1e-6, 20, lambda x : findNu(x, gamma) - nuM, False)


# print(findPMM(np.radians(35), 1.4))

def isenttot(M, g):
    return 1 + (g - 1) / 2 * pow(M, 2)

def isenptop(M, g):
    return pow(isenttot(M, g), (g / (g - 1)))

def isenrhotorho(M, g):
    return  pow(isenptop(M, g), (1 / g))

def findMu(Mach):
    return np.asin(1/Mach)

# class Node:
#     # mach and nu can be interchangeable (nu is function of Mach, vice versa)
#     # mu is function of mach (I think)
#     # theta, nu known for first line of points
#     # all edge points same as preceeding point
#     # theta + nu, theta - nu known for all points not on first line (with exception)
#     # points on centerline: know theta and theta + nu at each step
#     def __init__(self, kmin, kmax):
#         self.kmin = kmin
#         self.kmax = kmax

def findIntersectionAngles(x1, y1, phi1, x2, y2, phi2):
    # little geometry function to find the intersection of two lines
    # going through two points. Assumed that point 1 is above point 2,
    # phi 1 is angle of line to horizontal, with positive phi = downwards slope
    # phi 2 is angle of line to horizontal, with positive phi = upwards slope
    epsilon = y1 - y2
    delta = x2 - x1
    delta1 = (epsilon + delta*np.tan(phi2))/(np.tan(phi1)+np.tan(phi2))
    epsilon1 = delta1*np.tan(phi1)
    x3 = x1 + delta1
    y3 = y1 - epsilon1
    return x3, y3

## nozzle design function
# I'm not super happy with how I implemented this
# a much better strategy would have been to use a Node class that acts
# as a linked list to keep track of parent/child nodes, but ok whatever i guess

# currently takes inputs of design exit mach (controls initial slope)
# theta0 in radians (since the first characteristic should, but cannot, be at 0)
# number of charactersitic lines (integer)
# gamma, optional axes for graphing, and option to print tabulated values
# outputs the x and y coordinates of the wall contour

def designNozzleMoC(desMach, theta0, lines, gamma, axs = None, verbose = False):
    thetamax = findNu(desMach, gamma)/2 # calculate initial wall expansion angle
    dtheta = (thetamax - theta0)/(lines - 1) # find spacing between initial characteristics
    # TODO currently the characteristics are linearly spaced by theta, but it might be nicer to try linearly spacing them by mu
    nodes = [] # point number, kmin, kmax, theta, nu, mach, mu
    # big array that holds all the flow properties of the grid points
    # this thing gets... a little messy

    theta = theta0+dtheta*0
    nu = theta0 + dtheta*0
    kmin = theta + nu
    kmax = (theta - nu)
    mach = (findPMM(nu, gamma))
    mu = (findMu(mach))
    nodenum = 1
    # loop to initialize the first set  of grid points
    # all of these points are children of the expansion corner
    for i in range(lines):
        theta = theta0+dtheta*i # know theta and nu for initial points
        nu = theta0 + dtheta*i
        kmin = theta + nu # use that to find kmin, kmax, mach, and mu
        kmax = (theta - nu)
        mach = (findPMM(nu, gamma))
        mu = (findMu(mach))
        nodes.append([nodenum, kmin, kmax, theta, nu, mach, mu]) # put the valus in our big list
        nodenum += 1
    nodes.append([nodenum, kmin, kmax, theta, nu, mach, mu]) # wall points are just a copy of the last point we added
    nodenum += 1

    # this is the part where keeping track of which node we're on gets a little messy
    linesrem = lines
    while linesrem > 1: # the big loop cycle which characteristic line we're on
        linesrem -= 1 # first we start with the point on the centerline
        theta = 0 # perfect reflection across centerline, can assume theta = 0
        nu = nodes[-1-linesrem][1] # take the last point that is above the centerline, both parents have opposite nu so inherit that value
        kmin = theta + nu # find kmin, kmax, etc.
        kmax = theta - nu
        mach = findPMM(nu, gamma)
        mu = findMu(mach)
        indexCenterAnchor = len(nodes) # also it's nice to remember the position of the centerline point for this line
        nodes.append([nodenum, kmin, kmax, theta, nu, mach, mu])
        nodenum += 1
        for i in range(linesrem - 1): # This loop calculates all the point between the centerline and the wall
            kmin = nodes[-1-linesrem][1] # inherit kmin and kmax from parents
            kmax = nodes[indexCenterAnchor][2] # this is the part where list index arithmatic gets hard
            theta = 0.5*(kmin + kmax) # find theta, nu, etc from kmin and kmax
            nu = 0.5*(kmin - kmax)
            mach = findPMM(nu, gamma)
            mu = findMu(mach)
            nodes.append([nodenum, kmin, kmax, theta, nu, mach, mu])
            nodenum += 1
        nodes.append([nodenum, kmin, kmax, theta, nu, mach, mu]) # wall point is copy of last internal point (since no reflection off wall)
        nodenum += 1

    # start finding cartesian positions of all the points
    points = [[0],[1]] # throat is located at (0, 1)

    # find position of first point (centerline), this is just the mu angle of the characteristic
    y3 = 0
    x3 = (points[1][0]-y3) / np.tan(nodes[0][6] - nodes[0][3])
    points[0].append(x3)
    points[1].append(y3)

    # find remaining position of first set of points (where parent is throat expansion corner)
    for i in range(lines-1):
        x1 = points[0][0]
        y1 = points[1][0]
        x2 = points[0][-1]
        y2 = points[1][-1]
        phi1 = nodes[i+1][6] - nodes[i+1][3] # look at theta and mu to find mach lines (intersections are grid point locations)
        phi2 = (nodes[i+1][6] + nodes[i+1][3] + nodes[i][6])/2 # geometry!
        x3, y3 = findIntersectionAngles(x1, y1, phi1, x2, y2, phi2)
        points[0].append(x3)
        points[1].append(y3) # add coordinates to list
    wallslope = (thetamax + nodes[lines][3])/2 # find wall point from average slope of wall and characteristic
    x1 = points[0][0]
    y1 = points[1][0]
    x2 = points[0][-1]
    y2 = points[1][-1]
    phi1 = -wallslope
    phi2 = nodes[lines][6] + nodes[lines][3]
    x3, y3 = findIntersectionAngles(x1, y1, phi1, x2, y2, phi2)
    points[0].append(x3)
    points[1].append(y3)

    # big loop that handles all the points after the initial line
    linesrem = lines - 1
    while linesrem > 0:
        index = len(points[0])
        prevcenter = index - linesrem - 2 # keep track of what the centerline point on the previous line is
        y3 = 0 # add the point on the current line's centerline point
        x3 = (points[1][prevcenter+1]-y3) / np.tan(nodes[prevcenter][6]) + points[0][prevcenter+1]
        currNode = len(points[0]) # also keep track of what the current centerline is
        # use1 = nodes[prevcenter][0] # debug stuff
        # print("calculating {:} using {:}".format(currNode, use1))
        points[0].append(x3)
        points[1].append(y3)
        # find remaining internal points
        for i in range(linesrem-1):
            prevupper = index + i - linesrem - 1
            prevlower = index + i - 1
            # so... keeping track of which points we're referencing is getting harder
            # prevupper and prevlower keep track of the parents of the current node
            x1 = points[0][prevupper+1]
            y1 = points[1][prevupper+1]
            x2 = points[0][-1]
            y2 = points[1][-1]
            phi1 = (nodes[prevupper][6] - nodes[prevupper][3] + nodes[index+i+1][6] - nodes[index+i+1][3])/2
            phi2 = (nodes[prevlower][6] + nodes[prevlower][3] + nodes[index+i+1][6] + nodes[index+i+1][3])/2
            x3, y3 = findIntersectionAngles(x1, y1, phi1, x2, y2, phi2)
            currNode = len(points[0])
            # use1 = nodes[prevupper][0]
            # use2 = nodes[prevlower][0]
            # print("calculating {:} using {:}, {:}".format(currNode, use1, use2))
            points[0].append(x3)
            points[1].append(y3)
        # and then find the wall point using average theta between the current and previous points
        wallslope = (nodes[index - 2][3] + nodes[index + linesrem - 1][3])/2
        x1 = points[0][index - 1]
        y1 = points[1][index - 1]
        x2 = points[0][-1]
        y2 = points[1][-1]
        phi1 = -wallslope
        phi2 = nodes[index + linesrem - 2][6] + nodes[index + linesrem - 2][3]
        x3, y3 = findIntersectionAngles(x1, y1, phi1, x2, y2, phi2)
        points[0].append(x3)
        points[1].append(y3)
        linesrem -= 1

    # plot the locations of the grid points if given axes
    if axs != None:
        axs.axis('equal')
        axs.scatter(points[0],points[1])

    # if printing the tabulated values, convert to degrees first
    if verbose:
        for node in nodes:
            node[1] = np.degrees(node[1])
            node[2] = np.degrees(node[2])
            node[3] = np.degrees(node[3])
            node[4] = np.degrees(node[4])
            node[6] = np.degrees(node[6])
        table = tabulate(nodes, headers = ["#", "K-", "K+", "theta", "nu", "Mach", "mu"])
        print(table)

    # figure out what the characteristic lines are, and plot them
    # some of this probably could've been done in the previous loop but whatever
    centernodes = [] # will be keeping track of list of center nodes
    lastend = 0
    edgex = [points[0][0]]
    edgey = [points[1][0]] # add throat to edge points list
    for i in range(lines):
        # loop through characetristic lines
        xvals = []
        yvals = []
        xvals.append(points[0][0])
        yvals.append(points[1][0])
        for j in range(lines + 1):
            # and then loop through points on charactersitic lines
            # start by traversing across previous lines
            if j < i:
                # print("{:}, {:} ".format(i, j) + str(centernodes))
                xvals.append(points[0][centernodes[j]+i-j])
                yvals.append(points[1][centernodes[j]+i-j])
            elif j == i:
                # then add the centerline point
                centernodes.append(lastend + 1)
                xvals.append(points[0][centernodes[-1]])
                yvals.append(points[1][centernodes[-1]])
            else:
                # and lastly add all the points downstream of the centerline
                # print("{:}, {:} ".format(i, j) + str(centernodes))
                xvals.append(points[0][centernodes[-1]+j-i])
                yvals.append(points[1][centernodes[-1]+j-i])
        lastend = centernodes[-1] + lines - i
        # if plots enabled, add the lines to the plot as well
        if axs != None:
            axs.plot(xvals, yvals)
        # the last point is also an edge, keep track of that
        edgex.append(points[0][lastend])
        edgey.append(points[1][lastend])
    # add the wall line to the plot
    if axs != None:
        axs.plot(edgex, edgey, color='black')

    # print out edge contour
    if verbose:
        table = []
        for i in range(lines):
            table.append([edgex[i], edgey[i]])
        print(tabulate(table, headers = ["x", "y"]))
    # return edge countour
    return [edgex, edgey]

# GUI stuff
fig, axs = plt.subplots(figsize=(8,6))
fig.canvas.manager.set_window_title('MoC nozzle')
fig.tight_layout(pad=3.5)
fig.text(0.5,.94,"Method of Characteristics Nozzle Design", color = 'b',ha='center',fontsize=11)
fig.text(0.01, 0.01, "nmiliv 2024")
plt.subplots_adjust(bottom=0.3)

minmach = 1
maxmach = 10
defmach = findPMM(np.radians(18.375*2),1.4) # JDA uses a rounded design mach in example 11.1, so need to recreate that
mingamma = 1
maxgamma = 2
defgamma = 1.4
gamma = defgamma
mach = defmach
axgamma = plt.axes([0.2, 0.2, 0.6, 0.03])
axmach = plt.axes([0.2, 0.16, 0.6, 0.03])
gammaSlider = Slider(ax=axgamma, label="gamma", valmin=mingamma, valmax=maxgamma, valinit=defgamma, closedmin=False)
machSlider = Slider(ax=axmach, label='Mach', valmin=minmach, valmax=maxmach, valinit=defmach)
axlines = plt.axes([0.2, 0.12, 0.6, 0.03])
linesSlider = Slider(ax=axlines, label="# lines", valmin=1, valmax=30, valinit=7, valstep=1)
axtheta0 = plt.axes([0.2, 0.08, 0.6, 0.03])
theta0Slider = Slider(ax=axtheta0, label="Theta0", valmin=0, valmax=2, valinit=0.375, closedmin=False)
def update(val):
    gamma = gammaSlider.val
    mach = machSlider.val
    lines = linesSlider.val
    theta0 = theta0Slider.val
    axs.clear()
    designNozzleMoC(mach, np.radians(theta0), lines, gamma, axs = axs, verbose = False)
gammaSlider.on_changed(update)
machSlider.on_changed(update)
linesSlider.on_changed(update)
theta0Slider.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
resbutton = Button(resetax, 'Reset', color='gold', hovercolor='skyblue')

def resetSlider(event):
    gammaSlider.reset()
    machSlider.reset()
    linesSlider.reset()
    theta0Slider.reset()

resbutton.on_clicked(resetSlider)

printValUI = plt.axes([0.65, 0.025, 0.13, 0.04])
printbutton = Button(printValUI, 'Print Values', color='gold', hovercolor='skyblue')

def printVals(event):
    gamma = gammaSlider.val
    mach = machSlider.val
    lines = linesSlider.val
    theta0 = theta0Slider.val
    designNozzleMoC(mach, np.radians(theta0), lines, gamma, verbose = True)

printbutton.on_clicked(printVals)


def submitGamma(text):
    if text == "": return
    gamma = float(text)
    gammaSlider.set_val(gamma)
    gamma_box.set_val("")

axgammabox = plt.axes([0.9, 0.2, 0.05, 0.04])
gamma_box = TextBox(axgammabox, "")
gamma_box.on_submit(submitGamma)

def submitMach(text):
    if text == "": return
    mach = float(text)
    machSlider.set_val(mach)
    mach_box.set_val("")

axmachbox = plt.axes([0.9, 0.16, 0.05, 0.04])
mach_box = TextBox(axmachbox, "")
mach_box.on_submit(submitMach)

def submitLines(text):
    if text == "": return
    lines = int(round(float(text)))
    linesSlider.set_val(lines)
    lines_box.set_val("")

axlinesbox = plt.axes([0.9, 0.12, 0.05, 0.04])
lines_box = TextBox(axlinesbox, "")
lines_box.on_submit(submitLines)

def submitTheta0(text):
    if text == "": return
    theta0 = float(text)
    theta0Slider.set_val(theta0)
    theta0_box.set_val("")

axtheta0box = plt.axes([0.9, 0.08, 0.05, 0.04])
theta0_box = TextBox(axtheta0box, "")
theta0_box.on_submit(submitTheta0)

update(0)
# designNozzleMoC(findPMM(np.radians(18.375*2),1.4), np.radians(0.375), 7, 1.4, axs = axs, verbose = True)
plt.show()