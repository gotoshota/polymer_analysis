import os
import subprocess as sp
import glob
import numpy as np
import re
from scipy.optimize import curve_fit

####################################
## Parameters

trjnum = '1'
py_dir = '/lustre/save/users/bze/bending/#' + trjnum + '/'
shapes = ['ring']
#shapes = ['linear']
datadir = '/lustre/save/users/bze/bending/Data/taup_long'
datadir = '/lustre/save/users/bze/bending/Data/taup_scaled_test'
#bendings = ['e1','e2','e4','e8','e16','e32']
bendings = ['bend_1.5', 'bend_5']
#natoms = ['N10','N20','N40','N100','N200','N400']
natoms = ['N100','N200','N400']
rholist = ['0.1', '0.3', '0.5']
dt=0.01
step=100000000
freq=10000
delta = 0.1

prm = step *dt

################################

## Stretched exponential function
def func(x, a, b):
    y = np.exp(-(x/a)**b)
    return y
beta = 1.0

home_dir=os.path.abspath('.')
os.chdir(home_dir)
## Go to Working directory
for shape in shapes:
    for bending in bendings: 
        for natom in natoms:
            for rho in rholist:
                os.chdir(home_dir)
                target_dir = home_dir + '/' + shape + '/' + bending + '/' \
                + natom + '/rho_' + rho
                ## Rename old data to backup
                xvgfile = 'taup.xvg'
                cmd = 'mv ' + xvgfile + ' ' + xvgfile + '.old'
                sp.call(cmd.split())
                ## Get Xp_corre filenames from Xp_long
                
                os.chdir(target_dir)
                filelist = glob.glob('rmf*')
                filelist = sorted(filelist, key=lambda s: int(re.search(r'\d+', s).group()))
                for filename in filelist:
                    ## Read Xp_corre_long
                    os.chdir(py_dir_n)
                    p = filename.replace('rmf_','')
                    p = p.replace('.xvg','')
                    xdata, ydata = np.loadtxt(filename, unpack='True')
                    xdata = np.array(xdata)
                    ydata = np.array(ydata)
                    imax = len(xdata)
                    ## Determin the Initial value of a
                    limit = imax-1
                    for i in range(imax-1):
                        delta = abs(ydata[i] - 0.36)
                        if delta < 0.01:
                            prm = xdata[i]
                        elif ydata[i] <= 10**(-3):
                            limit = i
                            break
                    ## Make list for fitting
                    length = len(xdata)
                    xdata_fit = np.delete(xdata, range(limit-1,length-1))
                    ydata_fit = np.delete(ydata, range(limit-1,length-1))
                    ## Fitting
                    popt, pcov = curve_fit(func, xdata_fit, ydata_fit, p0 = [prm,beta], maxfev=1000000)
                    beta = popt[1] 
                    os.chdir(datadir)
                    n = natom.replace('N','')
                    n = int(n)
                    print(n)
                    rouse = float(p)*np.pi/n
                    rouse = np.sin(rouse)**2
                    with open (xvgfile, mode='a') as f:
                        f.write(str(p) + '  ' + str(popt[0]) + '   ' + str(popt[1]) + '  ' + str(rouse) + '\n')
