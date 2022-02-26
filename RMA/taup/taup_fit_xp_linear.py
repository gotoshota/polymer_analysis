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
#shapes = ['ring']
shapes = ['linear']
datadir = '/lustre/save/users/bze/bending/Data/taup_long'
datadir = '/lustre/save/users/bze/bending/Data/taup_scaled_test'
#bendings = ['e1','e2','e4','e8','e16','e32']
bendings = ['e1.5']
#natoms = ['N10','N20','N40','N100','N200','N400']
natoms = ['N100','N200','N400']

dt=0.01
step=100000000
freq=10000
delta = 0.1

delete_ini = 0
delete_mid_1 = 11
delete_mid_2 = 1000
delete_fin = 10001

delete_ini_2 = 100 
delete_fin_2 = 10001

prm = step *dt

################################

## Stretched exponential function
def func(x, a, b):
    y = np.exp(-(x/a)**b)
    return y

## Go to Working directory
for shape in shapes:
    py_dir_s = py_dir + shape
    for bending in bendings: 
        py_dir_e = py_dir_s + '/' + bending
        for natom in natoms:
            ## Rename old data to backup
            xvgfile = 'taup_' + shape + '_' + natom + '_' + bending + '_scaled.xvg'
            os.chdir(datadir)
            cmd = 'mv ' + xvgfile + ' ' + xvgfile + '.bk'
            sp.call(cmd.split())
            ## Get Xp_corre filenames from Xp_long
            py_dir_n = py_dir_e + '/' + natom + '/XVG/Xplong'
            os.chdir(py_dir_n)
            filelist = glob.glob('Xp*')
            filelist = sorted(filelist, key=lambda s: int(re.search(r'\d+', s).group()))
            for filename in filelist:
                ## Read Xp_corre_long
                os.chdir(py_dir_n)
                p = filename.replace('Xp_','')
                p = p.replace('.xvg','')
                xdata_long, ydata_long = np.loadtxt(filename, unpack='True')
                xdata_long = np.array(xdata_long)
                ydata_long = np.array(ydata_long)
                xdata_long = np.delete(xdata_long, 0)
                ydata_long = np.delete(ydata_long, 0)
                ## Read Xp_corre_short
                chdir = py_dir_e + '/' + natom + '/XVG'
                os.chdir(chdir)
                xdata_short, ydata_short = np.loadtxt(filename, unpack='True')
                xdata_short = np.array(xdata_short)
                ydata_short = np.array(ydata_short)
                xdata_short = np.delete(xdata_short, range(delete_ini,delete_mid_1))
                ydata_short = np.delete(ydata_short, range(delete_ini,delete_mid_1))
                xdata_short = np.delete(xdata_short, range(delete_mid_2,delete_fin))
                ydata_short = np.delete(ydata_short, range(delete_mid_2,delete_fin))
                ## Read Xp_corre very short
                xdata_short = np.delete(xdata_short, 0)
                ydata_short = np.delete(ydata_short, 0)
                #chdir = py_dir_e + '/' + natom + '/XVG/Xp.bak'
                chdir = py_dir_e + '/' + natom + '/XVG/Xpshort'
                os.chdir(chdir)
                xdata_veryshort, ydata_veryshort = np.loadtxt(filename, unpack='True')
                xdata_veryshort = np.array(xdata_veryshort)
                ydata_veryshort = np.array(ydata_veryshort)
                xdata_veryshort = np.delete(xdata_veryshort, range(delete_ini_2,delete_fin_2))
                ydata_veryshort = np.delete(ydata_veryshort, range(delete_ini_2,delete_fin_2))
                ## Conbine long and short data
                xdata = np.append(xdata_short, xdata_long)
                ydata = np.append(ydata_short, ydata_long)
                ## And very short
                xdata = np.append(xdata_veryshort, xdata)
                ydata = np.append(ydata_veryshort, ydata)
                ## Write out the combined Xp_corre data file
                dir_combine =  py_dir_e + '/' + natom + '/XVG/Xp_combined'
                #cmd = 'mkdir ' + dir_combine
                #sp.call(cmd.split() )
                os.chdir(dir_combine)
                imax = len(xdata)
                with open (filename, mode='w') as f:
                    for i in range(imax - 1):
                        f.write(str(xdata[i]) + '    ' + str(ydata[i]) + '\n')
                ## Determin the Initial value of a
                for i in range(imax-1):
                    delta = abs(ydata[i] - 0.36)
                    #print(delta)
                    if delta < 0.01:
                        prm = xdata[i]
                    elif ydata[i] <= 10**(-1):
                        limit = i
                        break
                ## Make list for fitting
                length = len(xdata)
                xdata_fit = np.delete(xdata, range(limit-1,length-1))
                ydata_fit = np.delete(ydata, range(limit-1,length-1))
                ## Fitting
                popt, pcov = curve_fit(func, xdata_fit, ydata_fit, p0 = [prm,1.0], maxfev=1000000)
                os.chdir(datadir)
                n = natom.replace('N','')
                n = int(n)
                print(n)
                rouse = float(p)*np.pi/2/n
                rouse = np.sin(rouse)**2
                with open (xvgfile, mode='a') as f:
                    f.write(str(p) + '  ' + str(popt[0]) + '   ' + str(popt[1]) + '  ' + str(rouse) + '\n')
