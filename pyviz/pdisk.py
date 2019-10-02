import numpy as np
import math
import matplotlib as mpl 
import matplotlib.tri as tri 
import matplotlib.pyplot as plt 
import configparser
import pylab 
from scipy import ndimage 
from scipy import integrate 
from scipy import interpolate 

#simulation hardwired parameters (unlikely to change)
nghost = 3 
r0     = 1.0 
bigG   = 1.0
Mstar  = 1.0 
omega0 = np.sqrt(bigG*Mstar/r0**3.0)
period0= 2.0*np.pi/omega0 

#plottig parameters 
fontsize= 24
nlev    = 128
nclev   = 6 
cmap    = plt.cm.inferno
nrad_lim   = 128
ntheta_lim = 64

#read sim parameters, by first converting to ini file, then use config parser
svar =  open("variables.par","r").read()
svar =  '[params]\n' + svar
svar = svar.replace("\t"," = \t")
file=open("variables.ini","w")
file.write(svar)
file.close()

config = configparser.ConfigParser()
config.read('variables.ini')
smallh = float(config['params']['ASPECTRATIO'])
flare  = float(config['params']['FLARINGINDEX'])
sigma0 = float(config['params']['SIGMA0'])
epsilon= float(config['params']['EPSILON'])
smallq = 3.0 - 2.0*(flare + 1.0);
NX     = int(config['params']['NX'])
NY     = int(config['params']['NY'])
NZ     = int(config['params']['NZ'])

#get co-ordinates (azimuth, radius, theta) 
data = np.loadtxt("domain_x.dat")
nx   = data.size
nphi = nx - 1
azi  = np.empty([nphi])
for i in range(0,nphi):
    beg = i
    end = beg+2
    azi[i] = np.mean(data[beg:end]) 
        #python loops/ranges don't include the last index, so this is equiv to data[i:i+1] in IDL 
    
data = np.loadtxt("domain_y.dat") 
ny   = data.size
nrad = ny - 2*nghost - 1
rad  = np.empty([nrad])
for i in range(0,nrad):
    beg = i+3 
    end = beg+2
    rad[i] = np.mean(data[beg:end]) 
rmin = data[nghost]/r0
rmax = data[ny - nghost-1]/r0
 
data   = np.loadtxt("domain_z.dat")
nz     = data.size
ntheta = nz - 2*nghost - 1
theta  = np.empty([ntheta])
for i in range(0,ntheta):
    beg = i+3 
    end = beg+2
    theta[i] = np.mean(data[beg:end]) 
thetamin = data[nghost]
thetamax = data[nz - nghost-1]


#get time axis 
data   = np.loadtxt("planet0.dat")
time   = data[:,8]

#get corresponding cylindrical domain 
#atm just take nRad = nrad and nZ = ntheta, but may need larger nZ to resolve inner disk 
Rmin  = rmin
Rmax  = rmax*np.sin(thetamin)
zmin  = rmax*np.cos(thetamax)
zmax  = rmax*np.cos(thetamin)
nRad  = nrad
nzcyl = ntheta*2
Raxis = np.linspace(Rmin, Rmax, nRad)
zaxis = np.linspace(zmin, zmax, nzcyl)

def rphi(var='dg',
             loc     = './', 
             zslice  = 0.5, 
             start   = 0,
             log     = None,
             pert    = None,
             title   = '',
    ):
      
    #prepare to read data 
    
    fname = loc+"gasdens"+str(start)+".dat"
    dens  = pylab.fromfile(fname).reshape(NZ,NY,NX) 

    fname = loc+"gasvy"+str(start)+".dat"
    vrad  = pylab.fromfile(fname).reshape(NZ,NY,NX)

    fname   = loc+"gasvz"+str(start)+".dat"
    vtheta  = pylab.fromfile(fname).reshape(NZ,NY,NX)

    if(var == 'vz'):
        data3d = get_vz(rad, theta, vrad, vtheta)
        title = '$v_{gz}/c_s$'

    if(var == 'dg'):
        fname = loc+"dust1dens"+str(start)+".dat"
        densd  = pylab.fromfile(fname).reshape(NZ,NY,NX)
        data3d = densd/dens
        title = r'$\rho_\mathrm{d}/\rho_\mathrm{g}$'
        if(log != None):
            data3d = np.log10(data3d)
            title = r'$\log{(\rho_\mathrm{d}/\rho_\mathrm{g})}$'

    if(var == 'gas'):
        data3d = dens 
        if(pert != None):
           fname  = loc+"gasdens0.dat"
           dens0  = pylab.fromfile(fname).reshape(NZ,NY,NX)
           data3d /= dens0
           title = r'$\rho_\mathrm{g}/\rho_\mathrm{g,i}$'
           if(log != None):
              data3d  = np.log10(data3d) 
              title = r'$\log{(\rho_\mathrm{g}/\rho_\mathrm{g,i})}$'
           else:
              data3d -= 1.0 
              title = r'$\Delta\rho_\mathrm{g}/\rho_\mathrm{g,i}$'
        else:
#            npl     = np.argmin(np.absolute(rad - r0))
#            rho_ref = dens0[NZ/2,npl,0]              
           rho_ref = sigma0/np.sqrt(2.0*np.pi)/(smallh*r0)
           data3d /= rho_ref
           if(log != None):
              data3d = np.log10(data3d)
              title = r'$\log{(\rho_\mathrm{g}/\rho_\mathrm{g,ref})}$'     
           else:
              data3d-=1.0
              title = r'$\Delta\rho_\mathrm{g}/\rho_\mathrm{g,ref}$'    
       
    if(var == 'dust'):
        fname = loc+"dust1dens"+str(start)+".dat"
        densd  = pylab.fromfile(fname).reshape(NZ,NY,NX)
        data3d = densd 
 
        if(pert != None):
           fname  = loc+"dust1dens0.dat"
           densd0  = pylab.fromfile(fname).reshape(NZ,NY,NX)
           data3d /= densd0
           title = r'$\rho_\mathrm{d}/\rho_\mathrm{d,i}$'
           if(log != None):
              data3d  = np.log10(data3d)
              title = r'$\log{(\rho_\mathrm{d}/\rho_\mathrm{d,i})}$'
           else:
              data3d -= 1.0
              title = r'$\Delta\rho_\mathrm{d}/\rho_\mathrm{d,i}$'
        else:
#            npl     = np.argmin(np.absolute(rad - r0))
#            rho_ref = dens0[NZ/2,npl,0]
           rho_ref  = sigma0/np.sqrt(2.0*np.pi)/(smallh*r0)
           rhod_ref = rho_ref*epsilon   
           data3d /= rhod_ref
           if(log != None):
              data3d = np.log10(data3d)
              title = r'$\log{(\rho_\mathrm{d}/\rho_\mathrm{d,ref})}$'
           else:
              data3d-=1.0
              title = r'$\Delta\rho_\mathrm{d}/\rho_\mathrm{d,ref}$'


    nzslice = int(zslice*ntheta)
    data2d  = data3d[nzslice,...].transpose()
    

    #plot data 

    plt.figure(figsize=(7,9))
    plt.subplots_adjust(left=0.2, right=0.9, top=0.95, bottom=0.1)

    minv = np.amin(data2d)
    maxv = np.amax(data2d)
    levels  = np.linspace(minv,maxv,nlev)
    clevels = np.linspace(minv,maxv,nclev)

    plt.rc('font',size=fontsize,weight='bold')
    plt.rc('axes',labelsize=fontsize)
    plt.ylim(0.0,1.0)
    plt.xlim(rmin, rmax)
    
    cp = plt.contourf(rad/r0, azi/(2.0*np.pi), data2d,
                      levels,
                      cmap=cmap
                      )

    plt.colorbar(cp,ticks=clevels,format='%.3f')
    plt.title(title)
    plt.xlabel('$r/r_0$')
    plt.ylabel('$\phi/2\pi$')

    fname = var+'2d_'+str(start).zfill(4)
    plt.savefig(fname,dpi=150)


    #polar contour plot 
    angles = np.repeat(azi[...,np.newaxis], nrad, axis=1)
    xaxis = (rad*np.cos(angles)).flatten()
    yaxis = (rad*np.sin(angles)).flatten()
    data2d_flat = data2d.flatten()

    triang = tri.Triangulation(xaxis, yaxis)

    xmid = xaxis[triang.triangles].mean(axis=1)
    ymid = yaxis[triang.triangles].mean(axis=1)
    mask = np.where(xmid*xmid + ymid*ymid < rmin*rmin, 1, 0)
    triang.set_mask(mask)
  
    plt.figure(figsize=(8,6))
    plt.gca().set_aspect('equal')
    plt.xlim(-rmax, rmax)

    cp = plt.tricontourf(triang, data2d_flat, 
                         levels,
                         cmap=cmap
                         )

    plt.colorbar(cp,ticks=clevels,format='%.3f')
    plt.title(title)
    plt.xlabel('$x/r_0$')
    plt.ylabel('$y/r_0$')
    
    fname = var+'xy_'+str(start).zfill(4)
    plt.savefig(fname,dpi=150)

    return


def get_cs(bigR):
    cs0 = smallh*np.sqrt(bigG*Mstar/r0)
    cs  = cs0*np.power(bigR/r0, -smallq/2.0)
    return cs

def get_vz(rad, 
           theta, 
           vrad, 
           vtheta,
           ):
    vz = np.zeros([ntheta,nrad,nphi])
    for j in range(0,nrad):
        for i in range(0,ntheta):
            bigR = rad[j]*np.sin(theta[i])
            cs   = get_cs(bigR)
            vz[i,j,:]   = vrad[i,j,:]*np.cos(theta[i])-vtheta[i,j,:]*np.sin(theta[i])
            vz[i,j,:]  /= cs #normalize vert velo by sound-speed
            
    return vz


def rz(var='gas',
             loc       = './', 
             azislice  = 0.0, 
             start     = 0,
             log       = None,
             pert      = None,
             plotrange = None,
             lowres    = None 
    ):
      

    fname = loc+"gasdens"+str(start)+".dat"
    dens  = pylab.fromfile(fname).reshape(NZ,NY,NX) 

    fname = loc+"gasvy"+str(start)+".dat"
    vrad  = pylab.fromfile(fname).reshape(NZ,NY,NX)

    fname   = loc+"gasvz"+str(start)+".dat"
    vtheta  = pylab.fromfile(fname).reshape(NZ,NY,NX)

    if(var == 'vz'):
        data3d = get_vz(rad, theta, vrad, vtheta)
        title = '$v_{gz}/c_s$'

    if(var == 'dg'):
        fname = loc+"dust1dens"+str(start)+".dat"
        densd  = pylab.fromfile(fname).reshape(NZ,NY,NX)
        data3d = densd/dens
        title = r'$\rho_\mathrm{d}/\rho_\mathrm{g}$'
        if(log != None):
            data3d = np.log10(data3d)
            title = r'$\log{(\rho_\mathrm{d}/\rho_\mathrm{g})}$'

    if(var == 'gas'):
        data3d = dens 
        if(pert != None):
           fname  = loc+"gasdens0.dat"
           dens0  = pylab.fromfile(fname).reshape(NZ,NY,NX)
           data3d /= dens0
           title = r'$\rho_\mathrm{g}/\rho_\mathrm{g,i}$'
           if(log != None):
              data3d  = np.log10(data3d) 
              title = r'$\log{(\rho_\mathrm{g}/\rho_\mathrm{g,i})}$'
           else:
              data3d -= 1.0 
              title = r'$\Delta\rho_\mathrm{g}/\rho_\mathrm{g,i}$'
        else:
#            npl     = np.argmin(np.absolute(rad - r0))
#            rho_ref = dens0[NZ/2,npl,0]              
           rho_ref = sigma0/np.sqrt(2.0*np.pi)/(smallh*r0)
           data3d /= rho_ref
           if(log != None):
              data3d = np.log10(data3d)
              title = r'$\log{(\rho_\mathrm{g}/\rho_\mathrm{g,ref})}$'     
           else:
              data3d-=1.0
              title = r'$\Delta\rho_\mathrm{g}/\rho_\mathrm{g,ref}$'    
       
    if(var == 'dust'):
        fname = loc+"dust1dens"+str(start)+".dat"
        densd  = pylab.fromfile(fname).reshape(NZ,NY,NX)
        data3d = densd 
 
        if(pert != None):
           fname  = loc+"dust1dens0.dat"
           densd0  = pylab.fromfile(fname).reshape(NZ,NY,NX)
           data3d /= densd0
           title = r'$\rho_\mathrm{d}/\rho_\mathrm{d,i}$'
           if(log != None):
              data3d  = np.log10(data3d)
              title = r'$\log{(\rho_\mathrm{d}/\rho_\mathrm{d,i})}$'
           else:
              data3d -= 1.0
              title = r'$\Delta\rho_\mathrm{d}/\rho_\mathrm{d,i}$'
        else:
#            npl     = np.argmin(np.absolute(rad - r0))
#            rho_ref = dens0[NZ/2,npl,0]
           rho_ref  = sigma0/np.sqrt(2.0*np.pi)/(smallh*r0)
           rhod_ref = rho_ref*epsilon   
           data3d /= rhod_ref
           if(log != None):
              data3d = np.log10(data3d)
              title = r'$\log{(\rho_\mathrm{d}/\rho_\mathrm{d,ref})}$'
           else:
              data3d-=1.0
              title = r'$\Delta\rho_\mathrm{d}/\rho_\mathrm{d,ref}$'


 

    tslice = time[start]/period0
    tstring = "{:.0f}".format(tslice)
    title +=', t='+tstring+r'$P_0$'

    #downsize data & coords  

    if(lowres != None):
        nrad_small   = np.amin([nrad, nrad_lim])
        ntheta_small = np.amin([ntheta, ntheta_lim])

        zoom   = np.empty(3)
        zoom   = [ntheta_small/ntheta, nrad_small/nrad, 1.0]
        
        rad_plot    = ndimage.zoom(rad, zoom[1])
        theta_plot  = ndimage.zoom(theta, zoom[0])
        nrad_plot   = rad_plot.size

        data3d = ndimage.zoom(data3d,zoom)
    else:
        rad_plot   = rad
        theta_plot = theta
        nrad_plot  = nrad 

    #take slice in azimuth or perform azi average
    if(azislice >= 0.0):
    	nazislice = int(azislice*nphi)
    	data2d  = data3d[...,nazislice]
    else:
        data2d  = np.average(data3d,axis=2)

    #polar contour plot 
    angles = np.repeat(np.pi/2.0-theta_plot[...,np.newaxis], nrad_plot, axis=1)
    xaxis = (rad_plot*np.cos(angles)).flatten()
    yaxis = (rad_plot*np.sin(angles)).flatten()
    data2d_flat = data2d.flatten()

    triang = tri.Triangulation(xaxis, yaxis)

    xmid = xaxis[triang.triangles].mean(axis=1)
    ymid = yaxis[triang.triangles].mean(axis=1)
    mask = np.where(xmid*xmid + ymid*ymid < rmin*rmin, 1, 0)
    triang.set_mask(mask)
  
    plt.figure(figsize=(9,4.5))
    plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.18)

    ymin = np.amin(yaxis)
    ymax = np.amax(yaxis)

    plt.ylim(ymin,ymax)
    plt.xlim(rmin, rmax)

    if(plotrange == None):
         minv = np.amin(data2d)
         maxv = np.amax(data2d)
    else:
         minv=plotrange[0]
         maxv=plotrange[1] 


    levels  = np.linspace(minv,maxv,nlev)
    clevels = np.linspace(minv,maxv,nclev)


    plt.rc('font',size=fontsize,weight='bold')
#    plt.rc('axes',labelsize=fontsize,labelweight='bold')



    cp = plt.tricontourf(triang, data2d_flat, 
                         levels,
                         cmap=cmap
                         )
    
    plt.colorbar(cp,ticks=clevels,format='%.2f')
    plt.title(title,weight='bold')

    plt.xticks(fontsize=fontsize,weight='bold')
    plt.xlabel('$R/r_0$',fontsize=fontsize)

    plt.yticks(fontsize=fontsize,weight='bold')
    plt.ylabel('$z/r_0$',fontsize=fontsize)
 
    fname = var+'rz_'+str(start).zfill(4)
    plt.savefig(fname,dpi=150)

    return


def cyl_to_sph_coords(z, R, phi):
    psi = math.atan2(z, R)
    t   = np.pi/2.0 - psi 
    r   = np.sqrt(z*z + R*R)
    p   = phi 
    return (t,r,p)

def sph_to_cyl_regrid(data_spherical):

    data_cylindrical = np.zeros([nzcyl, nRad, nphi])

    for k in range(0,nphi):
        phi = azi[k]
        data2d = data_spherical[...,k]
        f = interpolate.interp2d(rad, theta, data2d, kind='cubic')
        for j in range(0, nRad):
            R = Raxis[j]
            for i in range(0, nzcyl):
                z = zaxis[i]        
                t, r, p = cyl_to_sph_coords(z, R, phi)
                if(thetamin <= t <= thetamax):
                    data_cylindrical[i,j,k] = f(r, t)
    
    return data_cylindrical
    
def get_surfdens(data_spherical):
    data_vintegrated = np.zeros([nRad,nphi])

    data_cylindrical = sph_to_cyl_regrid(data_spherical)

    for j in range(0, nRad):
        for k  in range(0, nphi):
            data_vintegrated[j,k] = integrate.simps(data_cylindrical[:,j,k], zaxis)

    return data_vintegrated

'''
def test(start = 0,
         loc = './',
     ):

    fname = loc+"dust1dens"+str(start)+".dat"
    rhod  = pylab.fromfile(fname).reshape(NZ,NY,NX)
    
    rhod_cylindrical = sph_to_cyl_regrid(rhod)
    rhod_cylindrical_axi = np.mean(rhod_cylindrical, axis=2)

    fname = loc+"gasdens"+str(start)+".dat"
    rhog  = pylab.fromfile(fname).reshape(NZ,NY,NX)
    
    rhog_cylindrical = sph_to_cyl_regrid(rhog)
    rhog_cylindrical_axi = np.mean(rhog_cylindrical, axis=2)
  
    dg  = np.zeros([nzcyl, nRad])
    res = np.divide(rhod_cylindrical_axi, rhog_cylindrical_axi, out=dg, where=rhog_cylindrical_axi>0.0)

    cp = plt.contourf(Raxis, zaxis, dg)
    plt.colorbar(cp)

    plt.show()
'''

def metal(loc      = './', 
                start    = 0, 
                title    = '',
                plotrange=None, 
     ):

    '''
    measure the vertically-integrated dust-to-gas ratio
    we interpolate to a cylindrical grid, then integrate vertically 
    '''

    fname = loc+"gasdens"+str(start)+".dat"
    rhog  = pylab.fromfile(fname).reshape(NZ,NY,NX) 

    fname = loc+"dust1dens"+str(start)+".dat"
    rhod  = pylab.fromfile(fname).reshape(NZ,NY,NX) 

    sigma_d = get_surfdens(rhod)
    sigma_g = get_surfdens(rhog)

    sigd = np.mean(sigma_d,axis=1)
    sigg = np.mean(sigma_g,axis=1)

    metal = sigd/sigg
    
    plt.figure(figsize=(8,4.5))
    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    if(plotrange == None):
        ymin = np.amin(metal)
        ymax = np.amax(metal)
    else:
        ymin = plotrange[0]
        ymax = plotrange[1]

    plt.ylim(ymin,ymax)
    plt.xlim(rmin,rmax)

    plt.plot(rad,metal,linewidth=2)

    plt.rc('font',size=fontsize,weight='bold')

    plt.title(title,weight='bold')

    plt.xticks(fontsize=fontsize,weight='bold')
    plt.xlabel('$R/r_0$',fontsize=fontsize)
    
    plt.yticks(fontsize=fontsize,weight='bold')
    plt.ylabel('$\Sigma_d/\Sigma_g$',fontsize=fontsize)
 
    fname = 'metal_'+str(start).zfill(4)
    plt.savefig(fname,dpi=150)
    
    return

def dgmid(loc      = './', 
                start    = 0, 
                title    = '',
                plotrange=None, 
     ):

    '''
    measure the dust-to-gas ratio at the disk midplane 
    '''

    fname = loc+"gasdens"+str(start)+".dat"
    rhog  = pylab.fromfile(fname).reshape(NZ,NY,NX) 

    fname = loc+"dust1dens"+str(start)+".dat"
    rhod  = pylab.fromfile(fname).reshape(NZ,NY,NX) 

    t1    = np.argmin(np.absolute(theta - np.pi/2.0))

    rhod_mid = np.mean(rhod[t1,...],axis=1)
    rhog_mid = np.mean(rhog[t1,...],axis=1)

    dg = rhod_mid/rhog_mid
        
    
    plt.figure(figsize=(8,4.5))
    plt.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    
    if(plotrange == None):
        ymin = np.amin(dg)
        ymax = np.amax(dg)
    else:
        ymin = plotrange[0]
        ymax = plotrange[1]

    plt.ylim(ymin,ymax)
    plt.xlim(rmin,rmax)

    
    plt.plot(rad,dg,linewidth=2)

    plt.rc('font',size=fontsize,weight='bold')

    plt.title(title,weight='bold')

    
    plt.xticks(fontsize=fontsize,weight='bold')
    plt.xlabel('$R/r_0$',fontsize=fontsize)
    
    
    plt.yticks(fontsize=fontsize,weight='bold')
    plt.ylabel(r'$\rho_{d0}/\rho_{g0}$',fontsize=fontsize)
 
    
    fname = 'dgmid_'+str(start).zfill(4)
    plt.savefig(fname,dpi=150)
    

    return
