''' Utils for tidal processing
'''

import numpy as np
import xarray as xr
import dask.array as da
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
#from cmocean import cm
import time
import threading
from utide._ut_constants import ut_constants as utide

g=9.81
cpd=2.*np.pi/86400. # radian/s

#------------------------------ HRET reader -------------------------------------

def get_hret_ssh(constituents=['M2','N2','S2','K1','O1','P1'], lonb=None, latb=None, 
            hret='./Carrere_HRET_testing.nc', bathy=  './ETOPO2v2c_f4.nc'):
    ''' Load HRET ssh
    
    Parameters
    ----------
    consituents: list of str
        Consistuents considered
    lonb: None, list, tuple
        Bounds for longitude, e.g.: lonb=(10.,100.)
    latb: None, list, tuple
        Bounds for latitude, e.g.: latb=(-10.,10.)
    hret: str
        Path to HRET netcdf file
    bathy: str
        Path to Bathymetry netcdffile
    '''
    hret = xr.open_dataset(hret, chunks={'longitude': 500, 'latitude': 500})
    hret_constituents = ['M2','N2','S2','K1','O1','P1']
    for c in hret_constituents:
        if c not in constituents:
            del hret[c+'re']
            del hret[c+'im']        
    #
    omega = dict()
    for cst,o in zip(utide['const']['name'], utide['const']['freq']):
        if cst in constituents:
            omega[cst] = o*24. # cpd, input frequencies are cph
            print(cst+' omega=%e rad/s, %.3f cpd'%(o*2.*np.pi/3600., o*24.))
    #
    if lonb is not None:
        # should handle 360 wrapping
        hret = hret.where(hret['longitude']>=lonb[0], drop=True)
        hret = hret.where(hret['longitude']<=lonb[1], drop=True)
    if latb is not None:
        # should check conventions for longitude
        hret = hret.where(hret['latitude']>=latb[0], drop=True)
        hret = hret.where(hret['latitude']<=latb[1], drop=True)
    #
    if type(bathy) is str:
        h = load_bathy(lon=hret.longitude, lat=hret.latitude , bathy=bathy)
        hret = hret.assign(h=h)
    #
    return hret, constituents, omega


def get_hret_uv(**kwargs):
    ''' Load HRET currents
    
    Parameters
    ----------
    consituents: list of str
        Consistuents considered
    lonb: None, list, tuple
        Bounds for longitude, e.g.: lonb=(10.,100.)
    latb: None, list, tuple
        Bounds for latitude, e.g.: latb=(-10.,10.)
    hret: str
        Path to HRET netcdf file
    '''
    hret, constituents, omega = get_hret_ssh(**kwargs)
    #
    omega_K1 = 15.04107*np.pi/180./3600. # deg/h -> rad/s
    print(omega_K1)
    f = 2*omega_K1*cpd*np.sin(np.pi/180.*hret['latitude'])
    #
    U = xr.Dataset()
    V = xr.Dataset()
    for cst in constituents:
        eta = ri2c(hret[cst+'re'], hret[cst+'im'])
        detadx_c, detady_c = grad(eta)
        U[cst], V[cst] = mom_inv(detadx_c, detady_c, f, omega[cst]*cpd)
    return U, V, constituents, omega

#------------------------------ gradients ---------------------------------------

#
def grad(d, lon='longitude', lat='latitude'):
    """ Compute the gradient of data on lon/lat grid
    """
    dx = di_lonlat(d,lon)
    dy = di_lonlat(d,lat)
    return dx, dy

def di_lonlat(d, c):
    """ Compute the gradient of a variable laid out on a regular lon/lat grid
    """
    if c is 'longitude':
        # treats dateline data correctly
        di = lon_extension(d)
    else:
        di = d   
    #
    dx = di[c].diff(c,label='lower')*111.e3
    di = di.diff(c,label='lower')/dx
    di = (di + di.shift(**{c:1}))*.5
    #
    if c is 'longitude':
        di = di/np.cos(np.pi/180.*di['latitude'])
    return di

#
def lap(d, lon='longitude', lat='latitude'):
    """ Compute the laplacian of data on lon/lat grid
    """
    d2x = di2_lonlat(d,lon)
    d2y = di2_lonlat(d,lat)
    d2x, d2y = xr.align(d2x, d2y, join='outer')
    lap = d2x + d2y
    return lap

def di2_lonlat(d, c):
    """ Compute the second derivative of a variable laid out on a regular lon/lat grid
    """
    if c is 'longitude':
        # treats dateline data correctly
        di = lon_extension(d)
    else:
        di = d   
    #
    dx = di[c].diff(c,label='lower')*111.e3    
    di = di.diff(c,label='lower')/dx
    di = (di - di.shift(**{c:1}))/dx
    #
    if c is 'longitude':
        di = di/np.cos(np.pi/180.*di['latitude'])**2
    return di

#
def lon_extension(v):
    """ Extends data array longitudinally in order to include dateline
    points in a computation of gradients
    """
    if v['longitude'].min()==0. and v['longitude'].max()==359.95:
        v0 = v.sel(longitude=0.)
        v0['longitude']=360.
        return xr.concat([v,v0],dim='longitude')
    else:
        return v
#
def mom_inv(dvdx, dvdy, f, o, r=1./20./86400.):
    """ Inverse a linearized momentum equation
    """
    j = np.complex(0.,1.)
    _o = o + j*r
    uc = -g * (j*_o*dvdx - f*dvdy)/(_o**2-f**2)
    vc = -g * (f*dvdx + j*_o*dvdy)/(_o**2-f**2)
    return uc, vc

def ri2c(r,i):
    return r+np.complex(0,1.)*i

def c2ri(c):
    return np.real(c), np.imag(c)

#------------------------------ bathymetry ---------------------------------------

def load_bathy(lon=None, lat=None, b='etopo2', bathy = './ETOPO2v2c_f4.nc'):
    ''' Load bathymetry
    '''
    #
    if b is 'etopo2':
        #bfile = '/home2/pharos/othr/aponte/bathy/ETOPO2v2c_f4.nc'
        #bfile = './ETOPO2v2c_f4.nc'
        print(bathy)
        hb = -xr.open_dataset(bathy)['z']
    #
    if lon is not None and lat is not None:
        # should use xESMF
        from scipy.interpolate import RectBivariateSpline
        lonb = hb['x'].values
        latb = hb['y'].values
        #
        iroll = np.where(lonb<0)[0][-1]+1
        lonb = np.roll(lonb,-iroll)
        hb = np.roll(hb.values,-iroll,axis=1)
        lonb[lonb<0] = lonb[lonb<0] + 360.
        # should add a test for lon type
        hi = RectBivariateSpline(lonb, latb, hb.T, kx=1, ky=1)(lon,lat).T
        return xr.DataArray(hi, coords={'latitude': lat, 'longitude': lon}, dims=('latitude', 'longitude'))
    else:
        return hb
#-------------------------convert date YYMMdd into number day of the year------------
def date_to_nth_day(date, format='%Y%m%d'):
    import pandas as pd
    date = pd.to_datetime(date, format=format)
    new_year_day = pd.Timestamp(year=date.year, month=1, day=1)
    return (date - new_year_day).days + 1

def sobel_gradient(var, dx, dy):
    #Compute Sobel gradient
    from scipy.signal import convolve2d
    #Sobel Kernel
    sx = np.array([[-1/8,0,1/8],[-1/4,0,1/4],[-1/8,0,1/8]])
    sy = np.array([[1/8,1/4,1/8],[0,0,0],[-1/8,-1/4,-1/8]])
    #
    sobel_x = convolve2d(var, sx, boundary='symm', mode='same')/dx
    sobel_y = convolve2d(var, sy, boundary='symm', mode='same')/dy

    sobel_grad = np.sqrt(sobel_x**2+sobel_y**2)
    return sobel_y, sobel_x
        

def safe_make_folder(i):
    '''Makes a folder if not present'''
    try:  
        os.mkdir(i)
        print('creating: ', i)
    except:
        pass    

def plot_sst(sst, colorbar=False, title=None, vmin=None, vmax=None, savefig=None, offline=False, coast_resolution='110m', figsize=(10,10)):
    if vmin is None:
        vmin = sst.min()
    if vmax is None:
        vmax = sst.max()    
    MPL_LOCK = threading.Lock()
    with MPL_LOCK:
        if offline:
            plt.switch_backend('agg')
        #
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        try:
            im = sst.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax,x='lon', y='lat', add_colorbar=colorbar, cmap=cm.thermal)
            fig.colorbar(im)
            gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='k', 
                    alpha=0.5, linestyle='--')
            gl.xlabels_top = False
            ax.coastlines(resolution=coast_resolution, color='k')
        except:
            pass
        #
        if title is None:
            ax.set_title('HW sst')
        else:
            ax.set_title(title)
        #
        if savefig is not None:
            fig.savefig(savefig, dpi=150)
            plt.close(fig)
        #
        if not offline:
           plt.show()

def get_fes_uv(constituents=['M2','N2','S2','K1','O1','P1'], lonb=None, latb=None, 
            fes='/home2/pharos/othr/aponte/tides/FES2014/', bathy=  './ETOPO2v2c_f4.nc'):
    ''' Load FES currents
    
    Parameters
    ----------
    consituents: list of str
        Consistuents considered
    lonb: None, list, tuple
        Bounds for longitude, e.g.: lonb=(10.,100.)
    latb: None, list, tuple
        Bounds for latitude, e.g.: latb=(-10.,10.)
    fes: str
        Path to FES netcdf file
    '''

    uvdir = fes+'fes2014a_currents/'
    U = xr.Dataset()
    V = xr.Dataset()
    constituents_temp = [c.swapcase() for c in constituents]
    #fes_constituents = ['M2','N2','S2','K1','O1','P1']
    for c in constituents_temp:
        fname_u = uvdir+'eastward_velocity/'+c+'.nc'
        fes_u = xr.open_dataset(fname_u, chunks={'lon': 500, 'lat': 500})
        U[c.swapcase()+'Ua'] = fes_u['Ua']
        U[c.swapcase()+'Ug'] = fes_u['Ug']
        U['longitude'] = fes_u['lon']
        U['latitude'] = fes_u['lat']

        fname_v = uvdir+'northward_velocity/'+c+'.nc'
        fes_v = xr.open_dataset(fname_v, chunks={'lon': 500, 'lat': 500})
        V[c.swapcase()+'Va'] = fes_v['Va']
        V[c.swapcase()+'Vg'] = fes_v['Vg']
        V['longitude'] = fes_v['lon']
        V['latitude'] = fes_v['lat']
   
    #
    omega_K1 = 15.04107*np.pi/180./3600. # deg/h -> rad/s
    print(omega_K1)
    f = 2*omega_K1*cpd*np.sin(np.pi/180.*U['latitude'])
    #

    if lonb is not None:
        # should handle 360 wrapping
        U = U.where(U['longitude']>=lonb[0], drop=True)
        U = U.where(U['longitude']<=lonb[1], drop=True)
        V = V.where(V['longitude']>=lonb[0], drop=True)
        V = V.where(V['longitude']<=lonb[1], drop=True)
    if latb is not None:
        # should check conventions for longitude
        U = U.where(U['latitude']>=latb[0], drop=True)
        U = U.where(U['latitude']<=latb[1], drop=True)
        V = V.where(V['latitude']>=latb[0], drop=True)
        V = V.where(V['latitude']<=latb[1], drop=True)
    if type(bathy) is str:
        h = load_bathy(lon=U.longitude, lat=U.latitude , bathy=bathy)
        U = U.assign(h=h)
        V = V.assign(h=h)
    omega = dict()
    for cst,o in zip(utide['const']['name'], utide['const']['freq']):
        if cst in constituents:
            omega[cst] = o*24. # cpd, input frequencies are cph
            print(cst+' omega=%e rad/s, %.3f cpd'%(o*2.*np.pi/3600., o*24.))
    return U, V, constituents, omega

def make_cartopy(nrow,ncol,projection=ccrs.PlateCarree(), fig=plt.figure(), resolution='100m',nfig=1,title=None,lonticks = None):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import cartopy.feature as cfeature
    import matplotlib.ticker as mticker
    ax = fig.add_subplot(nrow,ncol,nfig,projection=projection)
    ax.coastlines(resolution=resolution, color='k')
    ax.add_feature(cfeature.LAND, facecolor = '0.75')
    
    #ax.gridlines(draw_labels = True)
    if type(title) is str:
        ax.set_title(title)
    gl = ax.gridlines(crs=projection, draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    if lonticks is not None:
        gl.xlabels_bottom = False
        ax.set_xticks(lonticks, crs=ccrs.PlateCarree())
    # Only PlateCarree and Mercator plots are currently supported.
    return ax

def get_slow_fast (v, omega = 2.*np.pi/86400.*1.923, omega2 = None):
    """ Extract the tidal signal with a harmonic analysis 
    The fit looks like: v = v0 + v1 x t + sum( vi_c cos(omega x t) + vi_s sin(omega x t))
    or:  v = v0 + v1 x t + sum( vi_c cos(omega x t) + vi_s sin(omega x t)) + sum( vi_c cos(omega2 x t) + vi_s sin(omega2 x t))
        Parameters
    ----------
        v:
            signal to decompose
        t:
            time line in seconds
        omega:
            tidal frequency in rad/s (M2 = 1.39844e-4)
        omega2:
            diurnal frequency in rad/s (M2/2 = 0,69922e-4)
    """
    t = v.time
    Nt = t.shape
    if (omega2 is not None):
        Xx = [t*0.+1., t , np.cos(omega*t), np.sin(omega*t), np.cos(omega2*t), np.sin(omega2*t)]
        X = np.vstack((np.ones_like(t),t,np.cos(omega*t),np.sin(omega*t),np.cos(omega2*t),np.sin(omega2*t)))
    else:
        X = np.vstack((np.ones_like(t),t,np.cos(omega*t),np.sin(omega*t)))
        Xx = [t*0.+1., t , np.cos(omega*t), np.sin(omega*t)]
    
    X=X.transpose()
    Mx=np.linalg.inv(X.transpose().dot(X))
    
   
    
    XtY = [(x*v*Nt).mean(dim='time') for x in Xx]

    B = []
    for i in range(len(Xx)):
        B.append(XtY[0]*0.)
        for j in range(len(Xx)):
            B[-1] += Mx[i,j]*XtY[j]
    
    vslow = B[0]*Xx[0] + B[1]*Xx[1]
    vfast = B[2]*Xx[2] + B[3]*Xx[3]
    
    return vslow, vfast

def skill(vtrue,v):
    """ compute the skill of a reconstruction (v) of the vtrue signal
    """
   
    dif_vtrue_v = vtrue - v
    skill = 1. - dif_vtrue_v.std('time')/vtrue.std('time')
    
    
    return skill
