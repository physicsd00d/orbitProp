import numpy as np
from numpy import arctan, arctan2, arcsin, arccos, sin, cos, sqrt
from numpy import mod, floor, pi
from numpy import array, dot, cross
from numpy.linalg import norm

from plotly.graph_objs import Mesh3d

R_equator = 6378.137;       #[km]
ecc_Earth = 0.081819221456;

def eccentricAnomaly(meanAnom, ecc, tol):
    """ meanAnom: array of mean anomalies to convert to eccentric anomalies
        ecc: the eccentricity of the orbit
        tol: numerical tolerance for error in the eccentricAnomaly
    """
    if (tol <= 0):
        raise('Your tolerance is a little unrealistic, huh?')
    
    eccAnom = np.zeros_like(meanAnom)
    for ix in range(len(eccAnom)):
        en = 0.
        mn = 0.
        err = 1E5
        
        while (err >= tol):
            eccAnom[ix] = en + (meanAnom[ix] - en + ecc*sin(en))/(1. - ecc*cos(en));
            err = abs(eccAnom[ix]-en);
            en = eccAnom[ix];
    
    return eccAnom


def Geodetic_To_ECEF(gdlat, lon, alt):
    """Takes oblate lat, lon, alt and converts to ECEF xyz"""

    Nlat = R_equator/sqrt(1-(ecc_Earth*sin(gdlat))**2);

    x = (Nlat + alt)*cos(gdlat)*cos(lon);
    y = (Nlat + alt)*cos(gdlat)*sin(lon);
    z = (Nlat*(1-ecc_Earth**2) + alt)*sin(gdlat);

    return array([x,y,z])


def ECEF_To_Geodetic(x, y, z):
    """Takes xyz in ECEF and converts to oblate lat,lon,alt"""
    
    # Get longitude and make sure it's not too close to zero
    tol = 1e-10;
    if ((abs(x) <= tol) & (abs(y) <= tol)):
        lon = 0;
    else:
        lon = arctan2(y,x);    #%[deg]
    
    # Initialize and run the iteration to find the geodetic latitude
    err = 50;
    rdelta = sqrt(x**2 + y**2);   
    gdlat = arcsin(z/sqrt(x**2 + y**2 + z**2));

    while( err > tol ):
        Nlat = R_equator/sqrt(1-(ecc_Earth*sin(gdlat))**2);

        gdlat_new = arctan((z+Nlat*sin(gdlat)*ecc_Earth**2)/rdelta);

        err = abs(gdlat_new - gdlat);
        gdlat = gdlat_new;

    alt = rdelta/cos(gdlat) - Nlat;

    return array([gdlat, lon, alt])

def PERI_C_ECI(Rasc, Inc, AoP):
    #radians!
    return np.array([
        [cos(Rasc)*cos(AoP) - cos(Inc)*sin(Rasc)*sin(AoP), 
         cos(AoP)*sin(Rasc) + cos(Rasc)*cos(Inc)*sin(AoP), sin(Inc)*sin(AoP)],
        [-cos(Rasc)*sin(AoP) - cos(Inc)*cos(AoP)*sin(Rasc), 
         cos(Rasc)*cos(Inc)*cos(AoP) - sin(Rasc)*sin(AoP), cos(AoP)*sin(Inc)],
        [sin(Rasc)*sin(Inc)                                 ,
         -cos(Rasc)*sin(Inc)                                , cos(Inc)]
        ])

def UTC_time(day, month, year, hours=0, minutes=0, seconds=0, offset=0):
    """
    % Updated May 28 2012
    % Good for any year
    % offset is in hours
    % output is in days    
    """
    assert(day > 0)
    assert(month > 0)

    RegYear  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    LeapYear = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if (mod(abs(year-1600),4) == 0):
        NumDaysPerMonth = LeapYear;
    else:
        NumDaysPerMonth = RegYear;

    UTC = 0;
    ix = 1;
    while ix < month:
        UTC = UTC + NumDaysPerMonth[ix-1]; # -1 b/c index from zero
        ix = ix+1;

    if day <= NumDaysPerMonth[month-1]:    # -1 b/c index from zero
        UTC = UTC + day
    else:
        raise('Not that many days in the month')

    return UTC + (hours+offset)/24 + minutes/(24*60) + seconds/(60*60*24);

def TrueAnomaly(EccAnom, MeanAnom, ecc):
    
    num = cos(EccAnom) - ecc;
    den = 1 - ecc*cos(EccAnom);

    half_plane = floor(MeanAnom/pi);  #which half plane are we in?

    TrueAnom = half_plane*pi + arccos((-1)**mod(half_plane,2) * num/den);

    return mod(TrueAnom,2*pi)

# For plotting orbits in plotly
def ellipsoidPlot(radius, nPts = 20):
    phi = np.linspace(0, 2*pi, nPts)
    theta = np.linspace(-pi/2, pi/2, nPts)
    phi, theta=np.meshgrid(phi, theta)

    x = cos(theta) * sin(phi) * radius
    y = cos(theta) * cos(phi) * radius
    z = sin(theta) * radius
    
    return Mesh3d({
                'x': x.flatten(), 
                'y': y.flatten(), 
                'z': z.flatten(), 
                'alphahull': 0})

def getOrbitalElements(r, v):
    """canonical units and column vectors please"""

    h = cross(r,v);

    p = norm(h)**2;

    ecc_vec = (norm(v)**2 - 1/norm(r))*r - dot(r,v)*v;

    ecc = norm(ecc_vec);

    if (ecc == 0):
        raise('Do Something about circular orbit')

    a = p/(1-ecc**2);

    meanmotion = (1/a)**(1.5);   #%rad/TU

    inc = arccos(h[2]/norm(h));

    if (inc == 0):
        raise('Do Something about equatorial orbit')
    elif (inc == pi):
        raise('Do Something about polar orbit')
    
    n_vec = cross(array([0, 0, 1]),h);

    raan = arccos(n_vec[0]/norm(n_vec));
    if (n_vec[1] < 0):
        raan = 2*pi - raan;

    aop = arccos(dot(n_vec,ecc_vec)/(norm(n_vec)*ecc));
    if (ecc_vec[2] < 0):
        aop = 2*pi - aop;

    nu0 = arccos(dot(ecc_vec,r)/(ecc*norm(r)));
    if (dot(r,v) < 0):
        nu0 = 2*pi - nu0;

    num = ecc + cos(nu0);
    den = 1 + ecc*cos(nu0);
    EccAnom0 = arccos(num/den);

    M0 = EccAnom0 - ecc*sin(EccAnom0);

    return np.array([a, ecc, inc, raan, aop, nu0, meanmotion, M0])

def siderealTime(UTC):
    #% Earth physical constants from Vallado
    rotEarthRad = 0.0000729211585530;  #% [rad/sec]
    rotEarthDay = 1.0027379093;  #%rev/day
    rotEarthDay = rotEarthRad*3600*24/(2*pi);

    #% Greenwich Sidereal Time at 0000h 1 January 2012 UTC
    #% from US Naval Observatory website
    gst2012start = 6.6706801; #% [sidereal hours] from vernal equinox

    #%convert to rad
    gst2012start = gst2012start * 2*pi/24;    #%[rad]

    #% Find angle of Greenwich meridian from vernal equinox [0, 2*pi)
    #% theta_g = mod(gst2012start + rotEarthDay*(UTC - 1)*2*pi,2*pi);
    theta_g = mod(gst2012start + rotEarthRad*(UTC - 1)*3600*24,2*pi);

    return theta_g



