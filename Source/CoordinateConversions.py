# SUPPORTING MODULE FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# by Josh Dillon

import numpy as np
import healpy as hp

#takes in a list with rows [RA, dec] in radians
def convertHorizToEqua(horizAngles):
	"testing"


test = np.arange(0,np.pi,.1)
#print test
#print np.cos(test)

print np.asarray(hp.pix2ang(16,np.arange(0,16*16*12))[:,0]
# 	horizPoint toHoriz(double LST){
# 		double lha = pi/12.0*LST - ra; //radians
# 		horizPoint hp;
#     	hp.alt = asin(sin(latInRad) * sin(dec) + cos(latInRad) * cos(dec) * cos(lha));
#     	hp.az = atan2(sin(lha) * cos(dec), cos(lha) * cos(dec) * sin(latInRad) - sin(dec) * cos(latInRad)) + pi;
#     	return hp;
# 	}



# struct horizPoint; //Predeclaration
# struct equaPoint; //Predeclaration

# //When appropriate, +x is east, +y is north,  +z is up
# struct cartVec{
# 	double x,y,z;	
# 	cartVec(){}
# 	cartVec(double xIn, double yIn, double zIn) : x(xIn), y(yIn), z(zIn) {}
# 	cartVec operator+(const cartVec& that){ return cartVec(x+that.x, y+that.y, z+that.z); }
# 	cartVec operator-(const cartVec& that){ return cartVec(x-that.x, y-that.y, z-that.z); }
# 	cartVec operator*(const double that){ return cartVec(x*that, y*that, z*that); }
# 	double dot(cartVec that){
# 		return (x*that.x + y*that.y + z*that.z);
# 	}
# 	cartVec cross(cartVec that){
# 		return cartVec(y*that.z - z*that.y, z*that.x - x*that.z, x*that.y - y*that.x);
# 	}
# 	cartVec normalize(){
# 		return (cartVec(x,y,z) * (1.0 / sqrt(x*x + y*y + z*z)));
# 	}
# 	horizPoint toHoriz();
# 	void print(){ cout << "["<< x << ", " << y << ", " << z << "]"; }
# };

# //Alt is 0 at the horizon, pi/2 at the zenith. Az is radians east of north
# struct horizPoint{
# 	double alt, az;
# 	horizPoint(){}
# 	horizPoint(double altIn, double azIn) : alt(altIn), az(azIn) {}
# 	cartVec toVec(){
# 		cartVec cCoord(sin(az)*cos(alt), cos(az)*cos(alt), sin(alt));
# 		return cCoord;
# 	}
# 	double greatCircleAngle(horizPoint hp){
# 		horizPoint thisHP(alt, az);
# 		cartVec here = thisHP.toVec();
# 		cartVec there = hp.toVec();
# 		double dist = sqrt(pow(here.x - there.x,2) + pow(here.y - there.y,2) + pow(here.z - there.z,2));
# 		return 2*asin(dist/2);
# 	}
# 	equaPoint toEqua(double LST);
# };

# //Converts cartesian vectors to altitude and azimuth. I needed to delcate it outside the struct because horizPoint hadn't been declared yet.
# horizPoint cartVec::toHoriz(){
# 	return horizPoint(asin(z/sqrt(x*x + y*y + z*z)), fmod(atan2(x,y)+2*pi,2*pi)); 
# }

# //Creates an equatorial pointing (RA and Dec) which can be converted into a horizontal pointing given an LST (assuming an array latitude as a global variable)
# struct equaPoint{
# 	double ra, dec;
# 	equaPoint(){}
# 	equaPoint(double raIn, double decIn) : ra(raIn), dec(decIn) {}
# 	horizPoint toHoriz(double LST){
# 		double lha = pi/12.0*LST - ra; //radians
# 		horizPoint hp;
#     	hp.alt = asin(sin(latInRad) * sin(dec) + cos(latInRad) * cos(dec) * cos(lha));
#     	hp.az = atan2(sin(lha) * cos(dec), cos(lha) * cos(dec) * sin(latInRad) - sin(dec) * cos(latInRad)) + pi;
#     	return hp;
# 	}
# 	pointing toHealpixPointing(){
# 		double phi = fmod(ra-facetRAinRad+pi, 2*pi);
# 		double theta = pi/2-dec+facetDecInRad;
# 		if ((theta < 0) || theta > pi) phi = fmod(phi + pi, 2*pi);
# 		if (theta < 0) theta = abs(theta);
# 		if (theta > pi) theta = (pi - (theta-pi));
# 		return pointing(theta, phi); //pointing is defined as the colatitude (0 to pi = N to S) and the longitude (0 to 2pi)
# 	}
# };

# //Converts cartesian vectors to altitude and azimuth. I needed to delcate it outside the struct because horizPoint hadn't been declared yet.
# equaPoint horizPoint::toEqua(double LST){
# 	equaPoint ep;
# 	ep.dec = asin(sin(alt)*sin(latInRad) - cos(alt)*cos(latInRad)*cos(az-pi));
# 	ep.ra = LST*2.0*pi/24 - atan2(sin(az-pi),(cos(az-pi)*sin(latInRad)+tan(alt)*cos(latInRad))); //The one I think is right
# 	return ep;
# }

# //Converts a pointing on the healpix sphere (which as been rotated so that the facet is centered on the equator) to an equatorial pointing
# equaPoint convertToEquaPoint(pointing& hpPoint){
# 	double RA = fmod(hpPoint.phi + facetRAinRad - pi, 2*pi);
# 	double dec = pi/2 - hpPoint.theta + facetDecInRad;
# 	if ((dec < -pi/2) || dec > pi/2) RA = fmod(RA + pi, 2*pi);
# 	if (dec < -pi/2) dec = -pi/2 + (-pi/2 - dec);
# 	if (dec > pi/2) dec = pi/2 - (dec - pi/2);
# 	return equaPoint(RA,dec);
# }

# //Stucture for keeping track of point sources
# struct pointSource{
# 	double flux, spectralIndex;
# 	equaPoint point;
# 	pointSource(){}
# 	pointSource(equaPoint pointIn, double fluxIn, double spectralIndexIn) : point(pointIn), flux(fluxIn), spectralIndex(spectralIndexIn) {}
# 	double fluxAt(double freq){ 
# 		return flux * pow(freq/pointSourceReferenceFreq, -spectralIndex);
# 	}
# };

# //Converts galactic coordinates to equatorial
# equaPoint convertGalToEquaPoint(pointing& galPointing){
# 	double b = pi/2.0 - galPointing.theta;
# 	double l = galPointing.phi;
# 	double pole_ra = 2.0*pi/360.0*192.859508;
#     double pole_dec = 2.0*pi/360.0*27.128336;
#     double posangle = 2.0*pi/360.0*(122.932-90.0); //32.932
#     double raHere = atan2((cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle))) + pole_ra;
#     double decHere = asin(cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec));
#     return equaPoint(raHere, decHere);
# }

# //Converts equatorial coordinates to galactic
# pointing convertEquaPointToGal(equaPoint& eq){
# 	double pole_ra = 2.0*pi/360.0*192.859508;
#     double pole_dec = 2.0*pi/360.0*27.128336;
#     double posangle = 2.0*pi/360.0*(122.932-90.0);
# 	double b = asin( sin(eq.dec)*sin(pole_dec) + cos(eq.dec)*cos(pole_dec)*cos(pole_ra - eq.ra) );
# 	double l = (posangle+pi/2) - atan2(cos(eq.dec)*sin(eq.ra-pole_ra), cos(pole_dec)*sin(eq.dec) - sin(pole_dec)*cos(eq.dec)*cos(eq.ra-pole_ra));	
# 	return pointing(pi/2-b, l);
# }