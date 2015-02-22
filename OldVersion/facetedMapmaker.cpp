
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include "fftw3.h"
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <healpix_base.h>
#include "healpix_map.h" 
#include "healpix_data_io.h" 
#include "gsl/gsl_fit.h"
#include "healpix_map_fitsio.h" 
#include "gsl/gsl_fit.h"
#include <gsl/gsl_linalg.h>
#include <map>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SETTINGS AND GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

string specsFilename = "Specifications.txt";

//Settings to load
double freq; //MHz
string polarization;
string baselinesFile, visibilitiesFolder, trueSkyFilename, finalMapFilename, PSFfilename, healpixPixelFilename, extendedHealpixPixelFilename, noiseCovFilename, DmatrixFilename, dataProductFilePrefix, pixelCoordinatesFilename, extendedPixelCoordiantesFilename;
double angularResolution; //Degrees
double arrayLat; //Degrees
double arrayLong; //Degrees
double facetRA; //Degrees
double facetDec; //Degrees;
double integrationTime; //seconds
double channelBandwidth; //MHz
double snapshotTime; //seconds
int NSIDE; //Sets HEALPIX resolution
double facetSize; //degrees
double PSFextensionBeyondFacetFactor; 
double maximumAllowedAngleFromPBCenterToFacetCenter; //degrees
double xpolOrientationDegreesEastofNorth;
double noiseStd; //Jy on each visibilility...this number is totally made up
bool gaussianPrimaryBeam; //if true, use guassian PB at zenith, else use MWA dipoles (ala omniscope)
double primaryBeamFWHM; //FWHM in degrees if using a guassian PB
bool overwriteVisibilitiesWithGSM;
int GSMres; //Resolution for simulating GSM visibilities
bool overwriteVisibilitiesWithPointSources; 
bool alsoComputePointSourcePSF; //This is used to compute the PSF at specific points
string pointSourceFile; //contains ras, decs, fluxes, and spectral indicies
double pointSourceReferenceFreq; //the frequency (in MHz) at which the flux is defined
double pointSourceFluxUpperLimit; //Only model point sources whose fluxes (times the beam at the center of the observation) are above this
double pointSourceFluxLowerLimit; //Only model point sources whose fluxes (times the beam at the center of the observation) are below this
string pointSourceOutputCatalogFilename;
string pointSourcePSFFilename;
bool throwOutFluxOutsideExtendedFacet; //if overwriting with the GSM and/or point sources, only include sources inside the extended facet


//Hard-coded settings and global variables
bool PSFpeaksAtOne = false; //Otherwise, the diagonal of the PSF matrix is set to be all ones
const double pi = 3.1415926535897932384626433832795;
const double c = 299792458;
int nAltBeamPoints = 45; //running from 0 to 90 in increments of 90/44
int nAzBeamPoints = 180; //running for 0 to 360 in increments of 360/179
double firstBeamFreq = 110; //MHz
double beamDeltaFreq = 10; //MHz
double lastBeamFreq = 190; //MHz
string beamDirectory = "../MWA_Primary_Beams/mwa_beam_";
int nFreqFiles = int((lastBeamFreq - firstBeamFreq)/10 + 1);
double deltaBeamAlt = 90.0 / (nAltBeamPoints-1);
double deltaBeamAz = 360.0 / (nAzBeamPoints-1);

//Hardcoded GSM variables
string FITS_directory = "../VisibilitySimulator/FITS_GSM/";
string GSMformat1 = "gsm";
string GSMformat2 = ".fits";
int componentsToFit = 9;

//Global variables to compute later
int nIntegrationsToUse = 0;
int nPixels, nPixelsExtended, nBaselines, nLSTs, nSnapshots, nFull;
int nPointSources = 1;
double k, facetDecInRad, facetRAinRad, angResInRad, latInRad, temperatureConversionFactor;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// USEFUL DATA STRUCTURES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Stores 2D coordinates a and d
struct coord{
	int i, j;
	coord(int iIn, int jIn) : i(iIn), j(jIn) {}
};

//Stores complex numbers (e.g. visibilities). Can multiply them by real (double precision) or other complex numbers
struct complex{
	double re, im;
	complex(){}
	complex(double reIn, double imIn) : re(reIn), im(imIn) {}
	double abs(){ return sqrt(re*re + im*im); }
	complex conj(){ return complex(re, -im); }
	complex operator*(const complex& that){ return complex(re*that.re - im*that.im, re*that.im + im*that.re); }
	complex operator*(const double that){ return complex(re*that, im*that); }
	complex operator+(const complex& that){ return complex(re + that.re, im + that.im); }
	complex operator-(const complex& that){ return complex(re - that.re, im - that.im); }
	void print(){ cout << re << " + " << im << "i"; }
};


struct horizPoint; //Predeclaration
struct equaPoint; //Predeclaration

//When appropriate, +x is east, +y is north,  +z is up
struct cartVec{
	double x,y,z;	
	cartVec(){}
	cartVec(double xIn, double yIn, double zIn) : x(xIn), y(yIn), z(zIn) {}
	cartVec operator+(const cartVec& that){ return cartVec(x+that.x, y+that.y, z+that.z); }
	cartVec operator-(const cartVec& that){ return cartVec(x-that.x, y-that.y, z-that.z); }
	cartVec operator*(const double that){ return cartVec(x*that, y*that, z*that); }
	double dot(cartVec that){
		return (x*that.x + y*that.y + z*that.z);
	}
	cartVec cross(cartVec that){
		return cartVec(y*that.z - z*that.y, z*that.x - x*that.z, x*that.y - y*that.x);
	}
	cartVec normalize(){
		return (cartVec(x,y,z) * (1.0 / sqrt(x*x + y*y + z*z)));
	}
	horizPoint toHoriz();
	void print(){ cout << "["<< x << ", " << y << ", " << z << "]"; }
};

//Alt is 0 at the horizon, pi/2 at the zenith. Az is radians east of north
struct horizPoint{
	double alt, az;
	horizPoint(){}
	horizPoint(double altIn, double azIn) : alt(altIn), az(azIn) {}
	cartVec toVec(){
		cartVec cCoord(sin(az)*cos(alt), cos(az)*cos(alt), sin(alt));
		return cCoord;
	}
	double greatCircleAngle(horizPoint hp){
		horizPoint thisHP(alt, az);
		cartVec here = thisHP.toVec();
		cartVec there = hp.toVec();
		double dist = sqrt(pow(here.x - there.x,2) + pow(here.y - there.y,2) + pow(here.z - there.z,2));
		return 2*asin(dist/2);
	}
	equaPoint toEqua(double LST);
};

//Converts cartesian vectors to altitude and azimuth. I needed to delcate it outside the struct because horizPoint hadn't been declared yet.
horizPoint cartVec::toHoriz(){
	return horizPoint(asin(z/sqrt(x*x + y*y + z*z)), fmod(atan2(x,y)+2*pi,2*pi)); 
}

//Creates an equatorial pointing (RA and Dec) which can be converted into a horizontal pointing given an LST (assuming an array latitude as a global variable)
struct equaPoint{
	double ra, dec;
	equaPoint(){}
	equaPoint(double raIn, double decIn) : ra(raIn), dec(decIn) {}
	horizPoint toHoriz(double LST){
		double lha = pi/12.0*LST - ra; //radians
		horizPoint hp;
    	hp.alt = asin(sin(latInRad) * sin(dec) + cos(latInRad) * cos(dec) * cos(lha));
    	hp.az = atan2(sin(lha) * cos(dec), cos(lha) * cos(dec) * sin(latInRad) - sin(dec) * cos(latInRad)) + pi;
    	return hp;
	}
	pointing toHealpixPointing(){
		double phi = fmod(ra-facetRAinRad+pi, 2*pi);
		double theta = pi/2-dec+facetDecInRad;
		if ((theta < 0) || theta > pi) phi = fmod(phi + pi, 2*pi);
		if (theta < 0) theta = abs(theta);
		if (theta > pi) theta = (pi - (theta-pi));
		return pointing(theta, phi); //pointing is defined as the colatitude (0 to pi = N to S) and the longitude (0 to 2pi)
	}
};

//Converts cartesian vectors to altitude and azimuth. I needed to delcate it outside the struct because horizPoint hadn't been declared yet.
equaPoint horizPoint::toEqua(double LST){
	equaPoint ep;
	ep.dec = asin(sin(alt)*sin(latInRad) - cos(alt)*cos(latInRad)*cos(az-pi));
	ep.ra = LST*2.0*pi/24 - atan2(sin(az-pi),(cos(az-pi)*sin(latInRad)+tan(alt)*cos(latInRad))); //The one I think is right
	return ep;
}

//Converts a pointing on the healpix sphere (which as been rotated so that the facet is centered on the equator) to an equatorial pointing
equaPoint convertToEquaPoint(pointing& hpPoint){
	double RA = fmod(hpPoint.phi + facetRAinRad - pi, 2*pi);
	double dec = pi/2 - hpPoint.theta + facetDecInRad;
	if ((dec < -pi/2) || dec > pi/2) RA = fmod(RA + pi, 2*pi);
	if (dec < -pi/2) dec = -pi/2 + (-pi/2 - dec);
	if (dec > pi/2) dec = pi/2 - (dec - pi/2);
	return equaPoint(RA,dec);
}

//Stucture for keeping track of point sources
struct pointSource{
	double flux, spectralIndex;
	equaPoint point;
	pointSource(){}
	pointSource(equaPoint pointIn, double fluxIn, double spectralIndexIn) : point(pointIn), flux(fluxIn), spectralIndex(spectralIndexIn) {}
	double fluxAt(double freq){ 
		return flux * pow(freq/pointSourceReferenceFreq, -spectralIndex);
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// LOADING AND INTIALIZATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Loads specifications from Specifications.txt or a specified filename;
void loadSpecs(){
	fstream infile(specsFilename.c_str(),fstream::in);
	string dummy;
	infile >> dummy >> dummy >> freq;
	infile >> dummy >> dummy >> polarization;
	infile >> dummy >> dummy >> baselinesFile;
	infile >> dummy >> dummy >> visibilitiesFolder;
	infile >> dummy >> dummy >> trueSkyFilename;
	infile >> dummy >> dummy >> finalMapFilename;
	infile >> dummy >> dummy >> PSFfilename;
	infile >> dummy >> dummy >> healpixPixelFilename;
	infile >> dummy >> dummy >> extendedHealpixPixelFilename;
	infile >> dummy >> dummy >> pixelCoordinatesFilename;
	infile >> dummy >> dummy >> extendedPixelCoordiantesFilename;
	infile >> dummy >> dummy >> noiseCovFilename;
	infile >> dummy >> dummy >> DmatrixFilename;
	infile >> dummy >> dummy >> dataProductFilePrefix;
	infile >> dummy >> dummy >> arrayLat;
	infile >> dummy >> dummy >> arrayLong;
	infile >> dummy >> dummy >> facetRA;
	infile >> dummy >> dummy >> facetDec;
	infile >> dummy >> dummy >> facetSize;
	infile >> dummy >> dummy >> NSIDE;
	infile >> dummy >> dummy >> PSFextensionBeyondFacetFactor;
	infile >> dummy >> dummy >> snapshotTime;
	infile >> dummy >> dummy >> integrationTime;
	infile >> dummy >> dummy >> channelBandwidth;
	infile >> dummy >> dummy >> gaussianPrimaryBeam;
	infile >> dummy >> dummy >> primaryBeamFWHM;
	infile >> dummy >> dummy >> xpolOrientationDegreesEastofNorth;
	infile >> dummy >> dummy >> maximumAllowedAngleFromPBCenterToFacetCenter;
	infile >> dummy >> dummy >> noiseStd;
	infile >> dummy >> dummy >> overwriteVisibilitiesWithGSM;
	infile >> dummy >> dummy >> GSMres;
	infile >> dummy >> dummy >> overwriteVisibilitiesWithPointSources;
	infile >> dummy >> dummy >> alsoComputePointSourcePSF;
	infile >> dummy >> dummy >> pointSourceFile;
	infile >> dummy >> dummy >> pointSourceReferenceFreq;
	infile >> dummy >> dummy >> pointSourceFluxUpperLimit;
	infile >> dummy >> dummy >> pointSourceFluxLowerLimit;
	infile >> dummy >> dummy >> pointSourceOutputCatalogFilename;
	infile >> dummy >> dummy >> pointSourcePSFFilename;
	infile >> dummy >> dummy >> throwOutFluxOutsideExtendedFacet;
	infile.close();

	//Other global variables to set
	//nFull = int(round(mappingFieldSizeFactor * mappingOverresolution * 2 / (sqrt(4*pi / (12.0*NSIDE*NSIDE))) / 2)*2);
	k = 2 * pi * freq * 1e6 / c; 
	facetDecInRad = facetDec*pi/180;
	facetRAinRad = facetRA*pi/180;
	latInRad = pi/180.0*arrayLat;
	temperatureConversionFactor = (c*c) / (2 * 1.3806488e-23 * pow(freq*1e6,2) * 1e26); //multiply by this number to convert Jy to K
	trueSkyFilename = dataProductFilePrefix + trueSkyFilename;
	healpixPixelFilename = dataProductFilePrefix + healpixPixelFilename;
	extendedHealpixPixelFilename = dataProductFilePrefix + extendedHealpixPixelFilename;
	pixelCoordinatesFilename = dataProductFilePrefix + pixelCoordinatesFilename;
	extendedPixelCoordiantesFilename = dataProductFilePrefix + extendedPixelCoordiantesFilename;
	finalMapFilename = dataProductFilePrefix + finalMapFilename;
	PSFfilename = dataProductFilePrefix + PSFfilename;
	noiseCovFilename = dataProductFilePrefix + noiseCovFilename;
	DmatrixFilename = dataProductFilePrefix + DmatrixFilename;
	pointSourceOutputCatalogFilename = dataProductFilePrefix + pointSourceOutputCatalogFilename;
	pointSourcePSFFilename = dataProductFilePrefix + pointSourcePSFFilename;
}

//Loads baselines (only one orientation, no complex conjugates, from the file under the global variable "baselinesFile", which is in "south east up" format)
vector<cartVec> loadBaselines(vector<int>& baselineRedundancy){
	cout << "Now loading all the baseline vectors." << endl << "[NOTE: THIS ASSUMES THAT THE DATA FORMAT IS SOUTH EAST UP, THOUGH THE INTERNAL FORMAT IS EAST NORTH UP (RIGHT HANDED SYSTEM)]" << endl;
	vector<cartVec> baselines;
	double south, east, up;
	int multiplicity;
	fstream infile(baselinesFile.c_str(),fstream::in);
	while(infile >> south >> east >> up >> multiplicity){
		cartVec thisBaseline;
		thisBaseline.y = -south; 
		thisBaseline.x = east;
		thisBaseline.z = up;
		baselines.push_back(thisBaseline);
		baselineRedundancy.push_back(multiplicity);
	}
	infile.close();
	nBaselines = baselines.size();
	return baselines;
}

//Loads the pre-computed principal axes of the plane of the array, which is, in general, slightly sloped relative to the the EW/NS plane, and a third axis perpendicular to the two and pointed mostly up
//The two axes must be perpendicular to one another. By convention, the first vector will be mostly East with no N/S component and the second will be mostly South
//The rotation can be paraterized by two cartesian rotations, one about the x (EW) axis and one about the y axis (SN)...this is completely general.
//Since rotations are not commutative, I'll first perform the rotation around the x axis and then the rotation about the y axis. This will preserve the property that x has no N/S component.
vector<cartVec> loadArrayPrincipalAxes(){
	cout << "Now loading array's principal axes...[NOTE: THESE ARE HARD CODED TO AN ARBITARY (SMALL) VALUE OR 0 FOR NOW.]" << endl;
	double thetaX = 0;//.02; //radians of rotation about the EW axis taking north toward up.
	double thetaY = 0;//.04; //radians of rotation about the NS axis taking east toward up.
	/*cartVec axisE(1, 0, 0); //Due east
	cartVec axisN(0, 1, 0); //Due south
	cartVec axisErotX(1, 0, 0);
	cartVec axisNrotX(0, cos(thetaX), sin(thetaX));*/
	cartVec axisErotXY(cos(thetaY), 0, sin(thetaY));
	cartVec axisNrotXY(-sin(thetaX)*sin(thetaY), cos(thetaX), sin(thetaX)*cos(thetaY)); 
	cartVec axisUrotXY = axisErotXY.cross(axisNrotXY);
	vector<cartVec> arrayPrincipalAxes;
	arrayPrincipalAxes.push_back(axisErotXY);		
	arrayPrincipalAxes.push_back(axisNrotXY);
	arrayPrincipalAxes.push_back(axisUrotXY);
	return arrayPrincipalAxes;
}

//This figures out the projection of the baselines into the best-fit plane of the array. This is needed for the 2D FFT.
//The projected baselines are written in terms of their components in terms of the array principal axes, not in terms of East, North, and Up.
vector<cartVec> calculateProjectedBaselines(vector<cartVec>& baselines, vector<cartVec>& arrayPrincipalAxes){
	vector<cartVec> projectedBaselines;
	for (int b = 0; b < nBaselines; b++){
		cartVec	projectedBaselineInENU = baselines[b] - arrayPrincipalAxes[2]*(baselines[b].dot(arrayPrincipalAxes[2]));
		cartVec projectedBaselineInXYZ(projectedBaselineInENU.dot(arrayPrincipalAxes[0]), projectedBaselineInENU.dot(arrayPrincipalAxes[1]), 0.0);
		projectedBaselines.push_back(projectedBaselineInXYZ);
	} 
	return	projectedBaselines;
}

//Loads visibilities, which have an LST and a real and imaginary component, from a file with the format I used for the omniscope
//Format is allVisibilites[baseline number same as baselines vector][different measurements].re/im
//It is assumed that all visibilities have the same set of LSTs, which is taken from the first baseline loaded
vector< vector<complex> > loadVisibilities(vector<cartVec>& baselines, vector<double>& LSTs){	
	vector< vector<complex> > allVisibilities;
	cout << "Now loading all visibilities for " << freq << " MHz and " << polarization << " polarization..." << endl;
	for (int n = 0; n < nBaselines; n++){
		cout << " " << floor(100.0 * n / nBaselines) << "% done. \r" << std::flush;

		//Format filename
		stringstream ss;
		ss << visibilitiesFolder  << -baselines[n].y << "_m_south_" << baselines[n].x << "_m_east_" << baselines[n].z << "_m_up_" << polarization << "_pol_" << freq << "_MHz.dat";
		//ss << visibilitiesFolder << "Visibilties_for_" << round(-baselines[n].y) << "_m_south_" << round(baselines[n].x) << "_m_east_" << round(baselines[n].z) << "_m_up_" << polarization << "_pol_" << freq << "_MHz.dat";
		string infilename = ss.str();
		fstream infile(infilename.c_str(),fstream::in);
		
		//Load all LST, re, and im into a vector of vectors of visibilities
		vector<complex> thisBaselinesVisibilities;
		double localSiderealTime, real, imaginary;
		while(infile >> localSiderealTime >> real >> imaginary){
			if (n == 0) LSTs.push_back(localSiderealTime);
			complex thisVisibility(real, imaginary);
			thisBaselinesVisibilities.push_back(thisVisibility);
		}
		infile.close();
		if (thisBaselinesVisibilities.size() == 0) cout << "WARNING: Cannot find " << infilename << endl;
		allVisibilities.push_back(thisBaselinesVisibilities);
	}
	nLSTs = LSTs.size();
	cout << "Done.                  " << endl;  
	return allVisibilities;
}

vector< vector<double> > loadNoiseVarianceOnEachVisibiltiy(vector<cartVec>& baselines, vector<double>& LSTs, vector<int> baselineRedundancy){
	//TODO: ENVENTUALLY I'LL WANT TO LOAD A MODEL FOR THE NOISE ON ANTENNA AND THEN COMPUTE WHAT THE PER VISIBILITY NOISE IS. FOR NOW I'LL JUST USE A SIMPLE MODEL.
	cout << "Now loading and computing the noise variance on each visibility...[NOTE: NOISE SET TO A SIMPLISTIC MODEL FOR NOW.]" << endl;
	vector< vector<double> > noiseVarianceOnEachVisibility(nBaselines, vector<double>(nLSTs,0));
	//cout << endl << "ALL NOISE VARIANCES SET TO 1" << endl << endl;
	for (int t = 0; t < nLSTs; t++){
		for (int b = 0; b < nBaselines; b++){
			noiseVarianceOnEachVisibility[b][t] = pow(noiseStd,2)/(integrationTime * channelBandwidth * 1e6 * baselineRedundancy[b]);
			//noiseVarianceOnEachVisibility[b][t] = 1;
		}
	}
	return noiseVarianceOnEachVisibility;
}

//This function loads the location of the primary beam in horizontal coordinates as a function of LST
vector<horizPoint> loadPBPointings(vector<double>& LSTs){
	//TODO: EVENTUALLY, WE WANT TO GENERALIZE THIS TO WHERE THE PRIMARY BEAM IS POINTED, BUT FOR NOW WE'LL ASSUME IT'S THE ZENITH
	cout << "Now loading all the primary beam pointings..." << endl;
	vector<horizPoint> PBpointings;
	for (int n = 0; n < nLSTs; n++){
		horizPoint thisPointing(pi/2, 0);
		PBpointings.push_back(thisPointing);
	}
	return PBpointings;
}

//This funciton loads the primary beam into a array of discretized values of alt and az.
//The data file, somewhat inconsistently, has azimuth 0 in the direction pointed by the XX polarization and increases CCW
vector< vector<double> > loadDiscretizedPrimaryBeam(){
	cout << "Now loading the primary beams for the nearest two frequencies and interpolating/extrapolating between them..." << endl;
	cout << "[NOTE: FOR NOW, PRIMARY BEAM IS ASSUMED TO BE THE SAME FOR ALL OBSERVATIONS]" << endl;
	vector< vector<double> > primaryBeamDiscretized(nAltBeamPoints, vector<double>(nAzBeamPoints, 0.0));
	double freq1 = -1;
	double freq2 = -1;
	for (double f = firstBeamFreq; f <= lastBeamFreq; f+=beamDeltaFreq){
		if (fabs(f - freq) < fabs(freq1 - freq)) freq1 = f;	
	} 
	for (double f = firstBeamFreq; f <= lastBeamFreq; f+=beamDeltaFreq){
		if ((fabs(f - freq) < fabs(freq2 - freq)) && (f != freq1)) freq2 = f;
	}
	stringstream freq1stream, freq2stream;
	freq1stream << beamDirectory << polarization << "_" << freq1 << ".dat";
	freq2stream << beamDirectory << polarization << "_" << freq2 << ".dat";
	string file1 = freq1stream.str();
	string file2 = freq2stream.str();
	fstream infile1(file1.c_str(),fstream::in);
	fstream infile2(file2.c_str(),fstream::in);
	double gain1, gain2;
	for (int alt = 0; alt < nAltBeamPoints; alt++){
		for (int az = 0; az < nAzBeamPoints; az++){
			infile1 >> gain1;
			infile2 >> gain2;
			primaryBeamDiscretized[alt][az] = (gain2 - gain1)/(freq2 - freq1) * (freq - freq1) + gain1;
		}
	}
	return primaryBeamDiscretized;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OUTPUT FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void exportMatrix(vector< vector<double> >& matrix, int dim1, int dim2, string filename){
	cout << "Now saving " << filename << " as a " << dim1 << " by " << dim2 << " matrix..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filename.c_str(), ios::trunc);		
	for (int i = 0; i < dim1; i++){
		for (int j = 0; j < dim2; j++){
			outfile << matrix[i][j] << " ";
		}
		outfile << endl;
	} 
	outfile.close();
}

void exportVector(vector<double>& vec, int dim, string filename){
	cout << "Now saving " << filename << " as a size " << dim << " vector..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filename.c_str(), ios::trunc);		
	for (int i = 0; i < dim; i++) outfile << vec[i] << endl;
	outfile.close();
}

void exportPixels(vector<int>& vec, string filename){
	cout << "Now saving " << filename << " as a list of " << vec.size() << " healpix pixels..." << endl;
	ofstream outfile;
	outfile.open(filename.c_str(), ios::trunc);		
	for (int i = 0; i < vec.size(); i++) outfile << vec[i] << endl;
	outfile.close();
}

void exportPartialHealpixMap(Healpix_Map<double>& map, vector<int>&  healpixIndices, string filename){
	cout << "Now saving " << filename << " as a size " << healpixIndices.size() << " portion of a HEALPIX map..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filename.c_str(), ios::trunc);		
	for (int i = 0; i < healpixIndices.size(); i++) outfile << map[healpixIndices[i]] << endl;
	outfile.close();
}

void exportCoordinates(vector<equaPoint>& pixelCoordinates, string filename){
	cout << "Now saving " << filename << " as a list of " << pixelCoordinates.size() << " healpix pixel equatorial coordinates..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filename.c_str(), ios::trunc);		
	for (int i = 0; i < pixelCoordinates.size(); i++) outfile << pixelCoordinates[i].ra << "     " << pixelCoordinates[i].dec << endl;
	outfile.close();
}

void exportComplexMatrix(vector< vector<complex> >& matrix, int dim1, int dim2, string filenameReal, string filenameImag){
	cout << "Now saving " << filenameReal << " as the " << dim1 << " by " << dim2 << " as the real matrix..." << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(filenameReal.c_str(), ios::trunc);		
	for (int j = 0; j < dim2; j++){
		for (int i = 0; i < dim1; i++){
			outfile << matrix[i][j].re << " ";
		}
		outfile << endl;
	} 
	outfile.close();

	cout << "Now saving " << filenameImag << " as the " << dim1 << " by " << dim2 << " as the imaginary matrix..." << endl;
	ofstream outfile2;
	outfile2.precision(16);	
	outfile2.open(filenameImag.c_str(), ios::trunc);		
	for (int j = 0; j < dim2; j++){
		for (int i = 0; i < dim1; i++){
			outfile2 << matrix[i][j].im << " ";
		}
		outfile2 << endl;
	} 
	outfile2.close();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GEOMETRY FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function computes the altitude and azimuth of the facet ceneter for all LSTs observed.
vector<horizPoint> computeFacetCenterPointings(vector<double>& LSTs){
	cout << "Now computing the horizontal angle to the facet center for all LSTs..." << endl;
	vector<horizPoint> facetCenterPointings;
	equaPoint facetCenter(facetRAinRad, facetDecInRad);
	for (int n = 0; n < nLSTs; n++) facetCenterPointings.push_back(facetCenter.toHoriz(LSTs[n]));
	return facetCenterPointings;
}

//This function determines whether the distance between the facet center and the center of the maximumAllowedAngleFromPBCenterToFacetCenter
vector<bool> determineIntegrationsToUse(vector<horizPoint>& PBpointings, vector<horizPoint>& facetCenterPointings){
	vector<bool> toUse;
	for (int n = 0; n < PBpointings.size(); n++){ 
		toUse.push_back(facetCenterPointings[n].greatCircleAngle(PBpointings[n]) < pi/180*maximumAllowedAngleFromPBCenterToFacetCenter);
		if (toUse[n]) nIntegrationsToUse++;
	}
	cout << nIntegrationsToUse << " of the " << PBpointings.size() << " total integrations have the facet center within " << maximumAllowedAngleFromPBCenterToFacetCenter << " degrees of the primary beam center." << endl;
	return toUse;
}

//This function assigns LST indices to snapshots, which are gridded together. All integrations in a snapshot must be perfectly sequential and each snapshot must be full (remaining integrations are discarded).
vector< vector<int> > assignLSTsToSnapshots(vector<bool> integrationsToUse){
	vector< vector<int> > snapshotLSTindices;	
	for (int t = 0; t < nLSTs; t++){	
		if (integrationsToUse[t]){ //if ready to start a snapshot
			bool sequential = true;
			for (int i = 0; i < snapshotTime/integrationTime; i++){
				if (!integrationsToUse[t+i]) sequential = false;
			}
			if (sequential){
				vector<int> integrationIndices;
				for (int i = 0; i < snapshotTime/integrationTime; i++){
					integrationIndices.push_back(t);
					t++;
				}
				snapshotLSTindices.push_back(integrationIndices);
				t--;
			}
		}
	}
	nSnapshots = snapshotLSTindices.size();
	cout << "These integrations have been grouped into " << nSnapshots << " snapshots of exactly " << snapshotTime << " seconds. " << nIntegrationsToUse - nSnapshots*snapshotTime/integrationTime << " integration(s) were discarded." << endl;
	return snapshotLSTindices;
}

//This function creates an empty healpix_map of doubles 
Healpix_Map<double> emptyHealpixMap(){
	Healpix_Map<double> emptyMap(int(round(log(NSIDE)/log(2))), RING);
	for (int n = 0; n < (12*NSIDE*NSIDE); n++) emptyMap[n] = 0;
	return emptyMap;
}

//This function figures out which HEALPIX indices are part of map and which are part of the extended map (used for computing the PSF)
vector<int> computeHealpixIndices(bool extended){
	Healpix_Map<double> sampleHealpixMap = emptyHealpixMap();
	equaPoint facetCenterEquaPoint(facetRAinRad, facetDecInRad);
	vector<int> healpixIndices;
	double angularRadius = facetSize/360*2*pi/2.0;
	if (extended) angularRadius *= PSFextensionBeyondFacetFactor;
	sampleHealpixMap.query_disc(facetCenterEquaPoint.toHealpixPointing(), angularRadius, healpixIndices);
	if (extended){
		nPixelsExtended = healpixIndices.size();
	} else {
		nPixels = healpixIndices.size();
	}
	pointing test = facetCenterEquaPoint.toHealpixPointing();

	return healpixIndices;
}

//This function creates a vector of pointings in equatorial coordinates for each pixel in the HEALPIX map
vector<equaPoint> computeEquaPointings(vector<int>& indices){
	Healpix_Map<double> sampleHealpixMap = emptyHealpixMap();
	vector<equaPoint> pixelEquaPointings;
	for (int n = 0; n < indices.size(); n++){
		pointing healpixPointing = sampleHealpixMap.pix2ang(indices[n]);
		equaPoint thisEquaPoint = convertToEquaPoint(healpixPointing);
		pixelEquaPointings.push_back(thisEquaPoint);
	}
	return pixelEquaPointings;
}

//Creates a healpix map of bools for whether or not a given healpix pixel is in the part of the sky we're trying to map.
Healpix_Map<bool> isInFacet(vector<int>& healpixIndices){
	int order = int(round(log(NSIDE)/log(2)));
	Healpix_Map<bool> mapOfPixelsInFacet(order, RING);
	for (int n = 0; n < NSIDE*NSIDE*12; n++) mapOfPixelsInFacet[n] = false;
	for (int i = 0; i < healpixIndices.size(); i++) mapOfPixelsInFacet[healpixIndices[i]] = true;
	return mapOfPixelsInFacet;
}

//Creates a vector that maps where map[nth index in healpixIndices] = mth index in extendedHealpixIndices
vector<int> mapIndicesFromFacetToExtendedFacet(vector<int>& healpixIndices, vector<int>& extendedHealpixIndices){
	vector<int> mapOfIndicesInExtendedIndexVector(nPixels, 0);
	int extendedCounter = 0;
	for (int n = 0; n < nPixels; n++){
		while (extendedHealpixIndices[extendedCounter] != healpixIndices[n]){
			extendedCounter++;
			if (extendedCounter >= nPixelsExtended){
				cout << "UNABLE TO MATCH EXTENDED FACET PIXEL INDICES TO FACET PIXEL INDICES!!!" << endl;
			}
		}
		mapOfIndicesInExtendedIndexVector[n] = extendedCounter;
	}
	return mapOfIndicesInExtendedIndexVector;
}

//Converts galactic coordinates to equatorial
equaPoint convertGalToEquaPoint(pointing& galPointing){
	double b = pi/2.0 - galPointing.theta;
	double l = galPointing.phi;
	double pole_ra = 2.0*pi/360.0*192.859508;
    double pole_dec = 2.0*pi/360.0*27.128336;
    double posangle = 2.0*pi/360.0*(122.932-90.0); //32.932
    double raHere = atan2((cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle))) + pole_ra;
    double decHere = asin(cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec));
    return equaPoint(raHere, decHere);
}

//Converts equatorial coordinates to galactic
pointing convertEquaPointToGal(equaPoint& eq){
	double pole_ra = 2.0*pi/360.0*192.859508;
    double pole_dec = 2.0*pi/360.0*27.128336;
    double posangle = 2.0*pi/360.0*(122.932-90.0);
	double b = asin( sin(eq.dec)*sin(pole_dec) + cos(eq.dec)*cos(pole_dec)*cos(pole_ra - eq.ra) );
	double l = (posangle+pi/2) - atan2(cos(eq.dec)*sin(eq.ra-pole_ra), cos(pole_dec)*sin(eq.dec) - sin(pole_dec)*cos(eq.dec)*cos(eq.ra-pole_ra));	
	return pointing(pi/2-b, l);
}

//This function retruns the value of the gain of the PB as a funciton of altitude and azimuth
//Unfortunately, the primary beam azimuth is stored with the XX polarization as azimuth zero and continues CCW. 
double primaryBeam(horizPoint& pointing, vector< vector<double> >& PB){
	//return 1;
	if (gaussianPrimaryBeam){
		if (pointing.alt > 0){
			double sigma = primaryBeamFWHM/360.0*2*pi/2.355;
			return exp(-pow(pi/2 -pointing.alt,2)/2/pow(sigma,2));///sigma/sqrt(2*pi); //Normalized so that PB(zenith)=1
		} else {
			return 0.0;
		}
	}
	if (pointing.alt <= 0) return 0.0;

	double altPixel = pointing.alt * 180.0 / pi / deltaBeamAlt;
	double azPixel = fmod(-pointing.az * 180.0 / pi + xpolOrientationDegreesEastofNorth + 360.0,360.0) / deltaBeamAz;
	int altIndex1 = int(floor(altPixel));
	int altIndex2 = int(ceil(altPixel));
	int azIndex1 = int(floor(azPixel));
	int azIndex2 = int(ceil(azPixel));
	
	//Handle integer pixel values to avoid getting 0/0
	if ((altIndex1 == altIndex2) && (azIndex1 == azIndex2)) return (PB[altIndex1][azIndex1]);
	if (altIndex1 == altIndex2) return ((PB[altIndex1][azIndex2] - PB[altIndex1][azIndex1]) * (azPixel - azIndex1) + PB[altIndex1][azIndex1]);
	if (azIndex1 == azIndex2) return ((PB[altIndex2][azIndex1] - PB[altIndex1][azIndex1]) * (altPixel - altIndex1) + PB[altIndex1][azIndex1]);

	double PBresult = (PB[altIndex1][azIndex1] * (altIndex2-altPixel) * (azIndex2-azPixel));
	PBresult += (PB[altIndex2][azIndex1] * (altPixel-altIndex1) * (azIndex2-azPixel));
	PBresult += (PB[altIndex1][azIndex2] * (altIndex2-altPixel) * (azPixel-azIndex1));
	PBresult += (PB[altIndex2][azIndex2] * (altPixel-altIndex1) * (azPixel-azIndex1));
	return PBresult;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS FOR OVERWRITING THE VISIBILITIES AND DEALING WITH POINT SOURCES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Load the components data and properly reweight and rescale the GSM principal componenets
Healpix_Map<double> computeGSM(){
	int principalComps = 3;
	vector< Healpix_Map<double> > GSMpcomps;
	for (int p = 1; p <= principalComps; p++){
		stringstream ss;
		ss << FITS_directory << GSMformat1 << p << GSMformat2 << GSMres;
		string GSMfile = ss.str();
		cout << "Now loading " << GSMfile << endl;
		Healpix_Map<double> map; 
		read_Healpix_map_from_fits(GSMfile,map); 
		GSMpcomps.push_back(map);
	}

	//Load in the componenets data
	stringstream ss;
	string componentsFile = "components.dat";
	ss << FITS_directory << componentsFile;
	fstream infile((ss.str()).c_str(),fstream::in);
	vector<double> pcompFreqs(componentsToFit,0);
	vector<double> pcompTemps(componentsToFit,0);
	vector< vector<double> > weights(3, vector<double>(componentsToFit,0));
	for (int n = 0; n < componentsToFit; n++){
		double f,T,p1,p2,p3;
		infile >> f >> T >> p1 >> p2 >> p3;
		pcompFreqs[n] = f;
		pcompTemps[n] = T;
		weights[0][n] = p1;
		weights[1][n] = p2;
		weights[2][n] = p3;
	}
	infile.close();
	
	//Determine the temperature though a power law fit
	double* logFreqs = new double[componentsToFit];
	double* logTemps = new double[componentsToFit];
	for (int n = 0; n < componentsToFit; n++){
		logFreqs[n] = log10(pcompFreqs[n]);
		logTemps[n] = log10(pcompTemps[n]);
	}
	double c0, c1, cov00, cov01, cov11, sumsq;
	gsl_fit_linear(logFreqs, 1, logTemps, 1, componentsToFit, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	delete[] logFreqs;
	delete[] logTemps;

	double overallTemperature = pow(10.0,c0)*pow(freq,c1);
	cout << "A f^" << c1 << " power law yields a GSM temperature normalization: " << overallTemperature << " K." << endl;

	//Determine the weights through linear interpolation
	int freqCounter = 1;
	while (freqCounter < componentsToFit){
		if (pcompFreqs[freqCounter] > freq) break;
		freqCounter++;
	}
	vector<double> w(3,0);
	for (int i = 0; i < 3; i++){
		double slope = (weights[i][freqCounter] - weights[i][freqCounter-1]) / (pcompFreqs[freqCounter] - pcompFreqs[freqCounter-1]);
		w[i] = weights[i][freqCounter-1] + slope*(freq - pcompFreqs[freqCounter-1]);
	} 

	//Combine the components and return the final map;
	Healpix_Map<double> GSM = GSMpcomps[0];
	double sum = 0;
	for (int n = 0; n < 12*GSMres*GSMres; n++){
		GSM[n] = (GSMpcomps[0][n]*w[0] + GSMpcomps[1][n]*w[1] + GSMpcomps[2][n]*w[2]) * overallTemperature;	
		sum += GSM[n];
	} 
	cout << "Average temperature before rebinning is " << sum / 12/GSMres/GSMres << " K." << endl;
	return GSM;
}

//Rebins the GSM into equatorial coordinates, making for an easier comparison between the true sky
Healpix_Map<double> computeGSMInEquatorialCoordinates(vector<int>& extendedHealpixIndices){
	int GSMresTemp = GSMres;
	if (GSMres > 512) GSMres = 512;
	Healpix_Map<double> galacticGSM = computeGSM();
	GSMres = GSMresTemp;

	Healpix_Map<double> equatorialGSM(int(round(log(GSMres)/log(2))), RING);
	for (int n = 0; n < 12*GSMres*GSMres; n++){
		pointing thisPoint = equatorialGSM.pix2ang(n);
		equaPoint thisEquaPoint = convertToEquaPoint(thisPoint);
		//equaPoint thisEquaPoint(thisPoint.phi, pi/2-thisPoint.theta);
		pointing thisGalPoint = convertEquaPointToGal(thisEquaPoint);
		equatorialGSM[n] = galacticGSM.interpolated_value(thisGalPoint);
	}

	double sum = 0;
	for (int n = 0; n < 12*GSMres*GSMres; n++) sum += equatorialGSM[n];
	cout << "Average temperature after rebinning is " << sum / 12 / GSMres / GSMres << " K." << endl;

	Healpix_Map<double> galacticGSMlowres = emptyHealpixMap();
	if (GSMres > NSIDE){
		galacticGSMlowres.Import_degrade(equatorialGSM);	
	} else {
		galacticGSMlowres = equatorialGSM;
	}
	
	/*vector<double> wholeGSM(NSIDE*NSIDE*12,0.0);
	for (int n = 0; n < NSIDE*NSIDE*12; n++) wholeGSM[n] = equatorialGSM[n];
	exportVector(wholeGSM, NSIDE*NSIDE*12, "wholeGSM.dat");	*/


	vector<double> facetMap(nPixelsExtended, 0.0);
	for (int n = 0; n < nPixelsExtended; n++) facetMap[n] = galacticGSMlowres[extendedHealpixIndices[n]];
	exportVector(facetMap, nPixelsExtended, trueSkyFilename);
	string onlyGSMFilename = dataProductFilePrefix + "trueSkyOnlyGSM.dat";
	exportVector(facetMap, nPixelsExtended, onlyGSMFilename);
	return equatorialGSM;
}


//Overwrite the visibilities with the true sky GSM, eliminating a possible source of error in the calculation. Eventually this will be unnecessary.
void overwriteVisibilitiesWithTheGSM(vector< vector<complex> >& allVisibilities, vector<cartVec>& baselines, vector< vector<double> >& PB, vector<double>& LSTs, vector<horizPoint>& facetCenterPointings,
	vector<int>& extendedHealpixIndices, vector<equaPoint>& extendedPixelEquaPointings, vector<bool>& integrationsToUse, Healpix_Map<bool>& mapOfPixelsInFacet){
	
	Healpix_Map<double> equatorialGSM = computeGSMInEquatorialCoordinates(extendedHealpixIndices);
	cout << endl << "Overwriting all visibilities with GSM-" << GSMres << "!!!" << endl << endl;

	for (int l = 0; l < LSTs.size(); l++){
		if (integrationsToUse[l]){
			cout <<  "Now working on LST = " << LSTs[l] << endl;// "\r                 " << std::flush;
			for (int b = 0 ; b < nBaselines; b++){
				allVisibilities[b][l].re = 0;
				allVisibilities[b][l].im = 0;
			}
			vector<cartVec> allPixelCartVecs;
			vector<horizPoint> allPixelHorizPoint;
			vector<equaPoint> pixelEquaPointings;
			for (int n = 0; n < 12*GSMres*GSMres; n++){
				pointing healpixPointing = equatorialGSM.pix2ang(n);
				equaPoint thisEquaPoint = convertToEquaPoint(healpixPointing);
				horizPoint thisHorizPoint = thisEquaPoint.toHoriz(LSTs[l]);
				cartVec thisCartVec = thisHorizPoint.toVec();
				allPixelCartVecs.push_back(thisCartVec);
				allPixelHorizPoint.push_back(thisHorizPoint);
			}
			for (int n = 0; n < 12*GSMres*GSMres; n++){
				if (allPixelHorizPoint[n].alt > 0){
					if (!throwOutFluxOutsideExtendedFacet || mapOfPixelsInFacet[n]){
						double constantAtThisPoint = equatorialGSM[n] * primaryBeam(allPixelHorizPoint[n], PB) * 4 * pi / (12 * GSMres * GSMres);
						for (int b = 0; b < nBaselines; b++){
							double b_dot_k = baselines[b].dot(allPixelCartVecs[n]) * (-k);
							allVisibilities[b][l].re += constantAtThisPoint * cos(b_dot_k);
							allVisibilities[b][l].im += constantAtThisPoint * sin(b_dot_k);
						}
					}
				}
			}
		}
	}
}

//Loads point sources into a vector, only keeping those above a certain flux * primary beam at central LST limit
void loadPointSources(vector<pointSource>& allPointSources, vector< vector<double> >& PB){	
	cout << "Now loading point sources from " << pointSourceFile << "." << endl;
	horizPoint zenith(pi/2, 0);
	double centralLST = facetRAinRad / (2*pi) * 24;
	fstream infile(pointSourceFile.c_str(),fstream::in);
	double ra, dec, flux, spectralIndex;
	while(infile >> ra >> dec >> flux >> spectralIndex){
		if (flux > pointSourceFluxLowerLimit){
			equaPoint psEquaPoint(ra*2*pi/360, dec*2*pi/360);
			horizPoint psHorizPoint = psEquaPoint.toHoriz(centralLST);
			double fluxTimesPB = flux * primaryBeam(psHorizPoint, PB);
			if (fluxTimesPB > pointSourceFluxLowerLimit && fluxTimesPB < pointSourceFluxUpperLimit){
				if (!throwOutFluxOutsideExtendedFacet || zenith.greatCircleAngle(psHorizPoint) < .8*PSFextensionBeyondFacetFactor*facetSize*2*pi/360/2){
					pointSource thisPointSource(psEquaPoint, flux, spectralIndex);
					allPointSources.push_back(thisPointSource);
				}
			}
		}
	}
	infile.close();
	nPointSources = allPointSources.size();
	cout << "Found " << nPointSources << " point sources between the beam times flux limits of " << pointSourceFluxLowerLimit << " Jy and " << pointSourceFluxUpperLimit << " Jy." << endl;

	cout << "Now saving all " << nPointSources << " point sources in " << pointSourceOutputCatalogFilename << endl;
	ofstream outfile;
	outfile.precision(16);	
	outfile.open(pointSourceOutputCatalogFilename.c_str(), ios::trunc);		
	for (int i = 0; i < nPointSources; i++){
		outfile << allPointSources[i].point.ra << " " << allPointSources[i].point.dec << " " << allPointSources[i].flux << " " << allPointSources[i].spectralIndex << " " << allPointSources[i].fluxAt(freq) << endl;
	} 
	outfile.close();
}

//This function takes a point source and saves the relevant information for a matlab code to reconstruct an optimal "true sky"
void savePointSourceParametersForDiscreteMap(vector<pointSource>& allPointSources, vector<cartVec>& baselines){
	
	string outfilename = "point_source_true_sky_parameters.dat";
	outfilename = dataProductFilePrefix + outfilename;
	ofstream outfile;
	outfile.open(outfilename.c_str(), ios::trunc);
	cout << "Now calculating and saving " << outfilename << "..." << endl;
	Healpix_Map<double> trueSky = emptyHealpixMap();
	for (int p = 0; p < nPointSources; p++){
		//find 4 nearest neighbors and figure out their k vectors
		fix_arr<int, 8> neighbors;
		pointing hpPoint = (allPointSources[p].point).toHealpixPointing();
		trueSky.neighbors(trueSky.ang2pix(hpPoint),neighbors);
		double flux = allPointSources[p].fluxAt(freq);
		vector<int> nearestPixels;
		nearestPixels.push_back(trueSky.ang2pix(hpPoint));
		
		double centralLST = facetRAinRad / (2*pi) * 24;
		vector<double> distances(8,0.0);
		horizPoint pointSourceHoriz = (allPointSources[p].point).toHoriz(centralLST);
		for (int n = 0; n < 8; n++){
			pointing pixelPointingHere = trueSky.pix2ang(neighbors[n]);
			equaPoint thisPixEqua = convertToEquaPoint(pixelPointingHere);
			horizPoint thisPixHoriz = thisPixEqua.toHoriz(centralLST);
			distances[n] = pointSourceHoriz.greatCircleAngle(thisPixHoriz);
		}
		for (int i = 1; i < 4; i++){
			double minDist = 1e10;
			int minDistN = -1;
			for (int n = 0; n < 8; n++){
				if (distances[n] < minDist){
					minDist = distances[n];
					minDistN = n;
				} 
			}
			nearestPixels.push_back(neighbors[minDistN]);
			distances[minDistN] = 1e10;
		}
		for (int i = 0; i < 4; i++) outfile << nearestPixels[i] << " ";
			
		//pick the first three baselines (TODO: these are assumed to be linearly independent)
		cartVec pointSourceVec = pointSourceHoriz.toVec();
		vector<cartVec> deltaKVectors;
		for (int i = 0; i < 4; i++){
			pointing pixelPointingHere = trueSky.pix2ang(nearestPixels[i]);
			equaPoint thisPixEqua = convertToEquaPoint(pixelPointingHere);
			horizPoint thisPixHoriz = thisPixEqua.toHoriz(centralLST);
			deltaKVectors.push_back(thisPixHoriz.toVec() - pointSourceVec);
		}
		for (int b = 0; b < 3; b++){
			for (int d = 0; d < 4; d++){
				outfile << deltaKVectors[d].dot(baselines[b].normalize()) << " ";
			}
		}
		outfile << endl;
	}
	outfile.close();

}


//Overwrite the visibilities with a set of bright point sources, loaded from a catalog
void overwriteVisibilitiesWithAListOfPointSources(vector< vector<complex> >& allVisibilities, vector<cartVec>& baselines, vector< vector<double> >& PB, vector<double>& LSTs, 
	vector<horizPoint>& facetCenterPointings, vector<pointSource>& allPointSources, vector<bool> integrationsToUse, vector<int>& extendedHealpixIndices){
	
	cout << endl << "Overwriting all visibilities with point sources from " << pointSourceOutputCatalogFilename << "!!!" << endl << endl;	
	for (int l = 0; l < LSTs.size(); l++){
		if (integrationsToUse[l]){
			cout <<  "Now working on LST = " << LSTs[l] << endl;// "\r                 " << std::flush;
			if (!overwriteVisibilitiesWithGSM){
				for (int b = 0 ; b < nBaselines; b++){
					allVisibilities[b][l].re = 0;
					allVisibilities[b][l].im = 0;
				}
			}
			for (int i = 0; i < nPointSources; i++){
				horizPoint thisHorizPoint = (allPointSources[i].point).toHoriz(LSTs[l]);
				cartVec thisCartVec = thisHorizPoint.toVec();
				double constantAtThisPoint = allPointSources[i].fluxAt(freq) * primaryBeam(thisHorizPoint, PB);
				for (int b = 0; b < nBaselines; b++){
					double b_dot_k = baselines[b].dot(thisCartVec) * (-k);
					allVisibilities[b][l].re += constantAtThisPoint * cos(b_dot_k);
					allVisibilities[b][l].im += constantAtThisPoint * sin(b_dot_k);
				}
			}
		}
	}

	savePointSourceParametersForDiscreteMap(allPointSources, baselines);

	//Calculate the true sky in the simplest way and save it
	Healpix_Map<double> trueSky = emptyHealpixMap();
	fix_arr<int, 4> pix;
	fix_arr<double, 4> wgt;
	for (int p = 0; p < nPointSources; p++){
		pointing hpPoint = allPointSources[p].point.toHealpixPointing();
		trueSky.get_interpol(hpPoint, pix, wgt);
		double flux = allPointSources[p].fluxAt(freq);
		for (int i = 0; i < 4; i++) trueSky[pix[i]] += wgt[i] * flux * (12 * NSIDE * NSIDE / 4 / pi) ;
	}
	vector<double> facetMap(nPixelsExtended, 0.0);
	for (int n = 0; n < nPixelsExtended; n++) facetMap[n] = trueSky[extendedHealpixIndices[n]];
	if (overwriteVisibilitiesWithGSM){
		fstream infile(trueSkyFilename.c_str(),fstream::in);
		double GSMpixel = 0;
		for (int n = 0; n < nPixelsExtended; n++){
			infile >> GSMpixel;
			facetMap[n] = facetMap[n] + GSMpixel;
		}
		infile.close();
	} 
	exportVector(facetMap, nPixelsExtended, trueSkyFilename);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// REPHASE AND RENORMALIZE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function calculates N^-1 y, the first step in the mapmaking algorithm
void noiseInverseVarianceWeight(vector< vector<complex> >& allVisibilities, vector< vector<double> >& noiseVarianceOnEachVisibiltiy, vector<bool>& integrationsToUse){
	cout << "Now multiplying all visibilities by the inverse noise variance...[NOTE: FOR NOW, EACH VISIBILITY HAS A NOISE VARIANCE OF 1]" << endl;
	for (int t = 0; t < nLSTs; t++){
		if (integrationsToUse[t]){
			for (int b = 0; b < nBaselines; b++) allVisibilities[b][t] = allVisibilities[b][t]*(1.0/noiseVarianceOnEachVisibiltiy[b][t]);
		}
	}
}

//This function multiplies each visibility by e^(i*b.k_0) where b is the baseline vector and k_0 points to the facet center
/*void rephaseVisibilities(vector< vector<complex> >& allVisibilities, vector<bool>& integrationsToUse, vector<cartVec>& baselines, vector<horizPoint>& facetCenterPointings){
	cout << "Now rephasing all visibilities to the facet center..." << endl;
	for (int t = 0; t < nLSTs; t++){
		if (integrationsToUse[t]){
			for (int b = 0; b < nBaselines; b++){
				double b_dot_k = baselines[b].dot(facetCenterPointings[t].toVec())*(k); 
				complex complexFactor = complex(cos(b_dot_k), sin(b_dot_k)); //In A we multiply by e^b.k so in A^t we multiply by e^-b.k
				allVisibilities[b][t] = allVisibilities[b][t]*(complexFactor);
			}
		}
	}
}*/

void rephaseVisibilities(vector< vector<complex> >& allVisibilities, vector< vector<int> >& snapshotLSTindices, vector<double> LSTs, vector<cartVec>& baselines, vector<horizPoint>& facetCenterPointings){
	for (int s = 0; s < snapshotLSTindices.size(); s++){
		int snapshotCentralLSTindex = snapshotLSTindices[s][int(round(snapshotLSTindices[s].size()/2.0-.5))];
		cartVec centralSnapshotFacetCenterVector = facetCenterPointings[snapshotCentralLSTindex].toVec();
		for (int i = 0; i < snapshotLSTindices[s].size(); i++){
			cartVec deltaTheta = centralSnapshotFacetCenterVector - (facetCenterPointings[snapshotLSTindices[s][i]]).toVec();
			for (int b = 0; b < nBaselines; b++){
				double argument = -k * baselines[b].dot(deltaTheta);
				allVisibilities[b][snapshotLSTindices[s][i]] = allVisibilities[b][snapshotLSTindices[s][i]]	* complex(cos(argument),sin(argument));
			}
		}
	}
}

void convertToTemperatureUnits(vector< vector<complex> >& allVisibilities, vector<bool>& integrationsToUse){
	cout << "Now converting the visibilities to temperature units..." << endl;
	for (int t = 0; t < nLSTs; t++){
		if (integrationsToUse[t]){
			for (int b = 0; b < nBaselines; b++) allVisibilities[b][t] = allVisibilities[b][t] * temperatureConversionFactor;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAPMAKING AND PSF FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//This function adds up all the visibilities corresponding to all baselines during a snapshot
//The negative baselines' visibiltiies are stores in indicies nBaselines through 2*nBaselines - 1
//TODO: save a factor of 2 in memoery and 4 in runtime by absorpbing the negative baselines into a simplified set of calculations
vector<complex> calculateNinvTimesy(vector< vector<complex> >& allVisibilities, vector<int>& LSTindices){
	vector<complex> NinvTimesy(nBaselines, complex(0,0));
	for (int n = 0; n < LSTindices.size(); n++){
		for (int b = 0; b < nBaselines; b++){
			NinvTimesy[b] = NinvTimesy[b] + allVisibilities[b][LSTindices[n]];
		}
	}
	return NinvTimesy;
}

//This function calculates the diagonal matrix Ninv used in the calculation of Atranspose * N^-1 * A for the PSF
vector<double> calculateNinv(vector< vector<double> >& noiseVarianceOnEachVisibiltiy, vector<int>& LSTindices){
	vector<double> Ninv(nBaselines, 0.0);
	for (int n = 0; n < LSTindices.size(); n++){
		for (int b = 0; b < nBaselines; b++){
			Ninv[b] = Ninv[b] + 1.0/noiseVarianceOnEachVisibiltiy[b][LSTindices[n]];
		}
	}
	return Ninv;
}

//This function computes A^t matrix for each snapshot. The action of this matrix is to convert unique baselines into a dirty map.
vector< vector<complex> > calculateKAtranspose(double LST, vector<equaPoint>& extendedPixelEquaPointings, vector<cartVec>& baselines, vector< vector<double> >& PB){

	//Calculate the real space diagonal part of Atranspose 
	double constantFactor = 4*pi / (12.0*NSIDE*NSIDE) * temperatureConversionFactor; 
	vector<double> realSpaceDiagonalPart(nPixelsExtended, constantFactor);
	
	vector<horizPoint> extendedPixelHorizPointsAtThisLST;
	vector<cartVec> extendedPixelCartVecsAtThisLST;
	for (int n = 0; n < nPixelsExtended; n++){
		extendedPixelHorizPointsAtThisLST.push_back(extendedPixelEquaPointings[n].toHoriz(LST));
		extendedPixelCartVecsAtThisLST.push_back(extendedPixelHorizPointsAtThisLST[n].toVec());
		realSpaceDiagonalPart[n] *= (primaryBeam(extendedPixelHorizPointsAtThisLST[n],PB));
	}
	
	//Calculate a the relevant Fourier transform matrix (A has e^ib.k, so A^t has e^-ib.k = e^i|k|b.theta_hat, as expressed mathematically here.)
	vector< vector<complex> > Atranspose(nPixelsExtended, vector<complex>(nBaselines, complex(0,0)));
	for (int n = 0; n < nPixelsExtended; n++){
		for (int b = 0; b < nBaselines; b++){
			double argument = k * extendedPixelCartVecsAtThisLST[n].dot(baselines[b]);
			Atranspose[n][b].re = cos(argument) * realSpaceDiagonalPart[n];
			Atranspose[n][b].im = sin(argument) * realSpaceDiagonalPart[n];
		}
	}
	return Atranspose;
}

//This computes A^t * Ninv * y for the snapshot and adds it to the overall map
void addSnapshotMap(Healpix_Map<double>& coaddedMap, vector<complex>& NinvTimesy, vector< vector<complex> >& KAtranspose, vector<int>& mapOfIndicesInExtendedIndexVector, vector<int>& healpixIndices){
	for (int n = 0; n < nPixels; n++){
		for (int b = 0; b < nBaselines; b++){
			coaddedMap[healpixIndices[n]] = coaddedMap[healpixIndices[n]] + 2*(KAtranspose[mapOfIndicesInExtendedIndexVector[n]][b] * NinvTimesy[b]).re;
		}
	}
}

//This function computes KAtNinvAKt for each snapshot and then adds them to the existing sum to get an overall PSF
void addSnapshotPSF(vector< vector<double> >& PSF, vector< vector<complex> >& KAtranspose, vector<double>& Ninv, vector<int>& mapOfIndicesInExtendedIndexVector){
	for (int n = 0; n < nPixels; n++){
		for (int nx = 0; nx < nPixelsExtended; nx++){
			for (int b = 0; b < nBaselines; b++){
				PSF[n][nx] = PSF[n][nx] + 2 *(KAtranspose[nx][b] * KAtranspose[mapOfIndicesInExtendedIndexVector[n]][b].conj() * Ninv[b]).re;
			}
		}
	}
}

//Calculates the matrix that maps baseline visibilities to the exact positions of the point sources.
vector< vector<complex> > calcualtePointSourceAmatrix(double LST, vector<cartVec>& baselines, vector< vector<double> >& PB, vector<pointSource>& allPointSources){
	vector<double> realSpaceDiagonalPart(nPointSources, 0.0);
	vector<cartVec> pointSoruceCartVecs;
	for (int n = 0; n < nPointSources; n++){
		horizPoint thisHorizPoint = allPointSources[n].point.toHoriz(LST);
		pointSoruceCartVecs.push_back(thisHorizPoint.toVec());
		realSpaceDiagonalPart[n] = primaryBeam(thisHorizPoint, PB);	
	} 

	vector< vector<complex> > Amatrix(nPointSources, vector<complex>(nBaselines, complex(0,0)));
	for (int n = 0; n < nPointSources; n++){
		for (int b = 0; b < nBaselines; b++){
			double argument = k * pointSoruceCartVecs[n].dot(baselines[b]);
			Amatrix[n][b].re = cos(argument) * realSpaceDiagonalPart[n];
			Amatrix[n][b].im = sin(argument) * realSpaceDiagonalPart[n];
		}
	}
	return Amatrix;
}

//This function computes KAtNinvA_pointsource for each snapshot and then adds them to the existing sum to get an overall PSF
void addSnapshotPointSourcePSF(vector< vector<double> >& pointSourcePSF, vector< vector<complex> >& KAtranspose, vector<double>& Ninv, vector< vector<complex> >& pointSourceAmatrix, vector<int>& mapOfIndicesInExtendedIndexVector){
	for (int n = 0; n < nPixels; n++){
		for (int p = 0; p < nPointSources; p++){
			for (int b = 0; b < nBaselines; b++){
				pointSourcePSF[n][p] = pointSourcePSF[n][p] + 2 *(pointSourceAmatrix[p][b] * KAtranspose[mapOfIndicesInExtendedIndexVector[n]][b].conj() * Ninv[b]).re;
			}
		}
	}
}


//This function normalizes the PSF to 1 at the desired point probed and it computes the diagonal normalization matrix needed to do that.
vector<double> computeNormalizationAndNormalizePSF(vector< vector<double> >& PSF, vector<int>& mapOfIndicesInExtendedIndexVector){
	cout << "Now normalizing the PSF..." << endl;
	vector<double> normalizationMatrix(nPixels,0.0);
	for (int n = 0; n < nPixels; n++){
		if (PSFpeaksAtOne) {
			double maxPSFValueInThisRow = 0.0;
			for (int nx = 0; nx < nPixelsExtended; nx++){
				if (PSF[n][nx] > maxPSFValueInThisRow) maxPSFValueInThisRow = PSF[n][nx];
			}
			normalizationMatrix[n] = 1.0/maxPSFValueInThisRow;
		} else {
			normalizationMatrix[n] = 1.0/PSF[n][mapOfIndicesInExtendedIndexVector[n]];	
		}	
	}
	for (int n = 0; n < nPixels; n++){
		for (int nx = 0; nx < nPixelsExtended; nx++) PSF[n][nx] = normalizationMatrix[n]*PSF[n][nx];
	}
	return normalizationMatrix;
}

//Apply the normal PSF renormalization to the point source PSF
void applyNormalizationToPointSourcePSF(vector<double>& normalizationMatrix, vector< vector<double> >& pointSourcePSF){
	for (int n = 0; n < nPixels; n++){
		for (int p = 0; p < nPointSources; p++) pointSourcePSF[n][p] = normalizationMatrix[n]*pointSourcePSF[n][p];
	}
}

//This function renormalizes the dirty map so that it is the true map convolved with the PSF (on average).
vector<double> renormalizeMap(vector<double>& normalizationMatrix, Healpix_Map<double>& coaddedMap, vector<int>& healpixIndices){
	vector<double> renormalizedMap(nPixels, 0.0);
	for (int n = 0; n < nPixels; n++) renormalizedMap[n] = coaddedMap[healpixIndices[n]] * normalizationMatrix[n];
	return renormalizedMap;
}

//This function computes the nosie covariance between pixels in the facet.
vector< vector<double> > computeNoiseCovariance(vector< vector<double> >& PSF, vector<double>& normalizationMatrix, vector<int>& mapOfIndicesInExtendedIndexVector){
	vector< vector<double> > noiseCovariance(nPixels, vector<double>(nPixels,0.0));
	for (int n1 = 0; n1 < nPixels; n1++){
		for (int n2 = 0; n2 < nPixels; n2++){
			noiseCovariance[n1][n2] = PSF[n1][mapOfIndicesInExtendedIndexVector[n2]] * normalizationMatrix[n2];
		}
	}
	return noiseCovariance;
}

//This function convolves the "true sky" with the PSF and compares that to the final map
void reportError(vector<double>& renormalizedMap, vector< vector<double> >& PSF){
	vector<double> trueSky(nPixelsExtended, 0.0);
	fstream infile(trueSkyFilename.c_str(),fstream::in);
	double trueSkyPixel = 0;
	for (int n = 0; n < nPixelsExtended; n++){
		infile >> trueSkyPixel;
		trueSky[n] = trueSkyPixel;
	}
	infile.close();
	vector<double> convolvedSky(nPixels, 0.0);
	for (int i = 0; i < nPixels; i++){
		for (int j = 0; j < nPixelsExtended; j++){
			convolvedSky[i] += PSF[i][j] * trueSky[j];
		}
	}
	double normOfMap = 0;
	double normOfDifference = 0; 
	for (int i = 0; i < nPixels; i++){
		normOfMap += pow(renormalizedMap[i],2);
		normOfDifference += pow(renormalizedMap[i] - convolvedSky[i],2);
	}
	double error = sqrt(normOfDifference) / sqrt(normOfMap);
	cout << "Map - Convolved Sky Error is " << error << endl << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){
	
	cout << endl << "Running Faceted Mapmaking..." << endl << endl;
	//Load relevant quanties and data
	if (argc == 2) specsFilename = argv[1];
	loadSpecs();
	vector<int> baselineRedundancy;
	vector<cartVec> baselines = loadBaselines(baselineRedundancy);
	vector<cartVec> arrayPrincipalAxes = loadArrayPrincipalAxes();
	vector<cartVec> projectedBaselines = calculateProjectedBaselines(baselines,arrayPrincipalAxes);
	vector<double> LSTs; //LSTs of each integration
	vector< vector<complex> > allVisibilities = loadVisibilities(baselines, LSTs);
	vector< vector<double> > noiseVarianceOnEachVisibiltiy = loadNoiseVarianceOnEachVisibiltiy(baselines, LSTs, baselineRedundancy);
	vector<horizPoint> PBpointings = loadPBPointings(LSTs);
	vector< vector<double> > discretizedPrimaryBeam = loadDiscretizedPrimaryBeam();

	//Geometric calculations 
	vector<horizPoint> facetCenterPointings = computeFacetCenterPointings(LSTs);
	vector<bool> integrationsToUse = determineIntegrationsToUse(PBpointings, facetCenterPointings);
	if (nIntegrationsToUse == 0) return 1;
	vector< vector<int> > snapshotLSTindices = assignLSTsToSnapshots(integrationsToUse);
	cout << "LSTs range from " << LSTs[snapshotLSTindices[0][0]] << " to " << LSTs[snapshotLSTindices[snapshotLSTindices.size()-1][snapshotLSTindices[snapshotLSTindices.size()-1].size()-1]] << "." << endl;
	
	vector<int> healpixIndices = computeHealpixIndices(false);
	vector<equaPoint> pixelEquaPointings = computeEquaPointings(healpixIndices);
	Healpix_Map<bool> mapOfPixelsInFacet = isInFacet(healpixIndices);
	exportPixels(healpixIndices,healpixPixelFilename);
	exportCoordinates(pixelEquaPointings, pixelCoordinatesFilename);
	
	vector<int> extendedHealpixIndices = computeHealpixIndices(true);
	vector<equaPoint> extendedPixelEquaPointings = computeEquaPointings(extendedHealpixIndices);
	Healpix_Map<bool> mapOfPixelsInExtendedFacet = isInFacet(extendedHealpixIndices);
	exportPixels(extendedHealpixIndices,extendedHealpixPixelFilename);
	exportCoordinates(extendedPixelEquaPointings, extendedPixelCoordiantesFilename);
	vector<int> mapOfIndicesInExtendedIndexVector = mapIndicesFromFacetToExtendedFacet(healpixIndices, extendedHealpixIndices);

	//Load Point Sources
	vector<pointSource> allPointSources;
	if (alsoComputePointSourcePSF || overwriteVisibilitiesWithPointSources) loadPointSources(allPointSources, discretizedPrimaryBeam);

	//Overwrite Visibilities
	if (overwriteVisibilitiesWithGSM) overwriteVisibilitiesWithTheGSM(allVisibilities, baselines, discretizedPrimaryBeam, LSTs, facetCenterPointings, extendedHealpixIndices, extendedPixelEquaPointings, integrationsToUse, mapOfPixelsInFacet);
	if (overwriteVisibilitiesWithPointSources) overwriteVisibilitiesWithAListOfPointSources(allVisibilities, baselines, discretizedPrimaryBeam, LSTs, facetCenterPointings, allPointSources, integrationsToUse, extendedHealpixIndices);

	//Rephase and renormalize
	convertToTemperatureUnits(allVisibilities, integrationsToUse);	
	noiseInverseVarianceWeight(allVisibilities, noiseVarianceOnEachVisibiltiy, integrationsToUse);
	rephaseVisibilities(allVisibilities, snapshotLSTindices, LSTs, baselines, facetCenterPointings);

	//Loop over all snapshots to make map and PSF contributions
	cout << "Now calculating the map and the PSF..." << endl;
	Healpix_Map<double> coaddedMap = emptyHealpixMap();
	vector< vector<double> > PSF(nPixels, vector<double>(nPixelsExtended,0.0));
	vector< vector<double> > pointSourcePSF(nPixels, vector<double>(nPointSources,0.0));
	for (int n = 0; n < nSnapshots; n++){
		cout << " " << floor(100.0 * n / nSnapshots) << "% done. \r" << std::flush;
		int snapshotCentralLSTindex = snapshotLSTindices[n][int(round(snapshotLSTindices[n].size()/2.0-.5))];
		vector<complex> NinvTimesy = calculateNinvTimesy(allVisibilities, snapshotLSTindices[n]);
		vector<double> Ninv = calculateNinv(noiseVarianceOnEachVisibiltiy, snapshotLSTindices[n]);
		vector< vector<complex> > KAtranspose = calculateKAtranspose(LSTs[snapshotCentralLSTindex], extendedPixelEquaPointings, baselines, discretizedPrimaryBeam);
		addSnapshotMap(coaddedMap, NinvTimesy, KAtranspose, mapOfIndicesInExtendedIndexVector, healpixIndices);
		addSnapshotPSF(PSF, KAtranspose, Ninv, mapOfIndicesInExtendedIndexVector);
		if (alsoComputePointSourcePSF){
			vector< vector<complex> > pointSourceAmatrix = calcualtePointSourceAmatrix(LSTs[snapshotCentralLSTindex], baselines, discretizedPrimaryBeam, allPointSources);
			addSnapshotPointSourcePSF(pointSourcePSF, KAtranspose, Ninv, pointSourceAmatrix, mapOfIndicesInExtendedIndexVector);
		}
	}

	cout << "Done.                  " << endl;  
	vector<double> normalizationMatrix = computeNormalizationAndNormalizePSF(PSF, mapOfIndicesInExtendedIndexVector);
	exportVector(normalizationMatrix,nPixels,DmatrixFilename);

	//Compute and save the final data products
	exportMatrix(PSF, nPixels, nPixelsExtended, PSFfilename);
	vector<double> renormalizedMap = renormalizeMap(normalizationMatrix, coaddedMap, healpixIndices);
	exportVector(renormalizedMap, nPixels, finalMapFilename);
	vector< vector<double> > noiseCovariance = computeNoiseCovariance(PSF, normalizationMatrix, mapOfIndicesInExtendedIndexVector);
	exportMatrix(noiseCovariance, nPixels, nPixels, noiseCovFilename);
	if (alsoComputePointSourcePSF){
		applyNormalizationToPointSourcePSF(normalizationMatrix, pointSourcePSF);
		exportMatrix(pointSourcePSF, nPixels, nPointSources, pointSourcePSFFilename);
	}

	reportError(renormalizedMap, PSF);
	cout << "Done making a map. " << endl << endl;
	return 0;
}
