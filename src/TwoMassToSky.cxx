/** @file
* TwoMassToSky extracts star positions and magnitudes from the
* 2MASS catalogue (on the user's file system) and generates
* an ASCII format of the stars distributed over the pixels
* in the field of view in the catalogue style of skymaker.
*
* For interactive infrequent access to basically the same
* features see the <a href="http://irsa.ipac.caltech.edu/applications/2MASS/IM/interactive.html">2MASS Image Service</a>.
*/

#include "config.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <cmath>
#include <cstdlib>

#include "TwoMassToSky.h"

using namespace std;

// #define DEBUG

/** Constructor.
* @param[in] detsiz Number of pixels of detector in horizontal and vertical direction 
* @param[in] pixsc Pixel scale in radians (on the sky) per pixel.
* @param[in] catdir Directory with the 2MASS catalog.
* @author Richard J. Mathar
* @since 2012-11-29
*/
TwoMassToSky::TwoMassToSky(int detsiz, float pixsc, const std::string catdir) 
	: detsize(detsiz), px (pixsc), catDir(catdir)
{
}

/** Create a Skymaker list of the field-of view on stdout.
* @param[in] ra Right ascension in the center of the chip [rad]
* @param[in] decl Declination in the center of the chip [rad]
* @param[in] maglim The minimum limiting magnitude.
*  Stars that are fainter are not added to the list of the output.
* @param[in] band One of "J", "H" or "K"
* @param[in] q quantum efficiency in the range 0 to 1.0.
* @author Richard J. Mathar
* @since 2012-11-29
* @since 2013-02-13 with additional geirs parameter
* @since 2025-06-26 adapted to the format of https://irsa.ipac.caltech.edu/2MASS/download/allsky/
*/
void TwoMassToSky::skymakeList(double ra, double decl,float maglim, char band, float q) const
{
	/* list of all possible bsc_??? file names of the point source catalog
	*/
	static const char * fnam[] = 
	{
	"aaa","aab","aac","aad","aae","aaf","aag","aah","aai","aaj","aak","aal","aam",
	"aan","aao","aap","aaq","aar","aas","aat","aau","aav","aaw","aax","aay","aaz",
	"aba","abb","abc","abd","abe","abf","abg","abh","abi","abj","abk","abl","abm",
	"abn","abo","abp","abq","abr","abs","abt","abu","abv","abw","abx","aby","abz",
	"aca","acb","acc","acd","ace",
	"baa","bab","bac","bad","bae","baf","bag","bah","bai","baj","bak","bal","bam",
	"ban","bao","bap","baq","bar","bas","bat","bau","bav","baw","bax","bay","baz",
	"bba","bbb","bbc","bbd","bbe","bbf","bbg","bbh","bbi"
	} ;

	/* decide which of the magnitudes, left to right, of the 2MASS catalogue
	* is to be considred for the limiting magnitude and output.
	*/
	int magidx ;
	switch( band)
	{
	case 'J':
		magidx = 0 ; break;
	case 'H':
		magidx = 1 ; break;
	case 'K':
	default:
		magidx = 2 ; break;
	}
	

	/* determine a spread over directories equivalent to the
	* fov up and down half a detector size from the DEC coordinate
	*/
	int spdLim[2] ;
	spdRange(decl, spdLim) ;

#pragma omp parallel for
	/* loop over the input files bsc_<file>
	*/
	for(int fco = 0 ; fco < static_cast<int> (sizeof(fnam)/sizeof(const char *)) ; fco++)
	{
		std::ostringstream fil ;
		fil << catDir << "/psc_" << fnam[fco] ;
#ifdef DEBUG
		std::cout << fil.str() << std::endl ;
#endif
		/* loop over contents of the file
		*/
		ifstream infil(fil.str().c_str()) ;
		if ( !infil)
			/* could not open..*/
			;
		else
		{
			/* loop over all lines in the file
			*/
			while ( !infil.eof() )
			{
				string lin ;
				getline(infil,lin) ;
#ifdef DEBUG
				std::cout << lin << endl;
#endif
				/* replace vertical bars by blanks; the elements we're interestede
				* in are always with the first roughly 110 bytes of a line.
				*/
				//char linstr[164] ;
				//sscanf(lin.c_str(), "%163s",linstr) ;
				//linstr[164] = '\0' ;

				std::replace(lin.begin(), lin.end(), '|', ' ') ;
				
				std::istringstream linstr(lin) ;

				/* ra (degrees) and dec (degrees)
				*/
				double radec[2] ;

				/* name, paraphrasing ra and dec in hours 0..24 and degrees -90..+90
				*/
				string dummy ;

				/* J H K magnitude
				*/
				float mags[3] ;

				/* The two flags strings are currently ignored.
				* First RA (degrees) and dec (degrees)
				*/
				linstr >> radec[0] >> radec[1] ;
				for(int skip=0 ; skip < 4 ; skip++)
					linstr >> dummy ;
				linstr >> mags[0] ;
				for(int skip=0 ; skip < 3 ; skip++)
					linstr >> dummy ;
				linstr >> mags[1] ;
				for(int skip=0 ; skip < 3 ; skip++)
					linstr >> dummy ;
				linstr >> mags[2] ;

				/* preselect by magnitude. Admit only mags[] which are smaller (arithmetically)
				* than maglim.
				*/
				bool admi = (mags[magidx] <= maglim) ? true : false;

				/* tangential projection; middle is ra, decl.
				*/
				if ( admi)
				{
					/* convert 2MASS object coordinates from degrees to radians
					*/
					radec[0] *= M_PI/180.0 ;
					radec[1] *= M_PI/180.0 ;

					/* xy[] are the FITS coordinates determined by tanProj() */
					int xy[2] ;
					admi = tanProj(ra,decl,radec,xy) ;
					if ( admi)
					{
						/* convert infrared magnitude to a virtual V magnitude
						*/
						const float vmag = mag2mag(mags[magidx],magidx,q) ;
						/* because the readobj() function of Skymaker ignores comment lines,
						* we insert the 2MASS catalog entry here, too.
						*/
#pragma omp critical
						cout << "# " << lin << endl 
							<< "100 " << xy[0] << " " << xy[1] << " " << vmag << endl ;
					}
				}
			}
		}
	}
} /* skymakeList */

/** Convert a sedecimal-String to degrees.
* @param[in] hexstr A string with optional sign, integer number, colon, integer number, colon
*   and integer or floating point number.
* @param[in] isdeg If true, assume that the string is in DD:MM:SS.ss format.
*   If false the string is in HH:MM:SS.ss format.
*   This implies an additional factor 15 to move on to degrees.
*   If the string contains a sign (plus or minus), this parameter is ignored
*   and the function assumes degrees.
* @return The value in degrees.
* @since 2012-11-13
*/
double TwoMassToSky::hex2deg(const string & hexstr, bool isdeg)
{
	/* find first colon, last colon, and optional sign */
	size_t fcol = hexstr.find_first_of(":") ;
	size_t lcol = hexstr.find_last_of(":") ;
	size_t sig = hexstr.find_first_of("+-") ;

	/* If a sign was found before the first colon, assume that the argument
	* was given in degrees, not hours.
	*/
	if ( sig != string::npos && sig < fcol)
		isdeg = true ;

	double deg ;
	// if ( fcol != string::npos && lcol != string::npos && dot != string::npos )
	if ( fcol != string::npos && lcol != string::npos )
	{
		int dm[2] ;

		/* first (optionally signed) number between start and fcol
		*/
		sscanf(hexstr.substr(0,fcol).c_str(),"%d",&dm[0]) ;

		/* second number between fcol and lcol
		*/
		sscanf(hexstr.substr(fcol+1,lcol-fcol-1).c_str(),"%d",&dm[1]) ;

		double s ;
		/* third number after lcol
		*/
		sscanf(hexstr.substr(lcol+1).c_str(),"%lf",&s) ;

		/* unsigned absolute value, hours or degrees
		*/
		deg  = fabs(dm[0]) + dm[1]/60.0 + s/3600.0 ;
		if ( ! isdeg)
			/* convert hours to degrees */
			deg *= 15.0 ;
		if ( dm[0] < 0)
			/* add sign back to result */
			deg *= -1 ;
	}
	else
	{
		/* if no or only 1 colon is found, assume this is a simple floatign point number, in degrees or hours
		*/
		sscanf(hexstr.c_str(),"%lf",&deg) ;
		/* convert from hours to degrees */
		if ( !isdeg)
			deg *= 15.0 ;
	}
	return deg ;
} /* hex2deg */


/** Convert postition into an index of the catalog
* @param[in] decl The pointing (declination) in the middle of the fov [rad].
*  The valid range is -1.57..1.57.
* @param[out] spd The lower and upper limit (upper exclusive) of the directory name.
*  Both in the range 0 to 180.
* @since 2012-12-10
*/
void TwoMassToSky::spdRange(double decl, int spd[2] ) const
{
	/* declination at the lower rim of the detector in radians. Note that we are only
	* illuminating non-rotated images with S=down, N= up, so a simple subtraction
	* suffices.
	*/
	double dec = decl-0.5*detsize*px ;
	/* 2MASS directories count from 0 to 179, which is dec in degrees with a bias of 90
	*/
	dec *= 180.0/M_PI ;
	dec += 90.0 ;
	/* the area may include one of the poles.. to wrap around correctly use a max() here
	*/
	spd[0] = max(0, (int)dec) ;

	dec = decl+0.5*detsize*px ;
	dec *= 180.0/M_PI ;
	dec += 90.0 ;
	spd[1] = min(180, 1+(int)dec) ;
}

/** Convert a sky coordinate (pair) to a pixel coordinate (pair) in tangential projection.
* @param[in] ra RA in the middle of the plate [rad]
* @param[in] decl Declination in the middle of the plate [rad]
* @param[in] radec The ra and dec position of the object on the sky [rads]
* @param[out] xy Set to the FITS image coordinates if the return value is true.
* @return true If the coordinate falls into the rectangular window of the detector.
*/
bool TwoMassToSky::tanProj(double ra, double decl, double radec[2], int xy[2] ) const
{
	/* If star coordinates are ra_star and dec_star and the plate center
	* is at ra_pl and dec_pl, the 3-dimensional cartesian coordinates of these two
	* positions in the topocentric system with the observer in the center are
	* x(star) = cos(ra_star)*cos(dec_star)
	* y(star) = sin(ra_star)*cos(dec_star)
	* z(star) = sin(dec_star)
	* x(pl) = cos(ra_pl)*cos(dec_pl)
	* y(pl) = sin(ra_pl)*cos(dec_pl)
	* z(pl) = sin(dec_pl)
	* The angular distance between these is by the vector dot product
	* cos(dist) = sin(dec_star)*sin(dec_pl)+cos(dec_star)*cos(dec_pl)*cos(ra_star-ra_pl)
	* This is angular distance between radec and the middle of the plate (radians)
	*/
	double andot = sin(radec[1])*sin(decl)+cos(radec[1])*cos(decl)*cos(radec[0]-ra) ;

	/* In the middle of the plate, the unit direction in alpha is the
	* partial derivate with respect to alpha,
	* d(plate)/dalpha = [-sin(ra_pl)*cos(dec_pl),cos(ra_pl)*cos(dec_pl),0] ;
	* Normalized this becomes
	* d(plate)/dalpha = [-sin(ra_pl),cos(ra_pl),0] =unit(alpha)

	* In the middle of the plate, the unit direction in delta is the
	* partial derivate with respect to delta, which has unit length:
	* d(plate)/ddelta = [-cos(ra_pl)*sin(dec_pl),-sin(ra_pl)*sin(dec_pl),cos(dec_pl)] =unit(delta).

	* These two partial derivatives span the flat coordinate space onto which
	* the star positions are projected.

	* The direction to the star is split into a vector to the plate center plus the two steps
	* along  these unit directions plus some (discarded) step back onto the celestial sphere
	* in the unit direction of the star:
	* star = plate + dalpha * unit(alpha) +ddelta*unit(delta)- t*star ;
	* where t>0 because walking in the tangential plane means increasing the
	* distance to the observer.
	* (1+t)*star = plate + dalpha * unit(alpha) +ddelta*unit(delta).
	* Multiply with the plate center vector and use orthogonality:
	* (1+t)*star*plate = 1.
	* =(1+t)*andot = oneplust*andot calculated above.
	*/
	double oneplust= 1./andot ;

	/* There may be cases where the stars are in the back of the observer,
	* where t>1, 1+t>2, which are not to be used...
	*/
	if ( oneplust >= 2.)
		return false;

	/* Take vector dot product with unit(alpha) and use orthogonality:
	* (1+t)*star*unit(alpha) = dalpha
	* The dot product between star and unit(alpha) is evaluated in cartesian coords.
	double dalpha = oneplust*(-cos(radec[0])*cos(radec[1])*sin(ra) + sin(radec[0])*cos(radec[1])*cos(ra) ) ;
	double dalpha = oneplust*cos(radec[1])*(-cos(radec[0])*sin(ra) + sin(radec[0])*cos(ra) ) ;
	*/
	const double dalpha = oneplust*cos(radec[1])*sin(radec[0]-ra) ;

	/* Multiply with unit(delta) and use orthogonality:
	* (1+t)*star*unit(delta) = ddelta
	const double ddelta = oneplust*(-cos(radec[0])*cos(radec[1])*cos(ra)*sin(decl) 
			- sin(radec[0])*cos(radec[1])*sin(ra)*sin(decl)
			+ sin(radec[1])*cos(decl)
			) ;
	*/
	const double ddelta = oneplust*(-cos(radec[1])*sin(decl) *cos(radec[0]-ra)
			+ sin(radec[1])*cos(decl)
			) ;

	/* dalpha and ddelta are the two cartesian coordinaes in units of radians.
	* Divide this through the pixel scale (also in rads/pix) to convert to pixels.
	*/
	double delt[2] ;
	delt[0] = dalpha/px ;
	delt[1] = ddelta/px ;

	/* Catch the cases where conversion to integer further down
	* leads to an overflow...
	*/
	if ( fabs(delt[0]) > detsize/2 || fabs(delt[1]) > detsize/2 )
		return false;

	xy[0] = (int)delt[0] ;
	xy[1] = (int)delt[1] ;
	if ( abs(xy[0]) <= detsize/2 && abs(xy[1]) <= detsize/2 )
	{
		/* flip ra positive axis to the left. Add 1 to obtain the 1-based FITS coords */
		xy[0] = -xy[0] + 1 + detsize/2 ;
		xy[1] = xy[1] + 1+detsize/2 ;
		return true;
	}
	else
		return false;
} /* tanProj */

/** Print the associated WCS keywords to the standard error output.
* The FITS image that is created by this invocation of the Skymaker
* is oriented with North up and East to the left. There is currently
* no framework to create a more general sky rotation implied by some
* instrument specific optics.

* In practise this means that the command line supports specification
* of the FITS images by equatorial coordinates, but not by some mix
* of altitudes, azimuths or hour angles and geographic latitudes
* and similar sets of parameters.

* @param[in] ra Right ascension in the center of the detector [rad]
* @param[in] decl Declination in the center of the detector [rad]
*/
void TwoMassToSky::wcs(double ra, double decl) const
{
	cerr << "CUNIT1 = 'deg' / [] unit along FITS axis 1" << endl ;
	cerr << "CUNIT2 = 'deg' / [] unit along FITS axis 2" << endl ;
	cerr << "CTYPE1 = 'RA---TAN' / [] tangential projection det. plane" << endl ;
	cerr << "CTYPE2 = 'DEC--TAN' / [] tangential projection det. plane" << endl ;
	cerr << "CRVAL1 = " << (ra*180.0/M_PI) << " / [deg] RA center of det .plane " << endl ;
	cerr << "CRVAL2 = " << (decl*180.0/M_PI) << " / [deg] DEC center of det .plane " << endl ;
	cerr << "CRPIX1 = " << ((detsize+1.)/2.) << " / [px]" << endl ;
	cerr << "CRPIX2 = " << ((detsize+1.)/2.) << " / [px]" << endl ;
	cerr << "CD1_2 = 0 / [deg/px] outer diagonal of trans. matrix" << endl ;
	cerr << "CD2_1 = 0 / [deg/px] outer diagonal of trans. matrix" << endl ;
	cerr << "CD1_1 = " << (-px*180.0/M_PI) << " / [deg/px] diagonal of WCS matrix" << endl ;
	cerr << "CD2_2 = " << (px*180.0/M_PI) << " / [deg/px] diagonal of WCS matrix" << endl ;
	cerr << "DETSIZE = '[1:" << detsize << ",1:" << detsize << "]' / [px]" << endl ;
	cerr << "PIXSCAL = " << (px*180.0*3600.0/M_PI) << " / [arcsec/px]" << endl ;
#if 1
	cerr << "RA = " << (ra*180.0/M_PI) << " / [deg] RA center of det .plane " << endl ;
	cerr << "DEC = " << (decl*180.0/M_PI) << " / [deg] RA center of det .plane " << endl ;
#endif
}

/** Convert a H, J or Ks magnitude to a V magnitude.
* This is only intended to be used for mimicking magnitudes for input to skymaker.
* The conversion computes a number of photons for the specified magnitude in
* the infrared and takes this number as the number of photons in the visible to
* obtain an "equivalent" magnitude in the visible.
* @param[in] magI The magnitude in the infrared band.
* @param[in] magidx The integer value of the band. 0 for J, 1 for H and 2 for K.
* @param[in] q A flux reducing quantum efficiency in the range 0 to 1.0.
* @return magI converted to the visible.
* @since 2012-12-22
*/
float TwoMassToSky::mag2mag(float magI, int magidx, float q) const
{
	/* table of band, lambda (micron) and reference flux at magnitude zero:
	* V, 0.55, 3640; J, 1.26, 1600 ; H, 1.6, 1080 ; K, 2.22, 670
	*/
	static double fluxx[4][2] = { {1.26, 1600}, {1.6, 1080}, {2.22,670}, {0.55,3640}} ;

	/* mag(b) = -2.5*log10(Flux(b)/reffl(b))  where b is a band.
	* Flux(b) = photon_count*h*nu(b)/(area*time).
        * photon_count = Flux(b)*area*time/(h*nu(b)), and lambda(b)*nu(b)= velocity of light.
        * So by construction of the conversion,
        * Flux(b)*area*time*lambda(b)/(h*c) = photon_count = Flux(V)*area*time*lambda(V)/(h*c).
        * Flux(V) = Flux(b)*lambda(b)/lambda(V),
        * Flux(V)/reffl(V) = Flux(b)*lambda(b)/lambda(V)/reffl(V),
        * log10[Flux(V)/reffl(V)] = log10[Flux(b)*lambda(b)/lambda(V)/reffl(V)],
        * log10[Flux(V)/reffl(V)] = log10[Flux(b)/reffl(b)*lambda(b)*reffl(b)/lambda(V)/reffl(V)],
        * log10[Flux(V)/reffl(V)] = log10[Flux(b)/reffl(b)] + log10[lambda(b)*reffl(b)/lambda(V)/reffl(V)],
        * -2.5*log10[Flux(V)/reffl(V)] = -2.5*log10[Flux(b)/reffl(b)] -2.5* log10[lambda(b)*reffl(b)/lambda(V)/reffl(V)],
        * mag(V) = mag(b) -2.5* log10[lambda(b)*reffl(b)/lambda(V)/reffl(V)],
	*/
	return magI -2.5*log10(q*fluxx[magidx][0]*fluxx[magidx][1]/fluxx[3][0]/fluxx[3][1]) ;
} /* mag2mag */

/** Emit a short usage reminder.
* @param[in] argv0 The name of the executable.
* @author R. J. Mathar
* @since 2012-11-29
*/
void usage(char *argv0)
{
        cout << "Usage:" << endl ;
        cout << argv0 << " [-D 2massdir] [-r RA/h] [-d DEC/deg] [-p px] [-s arcs] [-m mag] [-b {J,H,K}] [-q qeff] > sky.list 2> sky.hdr" << endl ;
        cout << "Option -D is the directory name of the 2MASS catalog." << endl ;
        cout << "Option -r is the pointing right ascension in hours." << endl ;
        cout << "Option -d is the pointing declination in degrees." << endl ;
        cout << "Option -p defines the number of pixels per FITS axis." << endl ;
        cout << "Option -s defines the scale, arcseconds per pixel." << endl ;
        cout << "Option -m defines the limiting magnitude." << endl ;
        cout << "Option -b defines the infrared band." << endl ;
        cout << "Option -q defines a quantum efficiency of the detector." << endl ;
        cout << "Standard output goes into a file for Skymaker" << endl ;
        cout << "Standard error output goes into a FITS header file" << endl ;
} /* usage */

/** The C++ program TwoMassToSky converts portions of the 2MASS catalog to a Skymaker's list file.
* The prerequisites of running the program are the regions of interest of the 2MASS catalog in
* the standard layout in the file system, which are files named ../???/t*.cat, where the three
* question marks are the three digits of the quantized declination (measured from 0 of the
* southern pole up to 179).
*
* This means the program will scan these directories, and if some of the files or their
* lines are missing, the stars that are not found will not be produced by
* the program either.
*
* The standard output contains lines in the Skymaker format, a 100 followed by the
* two FITS pixel locations and a magnitude. The standard error contains a snapshot
* that would be added to the FITS file header of what will be produced by Skymaker
* to have a useful WCS system across the FITS image. See the <code>Makefile</code>
* for an example of chaining the operations.
*
* @param[in] argc
* @param[in] argv
*
* The command line options are:
*
* <code>-D</code> is followed by the location of the directory of the 2MASS catalog.
*  This is without the <code>???/t*.cat</code> portion of the file names.
*
* <code>-r</code> is followed by the right ascension (in hours from 0 to 24) of the
*  pointing in the center of the FITS plate.
*  Instead of a floating point number in hours, the RA may also be provided by
*  the standard two-colon hex-format, HH:MM:SS.ss.
*
* <code>-d</code> is followed by the declination (in degrees from -90 to 90) of the
*  pointing in the center of the FITS plate.
*  Instead of a floating point number in hours, the DEC may also be provided by
*  the standard two-colon hex-format, +-DD:MM:SS.ss.
*
* <code>-p</code> defines the number of the pixels along x and along y.
*  (We are only dealing with quadratic detector areas.)
*  The product of this with the pixel scale is the two-sided field of view in which
*  stars of the 2MASS catalogue must reside to be copied to the output.
*  Warning: this *must* be the same as the IMAGE_SIZE in the sky.conf file .
*
* <code>-s</code> is the pixel scale in units of arcseconds per pixel.
*  Warning: this *must* be the same as the PIXEL_SIZE in the sky.conf file .
*
* <code>-m</code> clips the magnitude (in the infrared band, not refering to the visible)
*  of the star list that is put to the output. A number of 8.5, for example, means that
*  only objects brighter than 8.5 (numerically smaller than 8.5) are copied over.
*
* <code>-b</code> is followed by a single capital letter, one out of three J, H or K as expected.
*
* <code>-q</code> is followed by a number between 0 and 1, which scales the magnitude
*    of the star in the output according to that quantum efficiency.
*
*
* The syntax in overview:
*
* <code>TwoMassToSky [-D 2massdir] [-r RA/h] [-d DEC/deg] [-p px] [-s arcs] [-m mag] [-b {J,H,K}] > sky.list.in 2> sky.hdr</code>
* <code>rm sky.fits</code>
* <code>sky sky.list.in -c sky.conf</code>
* <code>rm sky.list</code>
* <code>fedithead sky.fits sky.hdr # or fmodhead sky.fits'[0]' sky.hdr</code>
* <code>ds9 sky.fits</code>
* @author Richard J. Mathar
* @since 2012-11-29
*/
int main(int argc, char *argv[])
{
	/* pixel scale in arcseconds/pix
	* Default is CAHA/PANIC H4-RG.
	*/
	float pixsc = 0.445*15.0/18.0 ;

	/* full detector pixel size. Default is Hawaii-2 mosaic or Hawaii-4RG like CAHA/PANIC */
	int detsize = 2*2048 ;

	/* directory with the 2MASS catalog. Default is to search in a parallel
	* directory of the working directory.
	*/
	std::string catdir("../2mass/") ;

	/* RA in degrees.
	*/
	double ra = 10.0 ;

	/* DEC in degrees . Close to the southern pole by default, assuming
	* that at least the first subcatalog of 2MASS has been downloaded.
	*/
	double dec = -89.0 ;

	/* limiting magnitude. By default arbitrarily faint objects are copied.
	*/
	int mag = 40 ;

	/* Generic quantum efficiency assume to be optimistically 90 percent.
	*/
	float q = 0.9 ;

	/* Infrared band, either K, H or J, as definedby the 2MASS entries
	*/
	char band = 'K' ;

	/* command line option character. Loop over all of the comand line arguments.
	*/
	char oc ;

	string radec ;
	while ( (oc=getopt(argc,argv,"hD:r:d:p:s:m:b:")) != -1 )
	{
		switch(oc)
		{
		case 'h' :
			usage(argv[0]) ;
			return 0 ;
			break;
		case 'D' :
			catdir = string(optarg) ;
			break;
		case 'r' :
			radec= optarg ;
			/* convert from hours, minutes, seconds to degrees
			*/
			ra = TwoMassToSky::hex2deg(radec,false) ;
			break;
		case 'd' :
			radec= optarg ;
			/* convert from degrees, minutes, seconds to degrees
			*/
			dec = TwoMassToSky::hex2deg(radec,true) ;
			break;
		case 'p' :
			detsize = atoi(optarg) ;
			break;
		case 's' :
			/* pixel scale in arcseconds per pixel
			*/
			pixsc = atof(optarg) ;
			break;
		case 'q' :
			q = atof(optarg) ;
			break;
		case 'm' :
			mag = atoi(optarg) ;
			break;
		case 'b' :
			band = optarg[0] ;
			break;
		case '?' :
                        cerr << "Invalid command line option " << optopt << endl ;
                        break ;
		}
	}

	/* pixel scale in radians per pix */
	pixsc *= M_PI/(180.0*3600.0) ;

	TwoMassToSky twom(detsize, pixsc,catdir) ;

	/* ra converted to radians */
	ra *= M_PI/180.0 ;

	/* dec converted to radians */
	dec *= M_PI/180.0 ;

	twom.skymakeList(ra,dec,mag,band,q) ;

	twom.wcs(ra,dec) ;

	return 0 ;
} /* main */
