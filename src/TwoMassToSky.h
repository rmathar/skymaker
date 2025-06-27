#pragma once

#include <string>

/*
* $Header: https://svn.mpia.de/gulli/geirs/src/trunk/extern/geirs2Panic/TwoMassToSky.h 781 2019-03-04 14:12:40Z mathar $
*/

class TwoMassToSky {
public:
	/** 2-sided field-of-view, expressed as the edge of
	* the full detector in units of pixels. We only consider
	* quadratic detectors. For a single Hawaii2 or Hawaii2RG chip, this is
	* 2048, for example.
	*/
	int detsize ;

	/** Pixel scale in the FITS image [rad/px].
	* Obtained from the [arcsec/px] number by multiplication with pi/(180*3600).
	*/
	float px ;

	/** Directory of 2MASS catalog with the xxx/t*.cat * files.
	* For example "tmc1" indicates that tmc1/000/t*.cat up to tmc1/189/t*.cat
	* are the files that have been processed by the tmcat expander.
	*/
	std::string catDir ;

	TwoMassToSky(int detsiz, float pixsc, const std::string catdir) ;

	void skymakeList(double ra, double decl, float maglim, char band, float q) const ;

	void wcs(double ra, double decl) const ;

	void spdRange(double decl, int spd[2]) const ;

	float mag2mag(float magI, int magidx, float q=1.) const ;

	static double hex2deg(const std::string & hexstr, bool isdeg) ;
protected:
	bool tanProj(double ra, double decl, double radec[2], int xy[2]) const ;

private:
} ;
