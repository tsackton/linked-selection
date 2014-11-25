#include <iostream>
#include <fstream>
#include <vector>         
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;

/////

// This code is designed to compute the effective of background selection in reducing levels of genetic diversity in a species including 
/////// local recombination rates
/////// local functional density

//// it has been successfully compiled and tested using g++ version 4.2.1 on mac OSX 10.6.8
////	g++ -o compute_gk compute_gk.cpp

//// Note that this code requires input to be sorted by species, chromosome, and position
///// The included R script will produce the appropriate output to replicate our analysis but to reuse this code
///// for other purposes please be aware of this requirement

/// Format of input is: species, chrom, window.id (can be the position, treated as a string), 
/// genetic position of start of window (Morgans), genetic position of midpoint of window, genetic position of end of window
/// functional density (fraction of total coding bp in genome that are in window)

/// Output includes species, chrom, and window.id (passed unchanged from input)
// assumed mutation rate (passed as an input parameter), sh, P, and G for window

class cmd_line {
	
public:

	string species ; 
	long double u_mut ; 			/// mutation rate specified on the command line as U per diploid

	void read_cmd_line ( int argc, char *argv[] ) ; 
	
} ;

void cmd_line::read_cmd_line ( int argc, char *argv[] ) {
	
	for (int i=1; i<argc; i++) {
		if ( strcmp(argv[i],"-u") == 0 ) {
			u_mut = atof(argv[++i]) ;
		}
		if ( strcmp(argv[i],"-s") == 0 ) {
			species = argv[++i] ;
		}
	}
	
}	

//// these are the input data which are computed by an r script also provided with this manuscript
class line {
	
public:
	
	string species ;
	string chrom ;
	string window ; 
	long double mi ;		/// mi is the start, mk is the window center, ml is the window end
	long double mk ;
	long double ml ; 
	long double fd ; 
	
	void print_line() ; 
	friend istream& operator>>( istream &is, line &in ) ;
	
} ;

void line::print_line() { 
	
	cout << species << "\t" <<  chrom << "\t" <<  window ;
	
} 

istream& operator>>(istream &is, line &in ) {
	
	is >> in.species >> in.chrom >> in.window >> in.mi >> in.mk >> in.ml >> in.fd ;
	
}

class chromosome {
	
public:
	
	string chrom ; 
	vector<line> windows ; 
	
} ;

/// this will read in the data for each chromosome of the species specified on the command line
void read_lines( vector<chromosome> &chr, cmd_line &options ) { 
	 
	while ( cin ) {
		line new_line ; 
		cin >> new_line ; 
		if ( new_line.species == options.species ) { 
			bool exists = false ; 
			for ( int i = 0 ; i < chr.size() ; i ++ ) { 
				if ( chr.at(i).chrom == new_line.chrom ) { 
					chr.at(i).windows.push_back(new_line) ; 
					exists = true ; 
				}
			}
			if ( exists == false ) {
				chromosome new_chrom ; 
				new_chrom.chrom = new_line.chrom ; 
				new_chrom.windows.push_back(new_line) ; 
				chr.push_back(new_chrom) ; 
			}
		}
	}
}

/// using the windowed dataset we compute the parameter gk. 
void compute_gk ( vector<chromosome> &chr, cmd_line &options ) {
	
	/// across a range of selective coefficents
	for ( float ex = 1 ; ex <= 5 ; ex += .04 ) { 
		float sh = pow( 10, -1 * ex ) ; 
		
		/// outbreeding coefficient over range 0 to 1 by 0.04 increments
		for ( float p = 1 ; p >= 0 ; p -= .04 ) {

			/// compute these sets for each chromosome separately.  
			for ( int c = 0 ; c < chr.size() ; c ++ ) { 
				for ( int w = 0 ; w < chr.at(c).windows.size() ; w ++ ) { 
					long double gk = 0 ;
					for ( int i = 0 ; i < chr.at(c).windows.size() ; i ++ ) { 
						long double numerator = ( options.u_mut * sh * chr.at(c).windows.at(i).fd ) ;
						long double denominator = 2 * ( sh + p * abs( chr.at(c).windows.at(w).mk - chr.at(c).windows.at(i).mi ) ) * ( sh + p * abs(  chr.at(c).windows.at(w).mk - chr.at(c).windows.at(i).ml ) ) ; 

						//cout << options.u_best << "\t" << sh << "\t" << chr.at(c).windows.at(i).fd << "\t" << chr.at(c).windows.at(w).mk << "\t" << chr.at(c).windows.at(i).mi << "\t" << chr.at(c).windows.at(i+1).mi << endl ; 
						//cout << best_numerator << "\t" << denominator << endl ; 
						
						gk += numerator/denominator ;
					}
					
					chr.at(c).windows.at(w).print_line() ; 
					
					cout << "\t" << options.u_mut << "\t" << sh << "\t" << p << "\t" << gk << endl ; 
				}
			}
		}
	}
}

int main ( int argc, char **argv ) {
	
	cmd_line options ; 
	options.read_cmd_line( argc, argv ) ;

	vector<chromosome> chr ; 
	read_lines( chr, options ) ; 
	
	compute_gk( chr, options ) ; 
	
}
	


















