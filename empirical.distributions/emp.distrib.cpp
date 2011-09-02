// Justin Ashworth 2011
// Simple empirical distributions and p-value estimates for comparing genome-wide positions
// Explicit estimation of random results that does not suffer from the invalid assumptions and weaknesses of trying to make the problem fit a hypergeometric distribution
// such weaknesses of the hypergeometric assumption include: no unbiased way to account for actual distances between peaks, or false assumptions of significance to large sequence windows, or simply larger numbers of genome-wide positions
// see Huen and Russell (2010) BMC Bioinformatics doi:10.1186/1471-2105-11-359

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

class Region {
public:
	Region() : id_(""), start_(0), end_(0), contains_positions_(false) {}
	Region(std::string id, int start, int end) : id_(id), start_(start), end_(end), contains_positions_(false) {}
public:
	std::string id_;
	int start_, end_;
	bool contains_positions_;
};

std::ostream & operator << ( std::ostream & out, Region const & rgn )
{
	out << rgn.id_ << " " << rgn.start_ << " " << rgn.end_ << " " << rgn.contains_positions_;
	return out;
}

typedef std::vector<Region> Regions;

void regions_contain_positions(Regions & regions, std::vector<int> const & positions)
{
	for ( Regions::iterator it(regions.begin()); it != regions.end(); ++it ) {
		it->contains_positions_ = false;
		for ( std::vector<int>::const_iterator pos(positions.begin()); pos != positions.end(); ++pos ) {
			if ( *pos >= it->start_ && *pos <= it->end_) {
				it->contains_positions_ = true;
				break;
			}
		}
	}
}

void regions_contain_both_positions(
	Regions & regions,
	std::vector<int> const & pos1,
	std::vector<int> const & pos2
)
{
	for ( Regions::iterator rgn(regions.begin()); rgn != regions.end(); ++rgn ) {
		rgn->contains_positions_ = false;
		bool contains1(false);
		for ( std::vector<int>::const_iterator pos(pos1.begin()); pos != pos1.end(); ++pos ) {
			if ( *pos >= rgn->start_ && *pos <= rgn->end_) {
				contains1 = true;
				break;
			}
		}
		if ( ! contains1 ) continue;
		for ( std::vector<int>::const_iterator pos(pos2.begin()); pos != pos2.end(); ++pos ) {
			if ( *pos >= rgn->start_ && *pos <= rgn->end_) {
				rgn->contains_positions_ = true;
				break;
			}
		}
	}
}

void option_error( std::string option ) {
	std::cerr << "error option " << option << std::endl;
	exit(EXIT_FAILURE);
}

// random number functor for faster(?) random vector creation with std::generate_n()
struct GenRand {
	int range_;
public:
	GenRand(int r) : range_(r) {
	//	std::cerr << "GenRand with range " << range_ << std::endl;
	}
	int operator()() { return rand() % range_ + 1; }
};

int main( int argc, char *argv[] ) {

	std::string type("occupancy"),pos1filename,pos2filename,regionsfilename;
	int iters(1),seqlen(0);
	for (int i(1); i<argc; ++i) {
		std::string arg( argv[i] );
		if (arg == "-t") {
			if ( ++i >= argc ) option_error(arg);
			type = argv[i];
		} else if (arg == "-p1" || arg == "-p" ) {
			if ( ++i >= argc ) option_error(arg);
			pos1filename = argv[i];
		} else if (arg == "-p2") {
			if ( ++i >= argc ) option_error(arg);
			pos2filename = argv[i];
			type == "overlap";
		} else if (arg == "-r") {
			if ( ++i >= argc ) option_error(arg);
			regionsfilename = argv[i];
		} else if (arg == "-i") {
			if ( ++i >= argc ) option_error(arg);
			std::istringstream ss(argv[i]);
			ss >> iters;
		} else if (arg == "-s") {
			if ( ++i >= argc ) option_error(arg);
			std::istringstream ss(argv[i]);
			ss >> seqlen;
		}
	}

	std::ifstream pos1file;
	pos1file.open( pos1filename.c_str() );
	if (!pos1file) {
		std::cerr << "error, couldn't open positions(1) file " << pos1filename << std::endl;
		exit(EXIT_FAILURE);
	} else std::cerr << "opened file " << pos1filename << std::endl;

	std::vector<int> pos1,pos2;
	std::string line;
	while( getline(pos1file,line) ){
		int val(0);
		std::istringstream ss(line);
		ss >> val;
		pos1.push_back(val);
	}

	if ( type == "overlap" ) {
		std::ifstream pos2file;
		pos2file.open( pos2filename.c_str() );
		if (!pos2file) {
			std::cerr << "error, couldn't open positions(2) file " << pos2filename << std::endl;
			exit(EXIT_FAILURE);
		} else std::cerr << "opened file " << pos2filename << std::endl;

		while( getline(pos2file,line) ){
			int val(0);
			std::istringstream ss(line);
			ss >> val;
			pos2.push_back(val);
		}
	}

	// if seqlen not given, assume max position is approximately seqlen
	if ( seqlen <= 0 ) {
		std::cerr << "seqlen (-s) not provided--assuming highest position is seqlen (probably wrong, but could be approximately right)" << std::endl;
		// careful to avoid NULL pointer segfaults with iterator algos
		if ( pos1.size() > 0 ) { seqlen = std::max( seqlen, *std::max_element(pos1.begin(),pos1.end()) ); }
		if ( pos2.size() > 0 ) { seqlen = std::max( seqlen, *std::max_element(pos2.begin(),pos2.end()) ); }
	}
	std::cerr << "seqlen " << seqlen << std::endl;

	std::ifstream regionsfile;
	regionsfile.open( regionsfilename.c_str() );
	if (!regionsfile) {
		std::cerr << "error, couldn't open file " << regionsfilename << std::endl;
		exit(EXIT_FAILURE);
	} else std::cerr << "opened file " << regionsfilename << std::endl;

	Regions regions;
	while( getline(regionsfile,line) ){
		std::string id;
		int start(0),end(0);
		std::istringstream ss(line);
		ss >> id >> start >> end;
//		std::cerr << id << " " << start << " " << end << std::endl;
		regions.push_back( Region(id,start,end) );
	}

	// regions that contain pos1 positions
	if ( type == "occupancy" ) {
		regions_contain_positions(regions,pos1);
	} else if ( type == "overlap" ) {
		regions_contain_both_positions(regions,pos1,pos2);
	} else {
		std::cerr << "Unknown mode " << type << std::endl;
		exit(EXIT_FAILURE);
	}

	int numcont(0);
	for (Regions::iterator it(regions.begin()); it != regions.end(); ++it) {if(it->contains_positions_) ++numcont;}
	if ( type == "occupancy" ) {
		std::cerr << numcont << " regions contain positions" << std::endl;
	} else if ( type == "overlap" ) {
		std::cerr << numcont << " regions contain positions from both sets of positions" << std::endl;
	}

	// now estimate significance by simulating random distribution of equal number of positions in integer range
	std::vector<int> randcont;
	int sum(0);
	for ( int i(0); i < iters; ++i ) {
		std::vector<int> randpos1( pos1.size() );
		std::generate_n(randpos1.begin(), pos1.size(), GenRand(seqlen));
//		std::cerr << "randpos1: ";
//		std::copy(randpos1.begin(),randpos1.end(),std::ostream_iterator<int>(std::cerr," "));
//		std::cerr << std::endl;
		std::cerr << '.';
		if ( type == "overlap" && pos2.size() > 0 ) {
			std::vector<int> randpos2( pos2.size() );
			std::generate_n(randpos2.begin(), pos2.size(), GenRand(seqlen));
			regions_contain_both_positions(regions,randpos1,randpos2);
		} else {
			regions_contain_positions(regions,randpos1);
		}

		int numrand(0);
		for (Regions::iterator it(regions.begin()); it != regions.end(); ++it) {if(it->contains_positions_)++numrand;}
		randcont.push_back(numrand);
		sum += numrand;
	}
	std::cerr << std::endl;

	std::cerr << "number of regions containing positions in empirical random distributions:" << std::endl;
	std::copy(randcont.begin(),randcont.end(),std::ostream_iterator<int>(std::cerr," "));
	std::cerr << std::endl;
	std::cerr << "avg counts in " << iters << " random distributions: " << sum/randcont.size() << std::endl;
	int ngte(0);
	for ( std::vector<int>::iterator it(randcont.begin()); it != randcont.end(); ++it ) {
		if ( *it >= numcont ) ++ngte;
	}
	std::cerr << "fraction of random distributions with >= " << numcont << " counts: " << ((double)ngte)/randcont.size() << std::endl;

	std::copy(randcont.begin(),randcont.end(),std::ostream_iterator<int>(std::cout," "));
	std::cout << std::endl;

	return 0;
}
