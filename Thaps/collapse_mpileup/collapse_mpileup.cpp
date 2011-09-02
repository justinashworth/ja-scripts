// parses samtools mpileup file, adds up counts on fwd and rvs strands, outputs these to separate files
// (custom solution that probably replicates some other existing one that I am not aware of ;)
// faster than Python (nice for huge files)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>

typedef std::string string;
typedef std::ifstream ifstream;
typedef std::ofstream ofstream;

class FilterSum {
public:
	FilterSum() : sum_(0) {}
	FilterSum(string prefix, string suffix, string pattern)
		: sum_(0), prefix_(prefix), suffix_(suffix), pattern_(pattern) {}

	// explicit copy constructor and assignment operator necessary due to ofstream member
	FilterSum( FilterSum const & f ) {
		*this = f;
	}
	FilterSum & operator= (FilterSum const & f ) {
		sum_ = f.sum_;
		prefix_ = f.prefix_;
		suffix_ = f.suffix_;
		pattern_ = f.pattern_;
		return *this;
	}

	~FilterSum() { of_.close(); }
	bool match(string const & match) const {
		return ( match.find(pattern_) != string::npos );
	}

	void add( int val ) { sum_ += val; }

	// open file if necesssary, give access to ofstream
	ofstream & out() {
		if ( !of_.is_open() ) {
			string name(prefix_ + suffix_);
			of_.open(name.c_str());
		}
		return of_;
	}

	void zero() { sum_ = 0; }
	int sum() const { return sum_; }

private:
	int sum_;
	string prefix_;
	string suffix_;
	string pattern_;
	std::ofstream of_;
};

//chr_1	1080	A	0	*	*	1	^F.	R	0	*	*	0	*	*	0	*	*	0	*	*	0	*	*
//chr_1	1081	T	0	*	*	1	.	[	0	*	*	0	*	*	0	*	*	0	*	*	0	*	*

void parse_file(string const & name, std::vector<FilterSum> & filters) {
	ifstream infile;
	infile.open(name.c_str());
	string line;
	while( getline(infile,line) ) {
		std::istringstream is(line);

		// get first three fields
		string seq; int pos; string base;
		is >> seq;
		is >> pos;
		is >> base;

		// get remaining fields in groups of three (counts,match,quality)
		while (is) {
			int counts; string match; string quality;
			is >> counts;
			is >> match;
			is >> quality;
			// each filter will sum matched counts on current line
			for (std::vector<FilterSum>::iterator f(filters.begin()); f != filters.end(); ++f) {
				if ( !f->match(match) ) continue;
				f->add(counts);
			}
		}

		// now output results of each filter
		for (std::vector<FilterSum>::iterator f(filters.begin()); f != filters.end(); ++f) {
			int sum( f->sum() );
			if ( sum == 0 ) continue;
			f->out() << seq << " " << pos << " " << sum << std::endl;
			// zero filter sum
			f->zero();
		}
	}
}

int main(int argc, char* argv[]) {

	std::vector<string> input;
	for (int i(1); i<argc; ++i) {
		string const & arg( argv[i] );
		// (option parsing would go here)
		input.push_back(arg);
	}

	string prefix( input.front() );

	// FilterSum class for cleanly specifying different symbols to match
	std::vector<FilterSum> filters;
	filters.push_back( FilterSum(prefix,".fwd",".") );
	filters.push_back( FilterSum(prefix,".rvs",",") );

	for (std::vector<string>::const_iterator in(input.begin()); in != input.end(); ++in) {
		parse_file(*in,filters);
	}

	return 0;
}
