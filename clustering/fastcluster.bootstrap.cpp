#include <string> // for std::string
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <cmath>
#include <math.h>
#include <vector>
#include <cstdlib> // exit, EXIT_FAILURE
#include <algorithm>
#include <iterator>
#include <iomanip>

#include "fastcluster_add.cpp"
//#include "fastcluster_R.cpp"

typedef std::vector<t_float> Values;
typedef std::vector<Values> Matrix;
typedef std::vector<std::string> Labels;
typedef std::vector<size_t> Indices;
typedef std::vector<int> Ints;

static t_float const MAGIC_BIGNUM = 111e11;

struct HclustResult {
	Indices indices;
	Labels ids;
	Ints merge;
	Values height;
	Ints order;
};
typedef std::vector<HclustResult> HclustResults;

std::ostream & operator << (std::ostream & out, Values const & vals )
{
	for(Values::const_iterator it(vals.begin()); it!=vals.end(); ++it){
		out << *it << " ";
	}
	return out;
}

t_float
pearson_correlation(Values const & row1, Values const & row2, bool const pairwise_complete_obs)
{
	size_t N(row1.size());
	t_float EX(0), EY(0), EXX(0), EYY(0), EXY(0), x(0), y(0);
	bool matches(false);
	for(size_t i(0); i<N; ++i){
		x = row1[i];
		y = row2[i];
		if(pairwise_complete_obs && (x==MAGIC_BIGNUM || y==MAGIC_BIGNUM)) continue;
		matches=true;
		EX += x;
		EY += y;
		EXX += x*x;
		EYY += y*y;
		EXY += x*y;
	}
	if(!matches) return(2.0);
	else return (EXY - EX*EY/N) / sqrt( (EXX - EX*EX/N)*(EYY - EY*EY/N) );
}

bool is_number(const std::string& s) {
	if(s.empty()) return false;
	for(std::string::const_iterator it(s.begin()); it!=s.end(); ++it){
		if(std::isdigit(*it) || (it==s.begin() && *it=='-')) continue;
		return false;
	}
	return true;
}

void
read_matrix_file(
	std::string const & filename,
	Matrix & rr,
	Labels & ids
)
{
	// open input file
	std::ifstream file;
	file.open( filename.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open file " << filename.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cerr << "Reading file " << filename << std::endl;

	// read input file
	std::string line;
	bool firstline(true);
	while ( getline( file, line ) ) {
		// skip the first line, which should be a header
		if(firstline){ firstline=false; continue; }
		std::istringstream linestream(line);
		Values ratios;
		std::string token;
		bool id(true);
		while(getline(linestream, token, ' ')){
//			std::cerr << token << ":";
			// first token in line should be an id
			if(id){
				ids.push_back(token);
				id=false;
				continue;
			}
			// the rest should be numbers, but sometimes "NA"
			// set any "NA"s to a magic big number to keep track later
			t_float value;
			if(token=="NA") value=MAGIC_BIGNUM;
			else value=std::stof(token);
//			std::cerr << value << std::endl;
			ratios.push_back(value);
		}

		rr.push_back(ratios);
	}
}

void resample_distances(
	t_float * const dist,
	t_float * dist_resampled,
	Indices const & indices
){
	std::ptrdiff_t p(0);
	size_t N(indices.size());
	for(size_t i(0); i<(N-1); ++i){
		for(size_t j(i+1); j<N; ++j){
			int ri(indices[i]), rj(indices[j]);
			// because bootstrapping involves resampling with replacement, this can happen
			if(ri==rj){
				dist_resampled[p] = 0;
//				std::cout << p << "," << i << "," << j << "," << ri << "," << rj << ",(self)" << std::endl;
			}else{
				// tricky diagonal array indexing arithmetic for the D indices corresponding to Ns i and j
				size_t linear_index(rj-ri-1);
				for(size_t k(N-ri); k<N; ++k) linear_index += k;
//				std::cout << p << "," << i << "," << j << "," << ri << "," << rj << "," << linear_index << std::endl;
				dist_resampled[p] = dist[linear_index];
			}
			++p;
		}
	}
//		std::cout << std::endl;
}

void head_of_result(HclustResult const & result, size_t nprint=3)
{
	std::cout << "Indices: ";
	std::copy(result.indices.begin(), result.indices.begin()+nprint, std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	std::cout << "Merge: ";
	std::copy(result.merge.begin(), result.merge.begin()+nprint, std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	std::cout << "Height: ";
	std::copy(result.height.begin(), result.height.begin()+nprint, std::ostream_iterator<t_float>(std::cout, " "));
	std::cout << "(min " << *std::min_element(result.height.begin(), result.height.end()) << ", max " << *std::max_element(result.height.begin(), result.height.end()) << ")";
	std::cout << std::endl;
	std::cout << "Order: ";
	std::copy(result.order.begin(), result.order.begin()+nprint, std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
}

void
get_k_clusters(HclustResult const &, std::vector<Labels> &, size_t const)
{
	// not sure of the best/cleanest/most efficient way to do this
}

void output_clusters(HclustResult const & result, std::string prefix="")
{
	// K-based cluster selection: recursive agglomeration until k is reached

	size_t n(result.indices.size()), by(1);
	if(n<=100) by=10;
	else by=100;

	size_t end(n/by);
	for(size_t i(1); i<=end; ++i){
		size_t k(i*by);
		std::cout << k << " clusters..." << std::endl;

		std::vector<Labels> clusters;
		get_k_clusters(result,clusters,k);

		std::ofstream of;
		std::ostringstream fname;
		fname << prefix << "." << std::setw(6) << std::setfill('0') << k;
		of.open(fname.str().c_str());
		for(size_t cl(0); cl<clusters.size(); ++cl){
			std::copy(clusters[cl].begin(), clusters[cl].end(), std::ostream_iterator<std::string>(of, " "));
			of << '\n';
		}
		of.close();
	}

}

void output_bootstrap_results(HclustResults const & bootstraps, std::string prefix="bootstrap.")
{
	// outputs that can be read into R
	std::ofstream of;
	std::string fname(prefix+"indices");
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->indices.begin(), bs->indices.end(), std::ostream_iterator<int>(of, " "));
		of << '\n';
	}
	of.close();

	fname = prefix+"labels";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->ids.begin(), bs->ids.end(), std::ostream_iterator<std::string>(of, " "));
		of << '\n';
	}
	of.close();

	fname = prefix+"merges";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->merge.begin(), bs->merge.end(), std::ostream_iterator<int>(of, " "));
		of << '\n';
	}
	of.close();

	fname = prefix+"heights";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->height.begin(), bs->height.end(), std::ostream_iterator<t_float>(of, " "));
		of << '\n';
	}
	of.close();

	fname = prefix+"orders";
	of.open(fname.c_str());
	for(HclustResults::const_iterator bs(bootstraps.begin()); bs!=bootstraps.end(); ++bs){
		std::copy(bs->order.begin(), bs->order.end(), std::ostream_iterator<int>(of, " "));
		of << '\n';
	}
	of.close();

}

void output_result(HclustResult const & result, std::string prefix="")
{
	// outputs that can be read into R
	std::ofstream of;
	std::string fname(prefix+"indices");
	of.open(fname.c_str());
	std::copy(result.indices.begin(), result.indices.end(), std::ostream_iterator<int>(of, " "));
	of << '\n';
	of.close();

	fname = prefix+"labels";
	of.open(fname.c_str());
	std::copy(result.ids.begin(), result.ids.end(), std::ostream_iterator<std::string>(of, " "));
	of << '\n';
	of.close();

	fname = prefix+"merges";
	of.open(fname.c_str());
	std::copy(result.merge.begin(), result.merge.end(), std::ostream_iterator<int>(of, " "));
	of << '\n';
	of.close();

	fname = prefix+"heights";
	of.open(fname.c_str());
	std::copy(result.height.begin(), result.height.end(), std::ostream_iterator<t_float>(of, " "));
	of << '\n';
	of.close();

	fname = prefix+"orders";
	of.open(fname.c_str());
	std::copy(result.order.begin(), result.order.end(), std::ostream_iterator<int>(of, " "));
	of << '\n';
	of.close();

//	output_clusters(result,prefix);
}

void run_fastcluster(HclustResult & bs_result, t_float * dist, int method)
{
	std::cout << "Fastcluster..." << std::endl;

	if (method<METHOD_METR_SINGLE || method>METHOD_METR_MEDIAN) {
		std::cerr << "ERROR: unknown clustering method requested" << std::endl;
		exit(EXIT_FAILURE);
	}

	const int N(bs_result.indices.size());
	// 'members': members in each node (for 'clustering in the middle of the tree')
	// here, just setting this whole thing to 1 (all data are for terminal branches)
	auto_array_ptr<t_float> members;
	members.init(N);
	for (int i(0); i<N; ++i) members[i] = 1;
	cluster_result hc(N-1);

	switch (method) {

		case METHOD_METR_SINGLE:
			std::cout << "Method " << method << " (SINGLE)" << std::endl;
			MST_linkage_core(N, dist, hc);
			break;
		case METHOD_METR_COMPLETE:
			std::cout << "Method " << method << " (COMPLETE)" << std::endl;
			NN_chain_core<METHOD_METR_COMPLETE, t_float>(N, dist, NULL, hc);
			break;
		case METHOD_METR_AVERAGE:
			std::cout << "Method " << method << " (AVERAGE)" << std::endl;
			NN_chain_core<METHOD_METR_AVERAGE, t_float>(N, dist, members, hc);
			break;
		case METHOD_METR_WEIGHTED:
			std::cout << "Method " << method << " (WEIGHTED)" << std::endl;
			NN_chain_core<METHOD_METR_WEIGHTED, t_float>(N, dist, NULL, hc);
			break;
		case METHOD_METR_WARD:
			std::cout << "Method " << method << " (WARD)" << std::endl;
			NN_chain_core<METHOD_METR_WARD, t_float>(N, dist, members, hc);
			break;
		case METHOD_METR_CENTROID:
			std::cout << "Method " << method << " (CENTROID)" << std::endl;
			generic_linkage<METHOD_METR_CENTROID, t_float>(N, dist, members, hc);
			break;
		case METHOD_METR_MEDIAN:
			std::cout << "Method " << method << " (MEDIAN)" << std::endl;
			generic_linkage<METHOD_METR_MEDIAN, t_float>(N, dist, NULL, hc);
			break;
		default:
			std::cerr << "ERROR: unknown clustering method requested" << std::endl;
			exit(EXIT_FAILURE);
	}

	std::cout << "Getting clustering results..." << std::endl;
	// here init vector sizes and pass array pointers to fastcluster
	bs_result.merge.resize(2*(N-1));
	bs_result.height.resize(N-1);
	bs_result.order.resize(N);

	if (method==METHOD_METR_CENTROID ||
	    method==METHOD_METR_MEDIAN)
		generate_R_dendrogram<true>(&bs_result.merge[0], &bs_result.height[0], &bs_result.order[0], hc, N);
	else
		generate_R_dendrogram<false>(&bs_result.merge[0], &bs_result.height[0], &bs_result.order[0], hc, N);

	members.free();
}

////////////////////////////////////////////////////////////////////////////////
void usage_error()
{
	std::cerr << "\n"
	 << " -r|--ratios     ratiosfile\n"
	 << " -b|--bootstraps     #              : number of bootstraps\n"
	 << " -m|--method     #              : clustering method (0. single, 1. complete, 2. average, 3. weighted, 4. ward, 5. centroid, 6. median)\n"
	 << "example: [executable] -r ratios.tab\n"
	 << "\n";
	exit(EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

std::srand(time(NULL));
//		std::srand(1);

	std::cout << std::endl;

	std::string ratiosfilename;
	unsigned bootstraps(100);
	int method(0);

	if (argc < 2) usage_error();

	// parse command line arguments
	for (int i(1); i < argc; ++i) {

		std::string arg(argv[i]);

		if (arg == "-r" || arg == "--ratios") {
			if (++i >= argc) usage_error();
			ratiosfilename = argv[i];

		} else if (arg == "-b" || arg == "--bootstraps") {
			if (++i >= argc) usage_error();
			bootstraps = atoi(argv[i]);

		} else if (arg == "-m" || arg == "--method") {
			if (++i >= argc) usage_error();
			method = atoi(argv[i]);

		} else if (arg == "-h" || arg == "--help") {
			usage_error();

		}
	}

	// read in matrix
	Matrix rr;
	Labels ids;
	read_matrix_file(ratiosfilename, rr, ids);
	size_t ntest(std::min(size_t(5),rr.size()));
	size_t nprint(ntest);

//	for(size_t i(0); i<ntest; ++i){ std::cout << ids[i] << " " << rr[i] << std::endl; }

	// test data and distance metric
	std::cout << "Testing distance metric (Pearson distance) for " << ntest << " genes:" << std::endl;
	for(size_t i(0); i<ntest; ++i){
		for(size_t j(i+1); j<ntest; ++j){
			t_float cor(pearson_correlation(rr[i],rr[j],true));
			std::cout << ids[i] << " vs " << ids[j] << "(Pearson correlation): " << cor << std::endl;
		}
	}

	// Parameter N: number of elements
	const t_index N(rr.size());

	// Parameter NN: number of non-redundant, non-self comparisons
	const std::ptrdiff_t NN = static_cast<std::ptrdiff_t>((N)*(N-1)/2);
	std::cout << "N is " << N << " and NN is " << NN << std::endl;
	std::cout << "Dissimilarity matrix..." << std::endl;
	auto_array_ptr<t_float> dist;
	dist.init(NN);
	std::ptrdiff_t p(0);
	for(size_t r1(0), end(rr.size()); r1<(end-1); ++r1){
		for(size_t r2(r1+1); r2<end; ++r2){
			dist[p] = 1.0 - pearson_correlation(rr[r1],rr[r2],true);
			// fastcluster takes squared distances for some metrics (i.e. Euclidean)
//			dist[p] = pow(1.0 - pearson_correlation(rr[r1],rr[r2],true),2);
			++p;
		}
	}

	// here: do an initial (non-bootstrapped) hclust
	HclustResult result;
	// fill indices
	for(size_t i(0); i<ids.size(); ++i){
		result.indices.push_back(i);
		result.ids.push_back(ids[i]);
	}

	run_fastcluster(result, dist, method);

	// bootstrap iterations
	HclustResults bootstrap_results;

	for(unsigned bs(0); bs<bootstraps; ++bs){
		std::cout << "Bootstrap " << bs << ":" << std::endl;

		HclustResult bs_result;

		// resample indices from [0,n) with replacement
    for(t_index i(0); i<N; ++i) bs_result.indices.push_back(rand() % N);
		std::cout << "Resampled indices: ";
		std::copy(bs_result.indices.begin(), bs_result.indices.begin()+nprint, std::ostream_iterator<int>(std::cout, " "));
		std::cout << std::endl;

		// store resampled id order
		for(size_t i(0); i<bs_result.indices.size(); ++i){
			bs_result.ids.push_back( ids[bs_result.indices[i]] );
		}

		// resample from the distance matrix, keeping track of resampled indices
		// this should be faster than recomputing all of the same distances again
		auto_array_ptr<t_float> dist_resampled;
		dist_resampled.init(NN);
		std::cout << "Resampling distances..." << std::endl;
		resample_distances(dist, dist_resampled, bs_result.indices);

		run_fastcluster(bs_result, dist_resampled, method);
		bootstrap_results.push_back(bs_result);

		dist_resampled.free();
	}
	dist.free();     // Free the memory now

	output_result(result,"hc.");

	output_bootstrap_results(bootstrap_results);

	// to do: height-based agglomeration to yield k clusters,
	// empirical distributions of pairwise cluster co-memberships for many choices of k

	return 0;
}

