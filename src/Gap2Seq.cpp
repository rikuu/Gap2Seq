/*****************************************************************************
 *  Gap2Seq
 *  Copyright (C) Leena Salmela, Kristoffer Sahlin, Veli Mäkinen,
 *  Alexandru Tomescu 2015
 *
 *  Contact: leena.salmela@cs.helsinki.fi
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/graphviz.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

#include <gatb/gatb_core.hpp>

#include <Gap2Seq.hpp>

//#define DEBUG

// Maximum number of paths for counting
#define MAX_PATHS (INT_MAX/2-1)

// Deafault parameter values
#define DEFAULT_K "31"
#define DEFAULT_SOLID "2"
#define DEFAULT_DIST_ERR "500"
#define DEFAULT_FUZ "10"
#define DEFAULT_MAX_MEMORY "20"  // Maximum memory usage of DP table per thread


/********************************************************************************/

// Constants for command line parameters
static const char* STR_KMER_LEN = "-k";
static const char* STR_SOLID_THRESHOLD = "-solid";
static const char* STR_READS = "-reads";
static const char* STR_SCAFFOLDS = "-scaffolds";
static const char* STR_FILLED_SCAFFOLDS = "-filled";
static const char* STR_DIST_ERROR = "-dist-error";
static const char* STR_FUZ = "-fuz";
static const char* STR_MAX_MEM = "-max-mem";
static const char* STR_SKIP_CONFIDENT = "-all-upper";
static const char* STR_UNIQUE = "-unique";

// For filling a single gap
static const char* STR_LEFT = "-left";
static const char* STR_RIGHT = "-right";
static const char* STR_LENGTH = "-length";

/*********************************************************************
Constructor for the tool.
*********************************************************************/
Gap2Seq::Gap2Seq() : Tool("Gap2Seq")
{
    // We add some custom arguments for command line interface
    getParser()->push_front (new OptionOneParam (STR_KMER_LEN, "kmer length",  false, DEFAULT_K));
    getParser()->push_front (new OptionOneParam (STR_SOLID_THRESHOLD, "Threshold for solid k-mers",  false, DEFAULT_SOLID));
    getParser()->push_front (new OptionOneParam (STR_READS, "FASTA/Q files of reads. For several files use a comma separated list.",  true));
    getParser()->push_front (new OptionOneParam (STR_SCAFFOLDS, "FASTA/Q file of scaffolds",  false, ""));
    getParser()->push_front (new OptionOneParam (STR_FILLED_SCAFFOLDS, "FASTA file of filled scaffolds",  false, ""));
    getParser()->push_front (new OptionOneParam (STR_DIST_ERROR, "Maximum error in gap estimates",  false, DEFAULT_DIST_ERR));
    getParser()->push_front (new OptionOneParam (STR_FUZ, "Number of nucleotides to ignore on gap fringes",  false, DEFAULT_FUZ));
    getParser()->push_front (new OptionOneParam (STR_MAX_MEM, "Maximum memory usage of DP table computation in gigabytes (excluding DBG)",  false, DEFAULT_MAX_MEMORY));
    getParser()->push_front (new OptionNoParam (STR_SKIP_CONFIDENT, "If specified, all filled bases are in upper case.",  false, false));
    getParser()->push_front (new OptionNoParam (STR_UNIQUE, "If specified, only gaps with a unique path of best length are filled.",  false, false));

    // TODO: Integrate all gap cutting and remove these
    getParser()->push_front (new OptionOneParam (STR_LEFT, "Left flank of a single gap", false, ""));
    getParser()->push_front (new OptionOneParam (STR_RIGHT, "Right flank of a single gap", false, ""));
    getParser()->push_front (new OptionOneParam (STR_LENGTH, "Length of a single gap", false, ""));
}

// Synchronizer for access to global variables
#ifndef SINGLE_THREAD
ISynchronizer *global_lock;
#endif

// Array of memory usage accessed by thread id
std::unordered_map<long long, long long> memuse;

// Print statistics for gap filling
void print_statistics(int filledStart, int gapStart, int gapEnd, int paths, char *buf, int fuz, int k, int left_fuz, int right_fuz, bool skip_confident, bool unique_paths, subgraph_stats substats, int gap, std::string comment) {
  #ifndef SINGLE_THREAD
    {
      LocalSynchronizer local(global_lock);
  #endif

      if (paths > 0 && (!unique_paths || paths == 1)) {
        int filledLen = (int)strlen(&buf[fuz-left_fuz])-k;

        int lower=0, upper=0;
        if (!skip_confident) {
          for (int j = 0; j < filledLen; j++) {
            if (isupper(buf[j+fuz-left_fuz])) {
              upper++;
            } else {
              lower++;
            }
          }

          std::cout << "Scaffold: " << comment << " GapStart: " << gapStart << " GapEnd: " << gapEnd <<
            " GapLength: " << gap << " PathsFound: " << paths << " FilledStart: " << filledStart <<
            " FilledEnd: " << filledStart + filledLen << " FilledGapLength: " <<  filledLen <<
            " LeftFuz: " << left_fuz << " RightFuz: " << right_fuz << " ConfidentBases: " << upper << " TotalBases: " << (upper+lower) << std::endl;
          std::cout << "SubgraphStats: Vertices: " << substats.vertices << " Edges: " <<
            substats.edges << " NontrivialStrongComponents: " << substats.nontrivial_components <<
            " SizeNontrivialStrongComponents: " << substats.size_nontrivial_components << " VerticesFinal: " <<
            substats.vertices_final << " EdgesFinal: " << substats.edges_final << std::endl;
        } else {
          std::cout << "Scaffold: " << comment << " GapStart: " << gapStart << " GapEnd: " << gapEnd <<
            " GapLength: " << gap << " PathsFound: " << paths << " FilledStart: " << filledStart <<
            " FilledEnd: " << filledStart + filledLen << " FilledGapLength: " <<  filledLen <<
            " LeftFuz: " << left_fuz << " RightFuz: " << right_fuz << std::endl;
        }
      } else {
        if (paths == -1) {
          std::cout << "Scaffold: " << comment << " GapStart: " << gapStart << " GapEnd: " << gapEnd <<
            " GapLength: " << gap << " PathsFound: 0 FilledStart: 0 FilledEnd: 0 FilledGapLength: 0 LeftFuz: " << left_fuz <<
            " RightFuz: " << right_fuz << " Memory limit exceeded" << std::endl;
        } else {
          std::cout << "Scaffold: " << comment << " GapStart: " << gapStart << " GapEnd: " << gapEnd <<
          " GapLength: " << gap << " PathsFound: " << paths << " FilledStart: 0 FilledEnd: 0 FilledGapLength: 0 LeftFuz: " << left_fuz <<
          " RightFuz: " << right_fuz << std::endl;
        }
      }
  #ifndef SINGLE_THREAD
    }
  #endif
}

/*********************************************************************
Main method for gap filling
*********************************************************************/
void Gap2Seq::execute ()
{
  // Initialize random number generation for choosing random paths
  srand(time(NULL));

  // Get the command line arguments
  int k = getInput()->getInt(STR_KMER_LEN);
  int solid = getInput()->getInt(STR_SOLID_THRESHOLD);
  std::string reads = getInput()->getStr(STR_READS);
  int d_err = getInput()->getInt(STR_DIST_ERROR);
  int fuz = getInput()->getInt(STR_FUZ);
  long long max_mem = (long long)(getInput()->getDouble(STR_MAX_MEM) * 1024*1024*1024);
  std::string readsGraph = reads + ".h5";
  bool skip_confident = (getInput()->get(STR_SKIP_CONFIDENT) != 0);
  bool unique_paths = (getInput()->get(STR_UNIQUE) != 0);

  std::cout << "k-mer size: " << k << std::endl;
  std::cout << "Solidity threshold: " << solid << std::endl;
  std::cout << "Reads file: " << reads << std::endl;
  // std::cout << "Scaffolds file: " << scaffolds << std::endl;
  // std::cout << "Filled scaffolds file: " << filled_scaffolds << std::endl;
  std::cout << "Distance error: " << d_err << std::endl;
  std::cout << "Fuz: " << fuz << std::endl;
  std::cout << "Max memory: " << max_mem << std::endl;
  std::cout << "Skip confident: " << skip_confident << std::endl;
  std::cout << "Unique: " << unique_paths << std::endl;

  // Create de bruijn graph
  Graph graph;
  if (is_readable(readsGraph)) {
    std::cout << "Loading from " << readsGraph << std::endl;
    graph = Graph::load(reads);
  } else {
    // Tokenize reads file (list of files separated by ,)
    std::vector<std::string> files;
    std::stringstream ss(reads.c_str());
    std::string file;
    while (std::getline(ss, file, ',')) {
      /*if (!is_readable(file)) {
        std::cerr << "Cannot access: " << file << std::endl;
        exit(EXIT_FAILURE);
      }*/

      files.push_back(file);
    }

    try {
      BankAlbum *b = new BankAlbum(files);
      graph = Graph::create(b, (char const *)"-kmer-size %d -abundance-min %d -debloom original -debloom-impl basic", k, solid);
    } catch (Exception& e) {
      std::cout << "DBG building failed: " << e.getMessage() << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  std::cout << graph.getInfo() << std::endl;

  if (getParser()->saw(STR_LEFT) && getParser()->saw(STR_RIGHT) && getParser()->saw(STR_LENGTH)) {
    std::string kmer_left = getInput()->getStr(STR_LEFT);
    std::string kmer_right = getInput()->getStr(STR_RIGHT);
    int length = getInput()->getInt(STR_LENGTH);

    // TODO: Use non-static flanks, i.e. flanks are >= k upto k+fuz
    assert(kmer_left.length() == k+fuz);
    assert(right_flank.length() == k+fuz);

    // Buffer to hold the sequence to fill the gap
    char *buf = new char[length+k+d_err+2*fuz+1+2];

    // How many k-mers on the edge of the gap are ignored by the paths
    int left_fuz = 0, right_fuz = 0;

    // Struct for subgraph statistics
    struct subgraph_stats substats;

    // Number of paths found
    int s = 0;
    s = fill_gap(graph, kmer_left, kmer_right, length+k, d_err, fuz, k, max_mem, &left_fuz, &right_fuz, buf, skip_confident, &substats);
    print_statistics(0, kmer_left.length(), kmer_left.length(), s, buf, fuz, k, left_fuz, right_fuz, skip_confident, unique_paths, substats, "", "");

    // At least one path found
	  if (s > 0 && (!unique_paths || s == 1)) {
      // Remove ignored edges from flanks, insert filled sequence and output
      // TODO: right flank is included?
      std::cout << kmer_left.substr(0, k+fuz-left_fuz) << &buf[fuz-left_fuz] << std::endl;
    }

    return;
  }

  // Open the input scaffold file
  std::string scaffolds = getInput()->getStr(STR_SCAFFOLDS);
  BankFasta scaffoldBank(scaffolds);
  BankFasta::Iterator itSeq(scaffoldBank);

  // Open the output scaffold file
  std::string filled_scaffolds = getInput()->getStr(STR_FILLED_SCAFFOLDS);
  BankFasta output(filled_scaffolds);

  // Count gaps and filled gaps
  int gapcount = 0;
  int filledgapcount = 0;

#ifndef SINGLE_THREAD
  // Synchronization for output (both file and stdout)
  global_lock = System::thread().newSynchronizer();
  IDispatcher *d = getDispatcher();
  // Default group size is too big (one thread will do all the work...)
  d->setGroupSize(1);

  max_mem = max_mem / d->getExecutionUnitsNumber();
#endif

  std::cout << "Max mem: " << max_mem << std::endl;

  // Parallel code starts here
#ifdef SINGLE_THREAD
  for(itSeq.first(); !itSeq.isDone();itSeq.next()) {
    // Original sequence
    std::string seq = itSeq->toString();
    std::string comment = itSeq->getComment();
#else
  // IDispatcher::Status status = d->iterate(&itSeq, [&] (const Sequence& sequence) {
  Range<int>::Iterator intIt(1,d->getExecutionUnitsNumber());
  itSeq.first();
  IDispatcher::Status status = d->iterate(intIt, [&] (int id) {
    bool more = true;

    while (more) {
    std::string seq = "";
    std::string comment = "";
    {
      LocalSynchronizer local(global_lock);
      if (!itSeq.isDone()) {
	seq = itSeq->toString();
	comment = itSeq->getComment();
	itSeq.next();
	more = true;
      } else {
	more = false;
      }
    }
    if (!more)
      break;
#endif

    // Gap filled sequence
    std::string filledSeq = "";

    // Position where previous gap ended. The filled sequence is
    // complete up to this position (in the original sequence)
    int prevGapEnd = 0;

    // Index over sequence positions
    size_t i = 0;

    // Iterate over the sequence and find gaps
    while (i < seq.size()) {
      // A gap starts here
      if (seq[i] == 'N' || seq[i] == 'n') {
#ifndef SINGLE_THREAD
	{
	  LocalSynchronizer local(global_lock);
#endif
	  gapcount++;
#ifndef SINGLE_THREAD
	}
#endif
	// The first possible starting position for left k-mer
	int kmer_start = i-k-fuz;
	// The length of gap
	int gap = 0;
	// Measure the gap length
	while (i < seq.size() && (seq[i] == 'N' || seq[i] == 'n')) {
	  i++;
	  gap++;
	}

	// Check that the end k-mer is complete, i.e. no N's and enough sequence
	bool ok = i+k+fuz <= seq.size() ? true : false;
	for (int j = 0; j < k+fuz && ok; j++) {
	  if (seq[i+j] == 'N' || seq[i+j] == 'n')
	    ok = false;
	}

	// Check that the start k-mer is complete
	if (kmer_start >= prevGapEnd && ok) {
    // Buffer to hold the sequence to fill the gap
    char *buf = new char[gap+k+d_err+2*fuz+1+2];

    // How many k-mers on the edge of the gap are ignored by the paths
    int left_fuz=0, right_fuz=0;

    // Struct for subgraph statistics
    struct subgraph_stats substats;

    // Number of paths found
    int s = 0;
    s = fill_gap(graph, seq.substr(kmer_start,k+fuz), seq.substr(i,k+fuz), gap+k, d_err, fuz, k, max_mem, &left_fuz, &right_fuz, buf, skip_confident, &substats);

    int filledStart = filledSeq.length() + kmer_start+k+fuz-left_fuz-prevGapEnd;
    int gapStart = kmer_start+k+fuz;
    print_statistics(filledStart, gapStart, i, s, buf, fuz, k, left_fuz, right_fuz, skip_confident, unique_paths, substats, gap, comment);

    // At least one path was found
    if (s > 0 && (!unique_paths || s == 1)) {
      filledgapcount++;

  #ifdef DEBUG
      std::cout << "Fill: " << &buf[fuz-left_fuz] << std::endl;
  #endif

      // Add seq from end of previous gap to gap start and the gap fill sequence
      filledSeq = seq.substr(prevGapEnd, kmer_start+k+fuz-left_fuz-prevGapEnd) + (string)&buf[fuz-left_fuz];

      // Remove the right kmer
      filledSeq = filledSeq.substr(0, filledSeq.length()-k);

      // adjust gap end index
      i += right_fuz;
    } else {
      // No paths found, keep the gap
      // Add seq from end of previous gap to end of this gap
      filledSeq = filledSeq + seq.substr(prevGapEnd, kmer_start+k+fuz+gap-prevGapEnd);
    }
    delete [] buf;
	} else {
	  // Either left or right k-mer is not complete
	  // Add seq from end of previous gap to end of this gap
	  filledSeq = filledSeq + seq.substr(prevGapEnd, kmer_start+k+fuz+gap-prevGapEnd);
	}
	// This gap has been processed
	prevGapEnd = i;
      } else {
	// No gap at position i
	i++;
      }
    }
    // Add the end of the seq
    filledSeq = filledSeq + seq.substr(prevGapEnd, seq.length()-prevGapEnd);

    // Write the filled scaffold to file
    Sequence s((char *)filledSeq.c_str());
    s._comment = comment;
#ifndef SINGLE_THREAD
    {
      LocalSynchronizer local(global_lock);
#endif
      output.insert(s);
      output.flush();
#ifndef SINGLE_THREAD
    }
#endif
#ifdef DEBUG
    std::cout << filledSeq << std::endl;
#endif
#ifdef SINGLE_THREAD
    }
#else
    }
    });
#endif
  output.flush();

  std::cout << "Filled " << filledgapcount << " gaps out of " << gapcount << std::endl;
}

// Custom memory allocator to be used with stl containers
template <typename T>
class count_allocator: public std::allocator<T> {
public:
  typedef size_t size_type;
  typedef T* pointer;
  typedef const T* const_pointer;

  template<typename _Tp1>
  struct rebind {
    typedef count_allocator<_Tp1> other;
  };

  pointer allocate(size_type n, const void *hint=0) {
#ifdef DEBUG
    fprintf(stderr, "Alloc %lu bytes.\n", n*sizeof(T));
#endif
    long long id = System::thread().getThreadSelf();
#ifndef SINGLE_THREAD
    {
      LocalSynchronizer local(global_lock);
#endif
      memuse[id] += n*sizeof(T);
#ifndef SINGLE_THREAD
    }
#endif
    return std::allocator<T>::allocate(n, hint);
  }

  void deallocate(pointer p, size_type n) {
#ifdef DEBUG
    fprintf(stderr, "Dealloc %lu bytes (%p).\n", n*sizeof(T), p);
#endif
    long long id = System::thread().getThreadSelf();
#ifndef SINGLE_THREAD
    {
      LocalSynchronizer local(global_lock);
#endif
      memuse[id] -= n*sizeof(T);
#ifndef SINGLE_THREAD
    }
#endif
    return std::allocator<T>::deallocate(p, n);
  }

  count_allocator() noexcept {
#ifdef DEBUG
    fprintf(stderr, "Hello allocator!\n");
#endif
  }

  count_allocator(const count_allocator &a) noexcept {}

  template <typename U>
  count_allocator(const count_allocator<U> &a) noexcept {}

  ~count_allocator() noexcept {}
};



// Data structure to hold the dp row for one k-mer (both forward and reverse)
class map_element {
private:
  // dp row, only nonzero entries stored
  std::vector<std::pair<int, int>, count_allocator<std::pair<int, int> > > s_f;
  std::vector<std::pair<int, int>, count_allocator<std::pair<int, int> > > s_r;

#ifdef STATS
  int nonzeros;
#endif

public:
  // Construct new dp row
  map_element(int dd) {
#ifdef STATS
    nonzeros=0;
#endif
  }

  // Construct new dp row
  map_element() {
#ifdef STATS
    nonzeros=0;
#endif
  }

  // Get the i'th value of the dp row for forward k-mer
  // Most accesses are to the last (or second last) item so check those first
  // If accessing other items, use binary search
  int getf(int i) {
    if (s_f.size() == 0)
      return 0;

    if (s_f[s_f.size()-1].first == i)
      return s_f[s_f.size()-1].second;

    if (s_f.size() > 1 && s_f[s_f.size()-2].first == i)
      return s_f[s_f.size()-2].second;

#ifdef DEBUG
    for(int j = 0; j < s_f.size(); j++) {
      std::cout << "(" << s_f[j].first << "," << s_f[j].second << ") ";
    }
    std::cout << std::endl;

    std::cout << i << std::endl;
#endif

    // Binary search
    int lo = 0;
    int hi = s_f.size()-1;
    int mid = (lo + hi)/2;
    while(lo <= hi) {
#ifdef DEBUG
      std::cout << mid << " ";
#endif
      if (s_f[mid].first < i) {
	lo = mid+1;
      } else if (s_f[mid].first > i) {
	hi = mid-1;
      } else {
	return s_f[mid].second;
      }
      mid = (lo+hi)/2;
    }
#ifdef DEBUG
    std::cout << std::endl;
#endif
    return 0;
  }

  // Get the i'th value of the dp row for reverse k-mer
  // Most accesses are to the last (or second last) item so check those first
  // If accessing other items, use binary search
  int getr(int i) {
    if (s_r.size() == 0)
      return 0;

    if (s_r[s_r.size()-1].first == i)
      return s_r[s_r.size()-1].second;

    if (s_r.size() > 1 && s_r[s_r.size()-2].first == i)
      return s_r[s_r.size()-2].second;

#ifdef DEBUG
    for(int j = 0; j < s_r.size(); j++) {
      std::cout << "(" << s_r[j].first << "," << s_r[j].second << ") ";
    }
    std::cout << std::endl;

    std::cout << i << std::endl;
#endif

    // Binary search
    int lo = 0;
    int hi = s_r.size()-1;
    int mid = (lo + hi)/2;
    while(lo <= hi) {
#ifdef DEBUG
      std::cout << mid << " ";
#endif
      if (s_r[mid].first < i) {
	lo = mid+1;
      } else if (s_r[mid].first > i) {
	hi = mid-1;
      } else {
	return s_r[mid].second;
      }
      mid = (lo+hi)/2;
    }
#ifdef DEBUG
    std::cout << std::endl;
#endif
    return 0;
  }

  // Supports only updating the last value and adding a new value
  // (with larger index i than previous entries)
  void setf(int i, int v) {
    if (s_f.size() == 0) {
      s_f.push_back(std::make_pair(i,v));
#ifdef STATS
      nonzeros++;
#endif
      return;
    }

    if (s_f[s_f.size()-1].first == i) {
      s_f[s_f.size()-1].second = v;
    } else {
      if (s_f[s_f.size()-1].first > i) {
	// Just checking for misuse
	std::cout << "Setf will fail: " << i << " < " << s_f[s_f.size()-1].first << std::endl;
      }
      s_f.push_back(std::make_pair(i,v));
#ifdef STATS
      nonzeros++;
#endif
    }
  }

  // Supports only updating the last value and adding a new value
  // (with larger index i than previous entries)
  void setr(int i, int v) {
    if (s_r.size() == 0) {
      s_r.push_back(std::make_pair(i,v));
#ifdef STATS
      nonzeros++;
#endif
      return;
    }

    if (s_r[s_r.size()-1].first == i) {
      s_r[s_r.size()-1].second = v;
    } else {
      if (s_r[s_r.size()-1].first > i) {
	// Just checking for misuse
	std::cout << "Setr will fail: " << i << " < " << s_r[s_r.size()-1].first << std::endl;
      }
      s_r.push_back(std::make_pair(i,v));
#ifdef STATS
      nonzeros++;
#endif
    }
  }

  std::vector<std::pair<int, int>, count_allocator<std::pair<int, int> > >::reverse_iterator getf_rbegin() {
    return s_f.rbegin();
  }
  std::vector<std::pair<int, int>, count_allocator<std::pair<int, int> > >::reverse_iterator getf_rend() {
    return s_f.rend();
  }
  std::vector<std::pair<int, int>, count_allocator<std::pair<int, int> > >::reverse_iterator getr_rbegin() {
    return s_r.rbegin();
  }
  std::vector<std::pair<int, int>, count_allocator<std::pair<int, int> > >::reverse_iterator getr_rend() {
    return s_r.rend();
  }



#ifdef STATS
  int get_nz() {
    return nonzeros;
  }
#endif
};


// Data structure to hold the dp row for one k-mer (both forward and reverse)
// No counting, just path existance
class map_element2 {
private:
  // dp row, only nonzero entries stored
  std::vector< int, count_allocator<int> > s_f;
  std::vector< int, count_allocator<int> > s_r;

#ifdef STATS
  int nonzeros;
#endif

public:
  // Construct new dp row
  map_element2(int dd) {
#ifdef STATS
    nonzeros=0;
#endif
  }

  // Construct new dp row
  map_element2() {
#ifdef STATS
    nonzeros=0;
#endif
  }

  // Get the i'th value of the dp row for forward k-mer
  // Most accesses are to the last (or second last) item so check those first
  // If accessing other items, use binary search
  bool getf(int i) {
    if (s_f.size() == 0)
      return false;

    if (s_f[s_f.size()-1] == i)
      return true;

    if (s_f.size() > 1 && s_f[s_f.size()-2] == i)
      return true;

#ifdef DEBUG
    for(int j = 0; j < s_f.size(); j++) {
      std::cout << s_f[j] << " ";
    }
    std::cout << std::endl;

    std::cout << i << std::endl;
#endif

    // Binary search
    int lo = 0;
    int hi = s_f.size()-1;
    int mid = (lo + hi)/2;
    while(lo <= hi) {
#ifdef DEBUG
      std::cout << mid << " ";
#endif
      if (s_f[mid] < i) {
	lo = mid+1;
      } else if (s_f[mid] > i) {
	hi = mid-1;
      } else {
	return true;
      }
      mid = (lo+hi)/2;
    }
#ifdef DEBUG
    std::cout << std::endl;
#endif
    return false;
  }

  // Get the i'th value of the dp row for reverse k-mer
  // Most accesses are to the last (or second last) item so check those first
  // If accessing other items, use binary search
  bool getr(int i) {
    if (s_r.size() == 0)
      return false;

    if (s_r[s_r.size()-1] == i)
      return true;

    if (s_r.size() > 1 && s_r[s_r.size()-2] == i)
      return true;

#ifdef DEBUG
    for(int j = 0; j < s_r.size(); j++) {
      std::cout << s_r[j] << ") ";
    }
    std::cout << std::endl;

    std::cout << i << std::endl;
#endif

    // Binary search
    int lo = 0;
    int hi = s_r.size()-1;
    int mid = (lo + hi)/2;
    while(lo <= hi) {
#ifdef DEBUG
      std::cout << mid << " ";
#endif
      if (s_r[mid] < i) {
	lo = mid+1;
      } else if (s_r[mid] > i) {
	hi = mid-1;
      } else {
	return true;
      }
      mid = (lo+hi)/2;
    }
#ifdef DEBUG
    std::cout << std::endl;
#endif
    return false;
  }

  // Supports only updating the last value and adding a new value
  // (with larger index i than previous entries)
  void setf(int i) {
    if (s_f.size() == 0) {
      s_f.push_back(i);
#ifdef STATS
      nonzeros++;
#endif
      return;
    }

    if (s_f[s_f.size()-1] == i) {
      return;
    } else {
      if (s_f[s_f.size()-1] > i) {
	// Just checking for misuse
	std::cout << "Setf will fail: " << i << " < " << s_f[s_f.size()-1] << std::endl;
      }
      s_f.push_back(i);
#ifdef STATS
      nonzeros++;
#endif
    }
  }

  // Supports only updating the last value and adding a new value
  // (with larger index i than previous entries)
  void setr(int i) {
    if (s_r.size() == 0) {
      s_r.push_back(i);
#ifdef STATS
      nonzeros++;
#endif
      return;
    }

    if (s_r[s_r.size()-1] == i) {
      return;
    } else {
      if (s_r[s_r.size()-1] > i) {
	// Just checking for misuse
	std::cout << "Setr will fail: " << i << " < " << s_r[s_r.size()-1] << std::endl;
      }
      s_r.push_back(i);
#ifdef STATS
      nonzeros++;
#endif
    }
  }

  std::vector<int, count_allocator<int> >::reverse_iterator getf_rbegin() {
    return s_f.rbegin();
  }
  std::vector<int, count_allocator<int> >::reverse_iterator getf_rend() {
    return s_f.rend();
  }
  std::vector<int, count_allocator<int> >::reverse_iterator getr_rbegin() {
    return s_r.rbegin();
  }
  std::vector<int, count_allocator<int> >::reverse_iterator getr_rend() {
    return s_r.rend();
  }



#ifdef STATS
  int get_nz() {
    return nonzeros;
  }
#endif
};


// Providing hashes for Nodes
struct node_hash {
  size_t operator()(const Node & node) const{
    return oahash(node.kmer);
  }
};

// Fill a gap. Returns the number of paths found and -1 if running out of memory.
// kmer_left: sequence of length k+fuz for extracting k-mers from the left end of the gap
// kmer_right: sequence of length k+fuz for extracting k-mers from the right end of the gap
// gap_len: estimated length of gap
// gap_err: upper bound on the error in gap length
// fuz: maximum number of gap bordering nucleotides that can be ignored
// k: k-mer size
// The following should be non-null if a (ad hoc random) filling sequence should be recovered by traceback.
// Otherwise they should all be NULL.
// left_fuz: The number of actual ignored nucleotides on the left end of the gap
// right_fuz: The number of actual ignored nucleotides on the right end of the gap
// fill: The sequence to fill the gap. Start offset of actual sequence in fill is fuz - left_fuz,
//       the fill seq always includes the right kmer
int Gap2Seq::fill_gap(Graph graph, std::string kmer_left, std::string kmer_right, int gap_len, int gap_err, int fuz, int k, long long max_mem, int *left_fuz, int *right_fuz, char *fill, bool skip_confident, struct subgraph_stats *substats) {

  // For simplicity, force gap_len and gap_err to be even
  if (gap_len % 2 == 1) {
    gap_len++;
  }
  if (gap_err % 2 == 1) {
    gap_err++;
  }

  long long id = System::thread().getThreadSelf();
#ifndef SINGLE_THREAD
  {
      LocalSynchronizer local(global_lock);
#endif
      memuse[id] = 0;
#ifndef SINGLE_THREAD
  }
#endif

  // Border sets for breadth first search
  std::unordered_set<Node, node_hash, equal_to<Node>, count_allocator<Node> > border;
  std::unordered_set<Node, node_hash, equal_to<Node>, count_allocator<Node> > nextBorder;

  std::unordered_map<Node, map_element2 *, node_hash, equal_to<Node>, count_allocator< std::pair <const Node, map_element2 *> > > reachableSetRight;

#ifdef DEBUG
  std::cout << "Search from right end" << std::endl;
  std::cout << "Kmer " << kmer_right << std::endl;
#endif

  // Rightmost starting k-mer
  std::string kmer = kmer_right.substr(kmer_right.length()-k,k);
  Node node = graph.buildNode(Data((char *)kmer.c_str()));

#ifdef DEBUG
  if (!graph.contains(node)) {
    std::cout << "Kmer " << kmer << " cannot be found in the graph!" << std::endl;
  }
#endif

  // Put the rightmost k-mer into border and reachable set
  int currentD = 0;

  if (graph.contains(node)) {
    border.insert(node);
    map_element2 *me = NULL;

    if (reachableSetRight.find(node) != reachableSetRight.end()) {
      me = reachableSetRight[node];
    }

    if (me == NULL) {
      me = new map_element2(gap_len/2+gap_err/2+fuz);
      reachableSetRight[node] = me;
    }

    // Initialize the dp row
    if (node.strand == STRAND_FORWARD) {
      me->setf(currentD);
    } else {
      me->setr(currentD);
    }
  }

  long long mymemuse;
#ifndef SINGLE_THREAD
  {
      LocalSynchronizer local(global_lock);
#endif
      mymemuse = memuse[id];
#ifndef SINGLE_THREAD
  }
#endif

  // BFS loop from the right until depth d/2+fuz is reached
  while(currentD <= gap_len/2+gap_err/2+fuz+1 && mymemuse < max_mem) {
    currentD++;

#ifdef DEBUG
    std::cout << "currentD: " << currentD << " Limit: " << gap_len/2+gap_err/2+fuz << std::endl;
    std::cout << "Memuse: " << mymemuse << std::endl;
#endif

    // Iterate over the nodes at current depth
#ifdef DEBUG
    std::cout << "Border size: " << border.size() << std::endl;
#endif

    for (std::unordered_set<Node, node_hash, equal_to<Node>, count_allocator<Node> >::iterator it = border.begin(); it != border.end(); ++ it) {
#ifndef SINGLE_THREAD
      {
	LocalSynchronizer local(global_lock);
#endif
	mymemuse = memuse[id];
#ifndef SINGLE_THREAD
      }
#endif

      if (mymemuse >= max_mem)
	break;

      Node n = (Node) *it;

#ifdef DEBUG
      std::cout << graph.toString(n) << (n.strand == STRAND_FORWARD ? " F": " R") <<  std::endl;
#endif

      // Get the neightbors of the node and update their dp rows
      Graph::Vector<Node> neighbors = graph.predecessors(n);
      map_element2 *node_me = reachableSetRight[n];
      // Number of paths to the parent node
      bool num_paths = n.strand == STRAND_FORWARD ? node_me->getf(currentD-1) : node_me->getr(currentD-1);

      for (size_t i = 0; i < neighbors.size(); i++) {
	map_element2 *me = NULL;

	if (reachableSetRight.find(neighbors[i]) != reachableSetRight.end()) {
	  me = reachableSetRight[neighbors[i]];
	}

	if (me == NULL) {
	  me = new map_element2(gap_len/2+gap_err/2+fuz);
	  reachableSetRight[neighbors[i]] = me;
	}
	if (neighbors[i].strand == STRAND_FORWARD) {
	  me->setf(currentD);
	} else {
	  me->setr(currentD);
	}

#ifdef DEBUG
	std::cout << "Neighbor: " << graph.toString(neighbors[i]) << (neighbors[i].strand == STRAND_FORWARD ? " F": " R") <<  std::endl;
#endif
	// Put the neighbor to the next border
	nextBorder.insert(neighbors[i]);
      }
    }

    // Swap the borders
    border.clear();
    border.swap(nextBorder);

    // Add next starting k-mer if one still exists
    if (currentD <= fuz) {
      std::string kmer = kmer_right.substr(kmer_right.length()-k-currentD, k);
      Node node = graph.buildNode(Data((char *)kmer.c_str()));
#ifdef DEBUG
      if (!graph.contains(node)) {
	std::cout << "Kmer " << kmer_right << " cannot be found in the graph!" << std::endl;
      }
#endif

#ifdef DEBUG
      std::cout << "Kmer: " << kmer << std::endl;
#endif

      if (graph.contains(node)) {
	border.insert(node);
	map_element2 *me = NULL;

	if (reachableSetRight.find(node) != reachableSetRight.end()) {
	  me = reachableSetRight[node];
	}

	if (me == NULL) {
	  me = new map_element2(gap_len+gap_err+2*fuz);
	  reachableSetRight[node] = me;
	}

	// Initialize the dp row. Note that here we initialize to 1
	// regardless of whether the node has already been reached. This
	// will prevent counting paths going via all the starting k-mers
	// several times
	if (node.strand == STRAND_FORWARD) {
	  me->setf(currentD);
	} else {
	  me->setr(currentD);
	}
      }
    }

#ifndef SINGLE_THREAD
    {
      LocalSynchronizer local(global_lock);
#endif
      mymemuse = memuse[id];
#ifndef SINGLE_THREAD
    }
#endif
  }

  // Search from the left
  border.clear();
  nextBorder.clear();

  // The number of paths
  int count = 0;

  // The set of nodes reached so far
  std::unordered_map<Node, map_element *, node_hash, equal_to<Node>, count_allocator< std::pair <const Node, map_element *> > > reachableSetLeft;

#ifdef DEBUG
  std::cout << "Search from left end" << std::endl;
  std::cout << "Kmer " << kmer_left << std::endl;
#endif

  // Leftmost starting k-mer
  kmer = kmer_left.substr(0,k);
  node = graph.buildNode(Data((char *)kmer.c_str()));

#ifdef DEBUG
  if (!graph.contains(node)) {
    std::cout << "Kmer " << kmer << " cannot be found in the graph!" << std::endl;
  }
#endif

  // Put the leftmost k-mer into border and reachable set
  currentD = 0;

  if (graph.contains(node)) {
    border.insert(node);
    map_element *me = NULL;

    if (reachableSetLeft.find(node) != reachableSetLeft.end()) {
      me = reachableSetLeft[node];
    }

    if (me == NULL) {
      me = new map_element(gap_len+gap_err+2*fuz);
      reachableSetLeft[node] = me;
    }

    // Initialize the dp row
    if (node.strand == STRAND_FORWARD) {
      me->setf(currentD, 1);
    } else {
      me->setr(currentD, 1);
    }
  }

  currentD++;

#ifndef SINGLE_THREAD
  {
      LocalSynchronizer local(global_lock);
#endif
      mymemuse = memuse[id];
#ifndef SINGLE_THREAD
  }
#endif

  // BFS loop until depth d/2+fuz is reached
  while(currentD <= gap_len+gap_err+2*fuz && mymemuse < max_mem) {
#ifdef DEBUG
    std::cout << "currentD: " << currentD << " Limit: " << gap_len/2+gap_err/2+fuz << std::endl;
    std::cout << "Memuse: " << mymemuse << std::endl;
#endif

    // Iterate over the nodes at current depth
#ifdef DEBUG
    std::cout << "Border size: " << border.size() << std::endl;
#endif

    for (std::unordered_set<Node, node_hash, equal_to<Node>, count_allocator<Node> >::iterator it = border.begin(); it != border.end(); ++ it) {
#ifndef SINGLE_THREAD
      {
	LocalSynchronizer local(global_lock);
#endif
	mymemuse = memuse[id];
#ifndef SINGLE_THREAD
      }
#endif

      if (mymemuse >= max_mem)
	break;

      Node n = (Node) *it;

#ifdef DEBUG
      std::cout << graph.toString(n) << (n.strand == STRAND_FORWARD ? " F": " R") <<  std::endl;
#endif

      // Get the neighbors of the node and update their dp rows
      Graph::Vector<Node> neighbors = graph.successors(n);
      map_element *node_me = reachableSetLeft[n];
      // Number of paths to the parent node
      int num_paths = n.strand == STRAND_FORWARD ? node_me->getf(currentD-1) : node_me->getr(currentD-1);

      for (size_t i = 0; i < neighbors.size(); i++) {
	if (currentD < gap_len/2+gap_err/2+fuz || reachableSetRight.find(neighbors[i]) != reachableSetRight.end()) {
	  map_element *me = NULL;

	  if (reachableSetLeft.find(neighbors[i]) != reachableSetLeft.end()) {
	    me = reachableSetLeft[neighbors[i]];
	  }

	  if (me == NULL) {
	    me = new map_element(gap_len+gap_err+2*fuz);
	    reachableSetLeft[neighbors[i]] = me;
	  }
	  if (neighbors[i].strand == STRAND_FORWARD) {
	    me->setf(currentD, me->getf(currentD) + num_paths > MAX_PATHS ? MAX_PATHS : me->getf(currentD) + num_paths);
	  } else {
	    me->setr(currentD, me->getr(currentD) + num_paths > MAX_PATHS ? MAX_PATHS : me->getr(currentD) + num_paths);
	  }

#ifdef DEBUG
	  std::cout << "Neighbor: " << graph.toString(neighbors[i]) << (neighbors[i].strand == STRAND_FORWARD ? " F": " R") <<  std::endl;
#endif
	  // Put the neighbor to the next border
	  nextBorder.insert(neighbors[i]);
	}
      }
    }

#ifndef SINGLE_THREAD
    {
      LocalSynchronizer local(global_lock);
#endif
      mymemuse = memuse[id];
#ifndef SINGLE_THREAD
    }
#endif
    if (mymemuse >= max_mem)
      break;

    // Swap the borders
    border.clear();
    border.swap(nextBorder);

    // Add next starting k-mer if one still exists
    if (currentD <= fuz) {
      std::string kmer = kmer_left.substr(currentD, k);
      Node node = graph.buildNode(Data((char *)kmer.c_str()));
#ifdef DEBUG
      if (!graph.contains(node)) {
	std::cout << "Kmer " << kmer_left << " cannot be found in the graph!" << std::endl;
      }
#endif

#ifdef DEBUG
      std::cout << "Kmer: " << kmer << std::endl;
#endif

      if (graph.contains(node)) {
	border.insert(node);
	map_element *me = NULL;

	if (reachableSetLeft.find(node) != reachableSetLeft.end()) {
	  me = reachableSetLeft[node];
	}

	if (me == NULL) {
	  me = new map_element(gap_len+gap_err+2*fuz);
	  reachableSetLeft[node] = me;
	}

	// Initialize the dp row. Note that here we initialize to 1
	// regardless of whether the node has already been reached. This
	// will prevent counting paths going via all the starting k-mers
	// several times
	if (node.strand == STRAND_FORWARD) {
	  me->setf(currentD, 1);
	} else {
	  me->setr(currentD, 1);
	}
      }
    }

    // Check if a path has been found
    if (currentD >= gap_len+2*fuz) {
      int err = currentD - gap_len - 2*fuz;

      // Iterate over the k-mers in the right edge of the gap
      for (int j = 0; j <= fuz && count == 0; j++) {
	std::string rkmer = kmer_right.substr(j,k);
	Node rnode = graph.buildNode(Data((char *)rkmer.c_str()));

#ifdef DEBUG
	std::cout << "Right k-mer: " << graph.toString(rnode) << std::endl;
	if (!graph.contains(rnode)) {
	  std::cout << "Kmer " << kmer_right << " cannot be found in the graph!" << std::endl;
	}
#endif

	map_element *right = NULL;

	if (reachableSetLeft.find(rnode) != reachableSetLeft.end())
	  right = reachableSetLeft[rnode];
	std::vector<int> pathLengths;

	// Check if this k-mer has been reached
	if (right == NULL) {
#ifdef DEBUG
	  std::cout << "right kmer not reached" << std::endl;
#endif
	} else {
	  // +fuz: path lengths are counted starting from the first allowed left k-mer
	  // +j: j nucleotides were ignored on the right side thus widening the gap
	  int actual_gap_len1 = gap_len + fuz + j + err;
	  int actual_gap_len2 = gap_len + fuz + j - err;

	  if (rnode.strand == STRAND_FORWARD) {
	    if (right->getf(actual_gap_len1) >= 1) {
	      count = count + right->getf(actual_gap_len1) > MAX_PATHS ? MAX_PATHS : count + right->getf(actual_gap_len1);
#ifdef DEBUG
	      std::cout << "Gap filled: " << actual_gap_len1 << "\t" << right->getf(actual_gap_len1) << std::endl;
#endif
	      pathLengths.push_back(actual_gap_len1);
	    }
	  } else {
	    if (right->getr(actual_gap_len1) >= 1) {
	      count = count + right->getr(actual_gap_len1) > MAX_PATHS ? MAX_PATHS : count + right->getr(actual_gap_len1);
#ifdef DEBUG
	      std::cout << "Gap filled: " << actual_gap_len1 << "\t" << right->getr(actual_gap_len1) << std::endl;
#endif
	      pathLengths.push_back(actual_gap_len1);
	    }
	  }

	  if (actual_gap_len2 != actual_gap_len1 && actual_gap_len2 >= 0) {
	    if (rnode.strand == STRAND_FORWARD) {
	      if (right->getf(actual_gap_len2) >= 1) {
		count = count + right->getf(actual_gap_len2) > MAX_PATHS ? MAX_PATHS : count + right->getf(actual_gap_len2);
#ifdef DEBUG
		std::cout << "Gap filled: " << actual_gap_len2 << "\t" << right->getf(actual_gap_len2) << std::endl;
#endif
		pathLengths.push_back(actual_gap_len2);
	      }
	    } else {
	      if (right->getr(actual_gap_len2) >= 1) {
		count = count + right->getr(actual_gap_len2) > MAX_PATHS ? MAX_PATHS : count + right->getr(actual_gap_len2);
#ifdef DEBUG
		std::cout << "Gap filled: " << actual_gap_len2 << "\t" << right->getr(actual_gap_len2) << std::endl;
#endif
		pathLengths.push_back(actual_gap_len2);
	      }
	    }
	  }
	}

	if (count > 0 && pathLengths.size() > 0 && fill != NULL) {
	  // Number of ignored nucleotides on the right edge
	  *right_fuz = j;

	  // Recover the subgraph of the DBG covered by all the paths
	  int currentD2;
	  Node current;

	  digraph subgraph(0);
	  std::unordered_map<Node, bnode, node_hash, equal_to<Node>, count_allocator< std::pair <const Node, bnode> > > node2boost;
	  int *branch = NULL;

	  if (!skip_confident) {
	    std::unordered_set<Node, node_hash, equal_to<Node>, count_allocator<Node> > backBorder;
	    std::unordered_set<Node, node_hash, equal_to<Node>, count_allocator<Node> > nextBackBorder;

	    backBorder.insert(rnode);
	    node2boost[rnode] = boost::add_vertex(subgraph);

	    bnode sink = node2boost[rnode];
	    bnode source = boost::add_vertex(subgraph);

#ifdef DEBUG
	    std::cout << "PathLengthSize: " << pathLengths.size() << std::endl;
#endif

	    for (size_t j = 0; j < pathLengths.size(); j++) {
	      backBorder.clear();
	      nextBackBorder.clear();
	      backBorder.insert(rnode);
	      current = rnode;
	      currentD2 = pathLengths[j];

	      while(currentD2 >= 0) {
		Node lnode;
		if (currentD2 <= fuz) {
		  lnode = graph.buildNode(Data((char *)kmer_left.substr(currentD2,k).c_str()));
		}
		for(auto it = backBorder.begin(); it != backBorder.end(); ++it) {
		  current = (Node) *it;
		  // Check for end condition
		  if (currentD2 > fuz || current != lnode) {
		    Graph::Vector<Node> neighbors = graph.predecessors(current);
		    for (size_t i = 0; i < neighbors.size(); i++) {
		      Node n = neighbors[i];
		      if (reachableSetLeft.find(n) != reachableSetLeft.end()) {
			if (n.strand == STRAND_FORWARD) {
			  if (reachableSetLeft[n]->getf(currentD2-1) > 0) {
			    nextBackBorder.insert(n);
			    if (node2boost.find(n) == node2boost.end())
			      node2boost[n] = boost::add_vertex(subgraph);
			    if (!boost::edge(node2boost[n], node2boost[current], subgraph).second) {
			      boost::add_edge(node2boost[n], node2boost[current], subgraph);
			    }
			  }
			} else {
			  if (reachableSetLeft[n]->getr(currentD2-1) > 0) {
			    nextBackBorder.insert(n);
			    if (node2boost.find(n) == node2boost.end())
			      node2boost[n] = boost::add_vertex(subgraph);
			    if (!boost::edge(node2boost[n], node2boost[current], subgraph).second) {
			      boost::add_edge(node2boost[n], node2boost[current], subgraph);
			    }
			  }
			}
		      }
		    }
		  } else {
		    if (!boost::edge(source, node2boost[current], subgraph).second) {
		      boost::add_edge(source, node2boost[current], subgraph);
		    }
		  }
		}

		backBorder.clear();
		backBorder.swap(nextBackBorder);
		currentD2--;
	      }
	    }

	    std::vector<size_t> components(boost::num_vertices(subgraph));
	    size_t num_components = strong_components(subgraph, make_iterator_property_map(components.begin(), get(boost::vertex_index, subgraph), components[0]));

#ifdef DEBUG
	    std::cout << "Size of the path subgraph: " << boost::num_vertices(subgraph) << " / " <<  boost::num_edges(subgraph)  <<std::endl;
	    std::cout << "Strongly connected components: " << num_components << std::endl;
#endif

	    int csize[num_components];

	    size_t num_real_vertices = boost::num_vertices(subgraph);
	    size_t num_real_edges = boost::num_edges(subgraph);
	    size_t num_nontrivial_components = 0;
	    size_t size_nontrivial_components = 0;

	    if (num_components != boost::num_vertices(subgraph)) {
#ifdef DEBUG
	      std::cout << "Subgraph contains cycles" << std::endl;
#endif
	      for (size_t i = 0; i < num_components; i++) {
		csize[i] = 0;
	      }
	      for (size_t i = 0; i < boost::num_vertices(subgraph); i++) {
		csize[components[i]]++;
	      }
	      bnode cnode[num_components];
	      for(size_t i = 0; i < num_components; i++) {
		if (csize[i] > 1) {
		  cnode[i] = boost::add_vertex(subgraph);
		  num_nontrivial_components++;
		}
	      }
	      for (size_t i = 0; i < num_real_vertices; i++) {
		if (csize[components[i]] > 1) {
		  size_nontrivial_components++;
		  boost::graph_traits<digraph>::in_edge_iterator e, end;
		  for(tie(e,end) = in_edges(i, subgraph); e != end; ++e) {
		    if (components[i] != components[boost::source(*e, subgraph)]) {
		      if (csize[components[boost::source(*e, subgraph)]] > 1) {
			if (boost::source(*e, subgraph) < i) {
			  boost::add_edge(cnode[components[boost::source(*e, subgraph)]], cnode[components[i]], subgraph);
			}
		      } else {
			boost::add_edge(boost::source(*e, subgraph), cnode[components[i]], subgraph);
		      }
		    }
		  }
		  boost::graph_traits<digraph>::out_edge_iterator e2, end2;
		  for(tie(e2,end2) = out_edges(i, subgraph); e2 != end2; ++e2) {
		    if (components[i] != components[boost::target(*e2, subgraph)]) {
		      if (csize[components[boost::target(*e2, subgraph)]] > 1) {
			if (boost::target(*e2, subgraph) < i)
			  boost::add_edge(cnode[components[i]], cnode[components[boost::target(*e2, subgraph)]], subgraph);
		      } else {
			boost::add_edge(cnode[components[i]], boost::target(*e2, subgraph), subgraph);
		      }
		    }
		  }
		}
	      }

	      for (size_t i = 0; i < num_real_vertices; i++) {
		if (csize[components[i]] > 1) {
		  boost::clear_vertex(i, subgraph);
		}
	      }
	    } else {
	      for (size_t i = 0; i < num_components; i++) {
		csize[i] = 1;
	      }
	    }

	    substats->vertices = num_real_vertices;
	    substats->edges = num_real_edges;
	    substats->nontrivial_components = num_nontrivial_components;
	    substats->size_nontrivial_components = size_nontrivial_components;
	    substats->vertices_final = (boost::num_vertices(subgraph) - size_nontrivial_components);
	    substats->edges_final = boost::num_edges(subgraph);

	    // Count for parallel branches
	    branch = new int[boost::num_vertices(subgraph)];
	    for(size_t i = 0; i < boost::num_vertices(subgraph); i++) {
	      branch[i] = 0;
	    }

	    // Topological sorting
	    container tsorted;
	    topological_sort(subgraph, std::back_inserter(tsorted));
	    int branchcount = 1;
	    for(container::reverse_iterator it=tsorted.rbegin(); it != tsorted.rend(); ++it) {
	      bnode n = *it;
	      if (boost::in_degree(n, subgraph)  >= 1 || boost::out_degree(n, subgraph) >= 1) {
		if (in_degree(n,subgraph) > 1) {
		  branchcount -= (in_degree(n, subgraph)-1);
		}
		branch[n] = branchcount;
		if (out_degree(n, subgraph) > 1) {
		  branchcount += (out_degree(n, subgraph)-1);
		}
	      }
	    }

#ifdef DEBUG
	    write_graphviz(std::cout, subgraph, boost::make_label_writer(branch));
#endif
	  }

	  // Recover the fill sequence

	  // Choose by random the length of the path to follow
	  currentD2 = pathLengths[rand() % pathLengths.size()];
	  int lastSolid = currentD2;

	  // Trace back in the dp matrix

	  current = rnode;
	  // Set of in neighbors of the current node
	  std::vector<Node> back;
	  fill[currentD2] = '\0';

	  while(currentD2 >= 0) {
	    // The current k-mer
	    std::string str = graph.toString(current);

#ifdef DEBUG
	    std::cout << str << std::endl;
#endif

	    // Check for end condition
	    if (currentD2 <= fuz) {
	      Node lnode = graph.buildNode(Data((char *)kmer_left.substr(currentD2,k).c_str()));
#ifdef DEBUG
	      std::cout << kmer_left << " " << graph.toString(lnode) << " " << graph.toString(current) << std::endl;
#endif
	      if (lnode == current) {
		*left_fuz = fuz-currentD2;
		break;
	      }
	    }

	    // Get the in neighbors
	    if (currentD2 > 0) {
	      if (skip_confident || branch[node2boost[current]] == 1) {
		lastSolid = currentD2;
	      }
	      if (currentD2 > lastSolid - k) {
		fill[currentD2-1] = toupper(str[str.length()-1]);
	      } else {
		fill[currentD2-1] = tolower(str[str.length()-1]);
	      }
	      Graph::Vector<Node> neighbors = graph.predecessors(current);
	      for (size_t i = 0; i < neighbors.size(); i++) {
		Node n = neighbors[i];
		if (reachableSetLeft.find(n) != reachableSetLeft.end()) {
		  if (n.strand == STRAND_FORWARD) {
		    if (reachableSetLeft[n]->getf(currentD2-1) > 0) {
		      back.push_back(n);
		    }
		  } else {
		    if (reachableSetLeft[n]->getr(currentD2-1) > 0) {
		      back.push_back(n);
		    }
		  }
		}
	      }
	      // There should be in-neighbors where the path originated from but check just in case...
	      if (back.size() == 0) {
		std::cout << "Unable to backtrace! " << currentD2 << " " << currentD << " " << graph.toString(rnode) <<  std::endl;
		// Free memory
		for(auto it = reachableSetLeft.begin(); it != reachableSetLeft.end(); ++it) {
		  delete it->second;
		}
		for(auto it = reachableSetRight.begin(); it != reachableSetRight.end(); ++it) {
		  delete it->second;
		}
		if (!skip_confident) {
		  delete [] branch;
		}
		return 0;
	      }
	      // Randomly choose one of the in neighbors
	      current = back[rand()%(back.size())];
	    }
	    currentD2--;
	    back.clear();
	  }
	  if (!skip_confident) {
	    delete [] branch;
	  }
	  break;
	}
      }
    }

    currentD++;

#ifndef SINGLE_THREAD
    {
      LocalSynchronizer local(global_lock);
#endif
      mymemuse = memuse[id];
#ifndef SINGLE_THREAD
    }
#endif
  }

#ifdef DEBUG
  std::cout << std::endl;
#endif

#ifdef DEBUG
    std::cout << "Reachable set: " << std::endl;
    for(auto it = reachableSetRight.begin(); it != reachableSetRight.end(); ++it) {
      std::cout << graph.toString(it->first) << std::endl;
      for(int i = 0; i <= gap_len+gap_err+2*fuz; i++) {
	if (it->first.strand == STRAND_FORWARD) {
	  std::cout << " " << (it->second)->getf(i);
	} else {
	  std::cout << " " << (it->second)->getr(i);
	}
      }
      std::cout << std::endl;
    }
#endif

#ifdef DEBUG
    if (mymemuse > max_mem) {
      std::cout << "Memory limit exceeded - giving up on this gap." << std::endl;
    }
#endif

    // Check for memory usage to set count accordingly
    if (mymemuse > max_mem)
      count = -1;

  // Free memory
#ifdef STATS
  long long nz = 0;
  long long n = 0;
#endif
  for(auto it = reachableSetLeft.begin(); it != reachableSetLeft.end(); ++it) {
#ifdef STATS
    n++;
    nz += it->second->get_nz();
#endif
    delete it->second;
  }
  for(auto it = reachableSetRight.begin(); it != reachableSetRight.end(); ++it) {
#ifdef STATS
    n++;
    nz += it->second->get_nz();
#endif
    delete it->second;
  }

#ifdef STATS
#ifndef SINGLE_THREAD
  {
      LocalSynchronizer local(global_lock);
#endif
      std::cout << "Nodes: " << n << " Nonzeros: " << nz << " Total memory: " << memuse[id] << std::endl;
#ifndef SINGLE_THREAD
  }
#endif
#endif

  return count;
}

std::set<Node> Gap2Seq::extract_reachable_nodes(Graph graph, std::string kmer, int d) {
  std::set<Node> border;
  std::set<Node> nextBorder;
  int currentD = 0;

  std::set<Node> reachableSet;

  Node node = graph.buildNode(Data((char *)kmer.c_str()));

  if (!graph.contains(node)) {
    std::cout << "Kmer " << kmer << " cannot be found in the graph!" << std::endl;
    exit(EXIT_FAILURE);
  }

  border.insert(node);

  while(currentD <= d) {
    for (std::set<Node>::iterator it = border.begin(); it != border.end(); ++ it) {
      Node n = (Node) *it;
      reachableSet.insert(n);
      Graph::Vector<Node> neighbors = graph.successors(n);
      for (size_t i = 0; i < neighbors.size(); i++) {
	nextBorder.insert(neighbors[i]);
      }
    }

    border.clear();
    border.swap(nextBorder);
    currentD++;
  }

  return reachableSet;
}
