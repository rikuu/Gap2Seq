/*****************************************************************************
 *  Gap2Seq
 *  Copyright (C) Leena Salmela, Kristoffer Sahlin, Veli Mäkinen,
 *  Alexandru Tomescu, Riku Walve 2017
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

#include <sys/types.h>

#include <string>
#include <vector>
#include <functional>

// GATB-core Bloom filter requires hash1 function for items
inline u_int64_t hash1(const std::string &key, u_int64_t seed=0) {
  return std::hash<std::string>{}(key);
}

#include <gatb/gatb_core.hpp>

#include <htslib/sam.h>

/*****************************************************************************/

static const char* STR_ALIGNMENT = "-bam";
static const char* STR_OUTPUT = "-reads";

static const char* STR_READ_LENGTH = "-read-length";
static const char* STR_MEAN = "-mean";
static const char* STR_STD_DEV = "-std-dev";

static const char* STR_SCAFFOLD = "-scaffold";
static const char* STR_GAP_BREAKPOINT = "-breakpoint";

static const char* STR_FLANK_LENGTH = "-flank-length";

static const char* STR_GAP_LENGTH = "-gap-length";
static const char* STR_THRESHOLD = "-unmapped";

static const char* STR_ONLY_UNMAPPED = "-unmapped-only";

/*****************************************************************************/

class io_t {
 public:
  samFile *sam = NULL;
  bam_hdr_t *header = NULL;
  hts_idx_t *idx = NULL;
  bool loaded = false;

  io_t(const std::string &);

  void unload() {
    bam_hdr_destroy(this->header);
    sam_close(this->sam);
    loaded = false;
  }
};

// Loads a bam/sam file into an IO object
io_t::io_t(const std::string &samFilename) {
  this->sam = sam_open(samFilename.c_str(), "r");
  if (this->sam == NULL) {
    // std::cerr << "ERROR: SAM file not found!" << std::endl;
    return;
  }

  this->header = sam_hdr_read(this->sam);
  if (this->header == NULL) {
    // std::cerr << "ERROR: SAM header error!" << std::endl;
    return;
  }

  this->idx = sam_index_load(this->sam, samFilename.c_str());
  if (this->idx == NULL) {
    // std::cerr << "ERROR: SAM index not found!" << std::endl;
    return;
  }

  this->loaded = true;
}

/*****************************************************************************/

inline uint8_t complement(const uint8_t n) {
  switch (n) {
    case 1: return 8; break;
    case 2: return 4; break;
    case 4: return 2; break;
    case 8: return 1; break;
    case 15:
    default: return 15; break;
  }
}

inline uint8_t querySequence(const uint8_t *query, const int32_t length,
    const int32_t index, const bool reverse) {
  if (!reverse)
    return bam_seqi(query, index);

  return complement(bam_seqi(query, length - 1 - index));
}

/*****************************************************************************/

// Converts an alignment to std::string. Handles reverse complements.
std::string convertToString(const uint8_t *query, const int32_t length,
    const bool reverse, char* buffer) {
  for (int32_t i = 0; i < length; i++) {
    switch (querySequence(query, length, i, reverse)) {
      case 0x1:
        buffer[i] = 'A';
        break;
      case 0x2:
        buffer[i] = 'C';
        break;
      case 0x4:
        buffer[i] = 'G';
        break;
      case 0x8:
        buffer[i] = 'T';
        break;
      case 0x15:
      default:
        buffer[i] = 'N';
        break;
    }
  }

  buffer[length] = '\0';
  return std::string(buffer);
}

/*****************************************************************************/

static inline const std::string bam2string(const bam1_t *bam) {
  return std::string(bam_get_qname(bam)) + (((bam->core.flag & BAM_FREAD1) != 0) ? "/1" : "/2");
}

static inline const std::string bam2string_mate(const bam1_t *bam) {
  return std::string(bam_get_qname(bam)) + (((bam->core.flag & BAM_FREAD1) != 0) ? "/2" : "/1");
}

/*****************************************************************************/

class sam_iterator {
private:
  samFile *m_sam;
  hts_itr_t *m_iter;

public:
  bam1_t *bam;

  sam_iterator(const io_t io, const int tid, const int start, const int end)
      : m_sam(io.sam), bam(bam_init1()) {
    m_iter = sam_itr_queryi(io.idx, tid, start, end);
    if (m_iter == NULL) {
      std::cerr << "WARNING: SAM iterator is NULL!" << std::endl;
    }
  }

  sam_iterator(const io_t io, const char *string)
      : m_sam(io.sam), bam(bam_init1()) {
    m_iter = sam_itr_querys(io.idx, io.header, string);
    if (m_iter == NULL) {
      std::cerr << "WARNING: SAM iterator is NULL!" << std::endl;
    }
  }

  ~sam_iterator() {
    hts_itr_destroy(m_iter);
    bam_destroy1(bam);
  }

  inline bool next() {
    if (m_iter == NULL) return false;
    return (sam_itr_next(m_sam, m_iter, bam) >= 0);
  }
};

// Counts the number of reads by iterating through the alignment file
uint64_t count_reads(const std::string &filename) {
  io_t io(filename);
  sam_iterator iter(io, ".");

  uint64_t count = 0;
  while (iter.next()) {
    count++;
  }

  io.unload();

  return count;
}

/*****************************************************************************/

class ReadFilter : public Tool {
 public:
    ReadFilter();
    void execute();

    // Prints an alignment in fasta format
    void print_fasta(const bam1_t*, char*, BankFasta*);

    // Prints all alignments in a region
    void process_region(const io_t&, const int, const int, const int, char*, IBloom<std::string>*, BankFasta*, int*, int*);
    void process_mates(const io_t&, const int, const int, const int, IBloom<std::string>*);
    void find_mates(const io_t&, char*, IBloom<std::string>*, BankFasta*, int*, int*);
    void process_unmapped(const io_t&, char*, IBloom<std::string>*, BankFasta*, int*);
};

ReadFilter::ReadFilter() : Tool("ReadFilter") {
  // Input / output
  getParser()->push_front(new OptionOneParam(STR_ALIGNMENT, "Aligned BAM file", true));
  getParser()->push_front(new OptionOneParam(STR_OUTPUT, "FASTA-formatted output file", true));

  // Read library parameters
  getParser()->push_front(new OptionOneParam(STR_READ_LENGTH, "Read length", true));
  getParser()->push_front(new OptionOneParam(STR_MEAN, "Mean insert size", true));
  getParser()->push_front(new OptionOneParam(STR_STD_DEV, "Insert size standard deviation", true));

  // Gap parameters
  getParser()->push_front(new OptionOneParam(STR_SCAFFOLD, "Scaffold name", true));
  getParser()->push_front(new OptionOneParam(STR_GAP_BREAKPOINT, "Gap position", true));

  getParser()->push_front(new OptionOneParam(STR_FLANK_LENGTH , "Flank length", false, "-1"));

  getParser()->push_front(new OptionOneParam(STR_GAP_LENGTH, "Gap length", false, "-1"));
  getParser()->push_front(new OptionOneParam(STR_THRESHOLD, "Threshold for using unmapped reads", false, "-1"));

  getParser()->push_front(new OptionNoParam(STR_ONLY_UNMAPPED, "Only output unmapped reads"));
}

/*****************************************************************************/

// Output a bam object into GATB fasta bank
void ReadFilter::print_fasta(const bam1_t *bam, char* buffer, BankFasta *bank) {
  const std::string sequence = convertToString(bam_get_seq(bam),
    bam->core.l_qseq, bam_is_rev(bam), buffer);

  Sequence seq(buffer);
  seq._comment = bam2string(bam);
  bank->insert(seq);
}

// Output reads that map to a region in a scaffold and not in the Bloom filter.
void ReadFilter::process_region(const io_t &io, const int tid, const int start,
    const int end, char* buffer, IBloom<std::string> *bloom, BankFasta *bank,
    int *seqlen, int *num_of_reads) {
  sam_iterator iter(io, tid, start, end);
  while (iter.next()) {
    if (!bloom->contains(bam2string(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      *seqlen += strlen(buffer);
      (*num_of_reads)++;
    }
  }
}

void ReadFilter::process_mates(const io_t &io, const int tid, const int start,
    const int end, IBloom<std::string> *bloom) {
  sam_iterator iter(io, tid, start, end);
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FMUNMAP) != 0) {
      bloom->insert(bam2string(iter.bam));
    }
  }
}

void ReadFilter::find_mates(const io_t &io, char *buffer, IBloom<std::string> *bloom,
    BankFasta *bank, int *seqlen, int *num_of_reads) {
  sam_iterator iter(io, ".");
  while (iter.next()) {
    if (bloom->contains(bam2string_mate(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      *seqlen += strlen(buffer);
      (*num_of_reads)++;
    }
  }
}

// Output all unmapped reads
void ReadFilter::process_unmapped(const io_t &io, char* buffer,
    IBloom<std::string> *bloom, BankFasta *bank, int* num_of_reads) {
  sam_iterator iter(io, ".");
  while (iter.next()) {
    if ((iter.bam->core.flag & BAM_FUNMAP) != 0 &&
          !bloom->contains(bam2string(iter.bam))) {
      print_fasta(iter.bam, buffer, bank);
      (*num_of_reads)++;
    }
  }
}

// Execute read ReadFilterion
void ReadFilter::execute() {
  const std::string alignment = getInput()->getStr(STR_ALIGNMENT);
  const std::string output = getInput()->getStr(STR_OUTPUT);

  const int read_length = static_cast<int>(getInput()->getInt(STR_READ_LENGTH));
  const int mean_insert = static_cast<int>(getInput()->getInt(STR_MEAN));
  const int std_dev = static_cast<int>(getInput()->getInt(STR_STD_DEV));

  const std::string scaffold = getInput()->getStr(STR_SCAFFOLD);
  const int breakpoint = static_cast<int>(getInput()->getInt(STR_GAP_BREAKPOINT));

  const int flank_length = static_cast<int>(getInput()->getInt(STR_FLANK_LENGTH));

  const int gap_length = static_cast<int>(getInput()->getInt(STR_GAP_LENGTH));
  const int threshold = static_cast<int>(getInput()->getInt(STR_THRESHOLD));

  const bool unmapped_only = getParser()->saw(STR_ONLY_UNMAPPED);

  // Load alignment file
  io_t io(alignment);
  if (!io.loaded) {
    std::cerr << "Error loading alignments" << std::endl;
    return;
  }

  // Allocate memory for string conversions
  char *buffer = new char[read_length+1];

  // Use basic Bloom filter from GATB
  const uint64_t num_of_reads = count_reads(alignment);
  IBloom<std::string> *bloom = new BloomSynchronized<std::string>(5 * num_of_reads);

  // Open output file
  BankFasta reads(output);
  int seqlen = 0, reads_extracted = 0;

  if (!unmapped_only) {
    // Compute scaffold id from scaffold name
    const int tid = bam_name2id(io.header, scaffold.c_str());

    // Extract pairs from the left mappings
    const int left_start = breakpoint - (mean_insert + 3*std_dev + 2*read_length);
    const int left_end = breakpoint - (mean_insert - 3*std_dev + read_length);
    process_mates(io, tid, left_start, left_end, bloom);

    // Extract pairs from the right mappings
    const int right_start = breakpoint + (mean_insert + 3*std_dev + read_length);
    const int right_end = breakpoint + (mean_insert - 3*std_dev + read_length);
    process_mates(io, tid, right_start, right_end, bloom);

    // Output reads and count length
    find_mates(io, buffer, bloom, &reads, &seqlen, &reads_extracted);

    // Output overlapping reads
    if (flank_length != -1) {
      const int start = breakpoint - flank_length;
      const int end = breakpoint + flank_length;
      process_region(io, tid, start, end, buffer, bloom, &reads, &seqlen, &reads_extracted);
    }
  }

  // Output unmapped reads
  if (unmapped_only || ((gap_length != -1 && threshold != -1) && ((seqlen / gap_length) < threshold))) {
    process_unmapped(io, buffer, bloom, &reads, &reads_extracted);
  }

  std::cout << "Extracted " << reads_extracted << " out of " << num_of_reads << " reads" << std::endl;

  // Cleanup
  reads.flush();
  delete[] buffer;
  delete bloom;
  io.unload();
}

int main(int argc, char* argv[]) {
  try {
    ReadFilter().run(argc, argv);
  } catch (Exception& e) {
    std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
