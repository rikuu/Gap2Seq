/*****************************************************************************
 *  GapCutter
 *  Copyright (C) Riku Walve 2015
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

#include <map>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include <gatb/gatb_core.hpp>

#define FILTER
#define SINGLE_THREAD

// Default parameter values
#define DEFAULT_K "31"
#define DEFAULT_FUZ "10"

/****************************************************************************/

// Constants for command line parameters
static const char* STR_KMER_LEN = "-k";
static const char* STR_SCAFFOLDS = "-scaffolds";
static const char* STR_CONTIGS = "-contigs";
static const char* STR_GAPS = "-gaps";
static const char* STR_FUZ = "-fuz";
static const char* STR_BED = "-bed";

static const char* STR_MASK = "-mask";
static const char* STR_SPLIT = "-split";
static const char* STR_SHORT = "-short";

/****************************************************************************/

static const char* STR_SCAFFOLD_MARKER = " scaffold ";
static const char* STR_GAP_MARKER = " gap ";

/****************************************************************************/

class GapCutter : public Tool
{
public:

    // Constructor
    GapCutter ();

    // Actual job done by the tool is here
    void execute ();

    void insertSequence(BankFasta &, const std::string &,
      const std::string &) const;

    size_t gapLength(const std::string &, const size_t) const;
    size_t distanceToNextGap(const std::string &, const size_t) const;
};


/*****************************************************************************
Constructor for the tool.
*****************************************************************************/

GapCutter::GapCutter ()  : Tool ("GapCutter")
{
    // We add some custom arguments for command line interface
    getParser()->push_front (new OptionOneParam (STR_SCAFFOLDS, "FASTA/Q file of scaffolds",  true));
    getParser()->push_front (new OptionOneParam (STR_CONTIGS, "FASTA file of contigs",  true));
    getParser()->push_front (new OptionOneParam (STR_GAPS, "FASTA file of gaps",  true));
    getParser()->push_front (new OptionOneParam (STR_BED, "BED file for gaps", true));
    getParser()->push_front (new OptionOneParam (STR_FUZ, "Number of nucleotides to ignore on gap fringes",  false, DEFAULT_FUZ));
    getParser()->push_front (new OptionOneParam (STR_KMER_LEN, "kmer length",  false, DEFAULT_K));

    getParser()->push_front (new OptionNoParam (STR_MASK, "Mask sequences too short to use",  false, false));
    getParser()->push_front (new OptionNoParam (STR_SPLIT, "Re-use flanks too short to split",  false, false));
    getParser()->push_front (new OptionNoParam (STR_SHORT, "Output flanks less than k+fuz",  false, false));
}

void GapCutter::insertSequence(BankFasta &bank, const std::string &comment,
  const std::string &sequence) const
{
  Sequence seq((char *) sequence.c_str());
  seq._comment = comment;
  bank.insert(seq);
}

size_t GapCutter::distanceToNextGap(const std::string &seq, const size_t start) const
{
  size_t distance = 0;
  for (size_t i = start+1; i < seq.size(); i++) {
    distance++;
    if (seq[i] == 'N' || seq[i] == 'n') {
      break;
    }
  }

  return distance;
}

size_t GapCutter::gapLength(const std::string &seq, const size_t start) const
{
  size_t gapLength = 0;
  for (size_t i = start; i < seq.size(); i++) {
    gapLength++;

    if (seq[i] != 'N' && seq[i] != 'n') {
      break;
    }
  }

  return gapLength;
}

void GapCutter::execute ()
{
  // Get the command line arguments
  const std::string scaffoldsFilename = getInput()->getStr(STR_SCAFFOLDS);
  const std::string contigsFilename = getInput()->getStr(STR_CONTIGS);
  const std::string gapsFilename = getInput()->getStr(STR_GAPS);
  const std::string bedFilename = getInput()->getStr(STR_BED);
  const int k = getInput()->getInt(STR_KMER_LEN);
  const int fuz = getInput()->getInt(STR_FUZ);

  const bool mask = (getInput()->get(STR_MASK) != 0);
  const bool split = (getInput()->get(STR_SPLIT) != 0);
  const bool nonstatic = (getInput()->get(STR_SHORT) != 0);

  std::cout << "Scaffolds file: " << scaffoldsFilename << std::endl;
  std::cout << "Contigs file: " << contigsFilename << std::endl;
  std::cout << "Gaps file: " << gapsFilename << std::endl;
  std::cout << "BED file: " << bedFilename << std::endl;
  std::cout << "k-mer size: " << k << std::endl;
  std::cout << "Fuz: " << fuz << std::endl;

  // Open the input scaffold file
  BankFasta scaffoldBank(scaffoldsFilename);
  BankFasta::Iterator itSeq(scaffoldBank);

  // Open the output files
  BankFasta contigBank(contigsFilename);
  BankFasta gapBank(gapsFilename);

  // Open output stream for BED file
  ofstream bedFile;
  bedFile.open(bedFilename);

  int contig = 0;
  int gap = 0;
  int scaffold = 0;
  int allGaps = 0;

  for (itSeq.first(); !itSeq.isDone(); itSeq.next()) {
    const std::string seq = itSeq->toString();
    const std::string comment = itSeq->getComment() + STR_SCAFFOLD_MARKER + std::to_string(scaffold);

    size_t i = 0;
    while (i < seq.size()) {
      const size_t d1 = distanceToNextGap(seq, i);
      const size_t l1 = gapLength(seq, i + d1);

      // Not enough sequence for left flank, skip to next possible flank position
      if (d1 < k) {
        if (d1 == 0) {
          i = seq.size();
        }

        i += d1 + l1;
        continue;
      }

      const size_t flank1 = min(d1, k+fuz);

      const size_t d2 = distanceToNextGap(seq, i + d1 + l1);
      const size_t l2 = gapLength(seq, i + d1 + l1 + d2);

      // Case 1: Gap has enough sequence on both sides
      if (d2 >= 2*k) {
        const size_t flank2 = min(d2-k, k+fuz);

        // Write the gap to file
        const std::string gapComment = comment + STR_GAP_MARKER + std::to_string(gap);
        insertSequence(gapBank, gapComment, seq.substr(i + d1 - flank1, flank1 + l1 + flank2));
        insertSequence(contigBank, gapComment, seq.substr(i, d1 - flank1));

        // Write the gap information to BED file
        bedFile << contigName << "\t" << i + d1 - flank1 << "\t" << i + d1 + l1 + flank2 << std::endl;

        gap++;
        contig++;

        i += d1 + l1 + flank2;
        continue;
      }

      const size_t d3 = distanceToNextGap(seq, i + d1 + l1 + d2 + l2);

      // Case 2: Right flank overlaps with next left flank
      if (d2 >= k) {
        assert(d3 >= k);
        const size_t flank3 = min(d3, k+fuz);

        const std::string gapComment = comment + STR_GAP_MARKER + std::to_string(gap);
        insertSequence(gapBank, gapComment + STR_SPLIT_MARKER + "1", seq.substr(i + d1 - flank1, flank1 + l1 + d2));
        insertSequence(gapBank, gapComment + STR_SPLIT_MARKER + "2", seq.substr(i + d1 + l1, d2 + l2 + flank3));
        insertSequence(contigBank, gapComment, seq.substr(i, d1 - flank1));

        bedFile << contigName << "\t" << i + d1 - flank1 << "\t" << i + d1 + l1 + d2 << std::endl;
        bedFile << contigName << "\t" << i + d1 + l1 << "\t" << i + d1 + l1 + d2 + l2 + flank3 << std::endl;

        gap++;
        contig++;

        i += d1 + l1 + d2 + l2 + flank3;
        continue;
      }

      // Case 3: Two (or more) gaps encircle a sequence(s) too short to be a flank(s)
      size_t lsum, dn = d1 + l1 + d2 + l2, d3;
      while (dn > 0 && dn < k) {
        lsum += dn + gapLength(seq, i + lsum + dn);
        dn = distanceToNextGap(seq, i + lsum);
      }

      assert(dn >= k);
      const size_t flank2 = min(dn, k+fuz);

      // Mask everything between flanks as a gap
      char *buffer = new char[lsum];
      for (size_t j = 0; j < lsum; j++) {
        buffer[j] = 'n';
      }

      const std::string gapComment = comment + STR_GAP_MARKER + std::to_string(gap);
      insertSequence(gapBank, gapComment, seq.substr(i + d1 - flank1, flank1) + std::to_string(buffer) + seq.substr(i + lsum, flank2));
      insertSequence(contigBank, gapComment, seq.substr(i, d1 - flank1));

      delete[] buffer;

      // Write the gap information to BED file
      bedFile << contigName << "\t" << i + d1 - flank1 << "\t" << i + lsum + flank2 << std::endl;

      gap++;
      contig++;

      i += lsum + flank2;
    }

    gapBank.flush();
    contigBank.flush();

    scaffold++;
  }

  std::cout << "Found " << allGaps << " gaps" << std::endl;

  std::cout << "Cut out " << scaffold << " scaffolds into " <<
    contig << " contigs and " << gap <<  " gaps" << std::endl;
}

/****************************************************************************/

int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.
        GapCutter().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
