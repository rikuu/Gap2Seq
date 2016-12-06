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
    getParser()->push_front (new OptionOneParam( STR_BED, "BED file for gaps", true));
    getParser()->push_front (new OptionOneParam (STR_FUZ, "Number of nucleotides to ignore on gap fringes",  false, DEFAULT_FUZ));
    getParser()->push_front (new OptionOneParam (STR_KMER_LEN, "kmer length",  false, DEFAULT_K));
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

    size_t previousGap = 0;

    for (size_t i = 0; i < seq.size(); i++) {
      if (seq[i] == 'N' || seq[i] == 'n') {
        allGaps++;

        #ifdef FILTER
          const size_t startPosition = i-k-fuz;
        #else
          const size_t startPosition = (i < k+fuz) ? i-k : i-k-fuz;
        #endif

        // The gap is too close to the beginning, consider it a contig
        if (previousGap > startPosition) {
          std::cout << "skipping " << comment << std::endl;
          i += gapLength(seq, i);
          continue;
        }

        // Count the gap length
        size_t length = gapLength(seq, i);

        // Merge gaps that are too close to each other
        const size_t distance = distanceToNextGap(seq, i + length);
        if (distance < (k+fuz)) {
          std::cout << "merging " << comment << std::endl;
          length += distance + gapLength(seq, i);
        }

        // The gap is too close to the end, consider it a contig
        #ifdef FILTER
          const size_t sequenceLength = (k+fuz)*2 + length;
          if ((startPosition + sequenceLength) >= seq.size()) {
            std::cout << "skipping " << comment << std::endl;
            break;
          }
        #else
          const size_t sequenceLength = ((startPosition + sequenceLength) >= seq.size()) ?
            k*2 + fuz + length : (k+fuz)*2 + length;
          if ((startPosition + sequenceLength - fuz) >= seq.size()) {
            break;
          }
        #endif

        // Write the gap to file
        const std::string gapComment = comment + STR_GAP_MARKER + std::to_string(gap);
        insertSequence(gapBank, gapComment, seq.substr(startPosition, sequenceLength));
        insertSequence(contigBank, gapComment, seq.substr(previousGap, startPosition - previousGap));

        // Write the gap information to BED file
        bedFile << contigName << "\t" << startPosition << "\t" << startPosition + sequenceLength << std::endl;

        gap++;
        contig++;

        previousGap = startPosition + sequenceLength;
        i = previousGap;
      }
    }

    // Write the ending contig of the scaffold to file
    if (previousGap < seq.size()) {
      // std::cout << "writing " << comment << std::endl;
      insertSequence(contigBank, comment, seq.substr(previousGap));
      contig++;
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
