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

#include <string>
#include <fstream>

#include <gatb/gatb_core.hpp>

// Default parameter values
#define DEFAULT_K "31"
#define DEFAULT_FUZ "10"
#define DEFAULT_MASK false
#define DEFAULT_NO_SPLIT false

/****************************************************************************/

// Constants for command line parameters
static const char *STR_KMER_LEN = "-k";
static const char *STR_SCAFFOLDS = "-scaffolds";
static const char *STR_CONTIGS = "-contigs";
static const char *STR_GAPS = "-gaps";
static const char *STR_FUZ = "-fuz";
static const char *STR_BED = "-bed";
static const char *STR_MASK = "-mask";
static const char *STR_NO_SPLIT = "-no-split";

/****************************************************************************/

static const char *STR_SCAFFOLD_MARKER = " scaffold ";
static const char *STR_CONTIG_MARKER = " contig ";
static const char *STR_GAP_MARKER = " gap ";
static const char *STR_SPLIT_MARKER = " split ";

/****************************************************************************/

class GapCutter : public Tool
{
public:
  GapCutter();
  void execute();

  void insertSequence(BankFasta &, const std::string &,
                      const std::string &) const;

  size_t gapLength(const std::string &, const size_t) const;
  size_t distanceToNextGap(const std::string &, const size_t) const;
};

/*****************************************************************************
Constructor for the tool.
*****************************************************************************/

GapCutter::GapCutter() : Tool("GapCutter")
{
  getParser()->push_front(new OptionOneParam(STR_SCAFFOLDS, "FASTA/Q file of scaffolds", true));
  getParser()->push_front(new OptionOneParam(STR_CONTIGS, "FASTA file of contigs", true));
  getParser()->push_front(new OptionOneParam(STR_GAPS, "FASTA file of gaps", true));
  getParser()->push_front(new OptionOneParam(STR_BED, "BED file for gaps", true));
  getParser()->push_front(new OptionOneParam(STR_FUZ, "Maximum number of nucleotides to ignore on gap fringes", false, DEFAULT_FUZ));
  getParser()->push_front(new OptionOneParam(STR_KMER_LEN, "k-mer length", false, DEFAULT_K));
  getParser()->push_front(new OptionNoParam(STR_MASK, "Mask sequences too short to use", false, DEFAULT_MASK));
  getParser()->push_front(new OptionNoParam(STR_NO_SPLIT, "Don't split flank sharing gaps", false, DEFAULT_NO_SPLIT));
}

void GapCutter::insertSequence(BankFasta &bank, const std::string &comment,
                               const std::string &sequence) const
{
  Sequence seq((char *)sequence.c_str());
  seq._comment = comment;
  bank.insert(seq);
}

size_t GapCutter::distanceToNextGap(const std::string &seq, const size_t start) const
{
  size_t distance = 0;
  for (size_t i = start; i < seq.size(); i++) {
    if (seq[i] == 'N' || seq[i] == 'n') {
      break;
    }

    distance++;
  }

  return distance;
}

size_t GapCutter::gapLength(const std::string &seq, const size_t start) const
{
  size_t gapLength = 0;
  for (size_t i = start; i < seq.size(); i++) {
    if (seq[i] != 'N' && seq[i] != 'n') {
      break;
    }

    gapLength++;
  }

  return gapLength;
}

void GapCutter::execute()
{
  const std::string scaffoldsFilename = getInput()->getStr(STR_SCAFFOLDS);
  const std::string contigsFilename = getInput()->getStr(STR_CONTIGS);
  const std::string gapsFilename = getInput()->getStr(STR_GAPS);
  const std::string bedFilename = getInput()->getStr(STR_BED);
  const size_t k = (size_t) getInput()->getInt(STR_KMER_LEN);
  const size_t fuz = (size_t) getInput()->getInt(STR_FUZ);
  const bool mask = (getInput()->get(STR_MASK) != 0);
  const bool no_split = (getInput()->get(STR_NO_SPLIT) != 0);

  std::cout << "Scaffolds file: " << scaffoldsFilename << std::endl;
  std::cout << "Contigs file: " << contigsFilename << std::endl;
  std::cout << "Gaps file: " << gapsFilename << std::endl;
  std::cout << "BED file: " << bedFilename << std::endl;
  std::cout << "k-mer size: " << k << std::endl;
  std::cout << "Fuz: " << fuz << std::endl;
  std::cout << "Mask: " << mask << std::endl;
  std::cout << "Split: " << !no_split << std::endl;

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

  for (itSeq.first(); !itSeq.isDone(); itSeq.next()) {
    const std::string seq = itSeq->toString();

    // Get contig identifier from the sequence comment
    const size_t contigNamePos = itSeq->getComment().find(" ");
    assert(contigNamePos != std::string::npos);
    const std::string contigName = itSeq->getComment().substr(0, contigNamePos);

    size_t i = 0;
    while (i < seq.size()) {
      const std::string comment = itSeq->getComment() +
                                  STR_SCAFFOLD_MARKER + std::to_string(scaffold) +
                                  STR_CONTIG_MARKER + std::to_string(contig);

      const std::string gapComment = comment +
                                     STR_GAP_MARKER + std::to_string(gap);

      const size_t d1 = distanceToNextGap(seq, i);
      const size_t l1 = gapLength(seq, i + d1);

      // No gaps left, write contig
      if (d1 > 0 && l1 == 0) {
        insertSequence(contigBank, comment, seq.substr(i));
        contig++;

        i = seq.size();
        continue;
      }

      // Not enough sequence for left flank, skip to next possible flank position
      if (d1 < k) {
        insertSequence(contigBank, comment, seq.substr(i, d1 + l1));
        contig++;

        i += d1 + l1;
        continue;
      }

      const size_t flank1 = std::min(d1, k + fuz);

      const size_t d2 = distanceToNextGap(seq, i + d1 + l1);
      const size_t l2 = gapLength(seq, i + d1 + l1 + d2);

      const size_t d3 = distanceToNextGap(seq, i + d1 + l1 + d2 + l2);

      // Case 1: Gap has enough sequence on both sides
      if (d2 >= 2 * k || (d2 >= k && d3 == 0)) {
        // If last gap, add the rest of the scaffold as right flank
        const size_t flank2 = (d2 >= 2 * k) ? std::min(d2, k + fuz) : d2;

        // Write the gap to file
        insertSequence(gapBank, gapComment, seq.substr(i + d1 - flank1, flank1 + l1 + flank2));
        insertSequence(contigBank, gapComment, seq.substr(i, d1 - flank1));

        // Write the gap information to BED file
        bedFile << contigName << "\t" << i + d1 - flank1 << "\t" << i + d1 + l1 + flank2 << std::endl;

#ifdef DEBUG
        std::cout << "Case1: " << comment << " d1: " << d1 << " d2: " << d2 << " l1: " << l1 << " l2: " << l2 << " flank1: " << flank1 << " flank2: " << flank2 << " start: " << i + d1 - flank1 << " end: " << i + d1 + l1 + flank2 << " length: " << d1 + l1 + flank2 << std::endl;
#endif

        gap++;
        contig++;

        i += d1 + l1 + flank2;
        continue;
      }

      // Case 2: Two gaps have overlapping flanks
      if (d2 >= k) {
        if (!no_split && d3 >= k) {
          const size_t flank3 = std::min(d3, k + fuz);

          // First gap gets entire middle section as right flank, second gap gets
          // only k length right flank. Unfair and greedy, but makes merging them
          // back easier.
          insertSequence(gapBank, gapComment + STR_SPLIT_MARKER + "1",
                         seq.substr(i + d1 - flank1, flank1 + l1 + d2));
          insertSequence(gapBank, gapComment + STR_SPLIT_MARKER + "2 " + std::to_string(k),
                         seq.substr(i + d1 + l1 + d2 - k, k + l2 + flank3));
          insertSequence(contigBank, gapComment, seq.substr(i, d1 - flank1));

          bedFile << contigName << "\t" << i + d1 - flank1 << "\t" << i + d1 + l1 + d2 << std::endl;
          bedFile << contigName << "\t" << i + d1 + l1 << "\t" << i + d1 + l1 + d2 + l2 + flank3 << std::endl;

#ifdef DEBUG
          std::cout << "Case2: " << comment << " d1: " << d1 << " d2: " << d2 << " d3: " << d3 << " l1: " << l1 << " l2: " << l2 << " flank1: " << flank1 << " flank3: " << flank3 << std::endl;
#endif

          gap++;
          contig++;

          i += d1 + l1 + d2 + l2 + flank3;
          continue;
        } else {
          const size_t flank2 = std::min(d2, k + fuz);

          insertSequence(gapBank, gapComment, seq.substr(i + d1 - flank1, flank1 + l1 + flank2));
          insertSequence(contigBank, gapComment, seq.substr(i, d1 - flank1));

          bedFile << contigName << "\t" << i + d1 - flank1 << "\t" << i + d1 + l1 + flank2 << std::endl;

          gap++;
          contig++;

          i += d1 + l1 + flank2;
          continue;
        }
      }

      // Case 3: Two (or more) gaps surround a sequence(s) too short to be a
      // flank(s)
      size_t lsum = l1 + d2 + l2, dn = d3;
      while (dn > 0 && dn < k) {
        lsum += dn + gapLength(seq, i + d1 + lsum + dn);
        dn = distanceToNextGap(seq, i + d1 + lsum);
      }

      // Last sequence not long enough to be a flank
      if (dn < k) {
        insertSequence(contigBank, comment, seq.substr(i));
        contig++;

        i = seq.size();
        continue;
      }

      // Mask everything between flanks as a gap
      if (mask) {
        char *buffer = new char[lsum + 1];
        for (size_t j = 0; j < lsum; j++) {
          buffer[j] = 'n';
        }
        buffer[lsum] = '\0';

        const size_t flank2 = std::min(dn, k + fuz);
        insertSequence(gapBank, gapComment, seq.substr(i + d1 - flank1, flank1) + std::string(buffer) + seq.substr(i + d1 + lsum, flank2));
        insertSequence(contigBank, gapComment, seq.substr(i, d1 - flank1));

        bedFile << contigName << "\t" << i + d1 - flank1 << "\t" << i + d1 + lsum + flank2 << std::endl;

        delete[] buffer;

#ifdef DEBUG
        std::cout << "Case3: " << comment << " d1: " << d1 << " d2: " << d2 << " d3: " << d3 << " l1: " << l1 << " l2: " << l2 << " lsum: " << lsum << " flank1: " << flank1 << " flank2: " << flank2 << std::endl;
#endif

        gap++;
        contig++;

        i += d1 + lsum + flank2;
      } else {
        insertSequence(contigBank, comment, seq.substr(i, d1 + lsum));
        contig++;

        i += d1 + lsum;
      }
    }

    gapBank.flush();
    contigBank.flush();

    scaffold++;
  }

  std::cout << "Cut " << scaffold << " scaffolds into " << contig << " contigs and " << gap << " gaps" << std::endl;
}

/****************************************************************************/

int main(int argc, char *argv[])
{
  try {
    GapCutter().run(argc, argv);
  } catch (Exception &e) {
    std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
