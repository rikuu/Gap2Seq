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

/****************************************************************************/

// Constants for command line parameters
static const char *STR_SCAFFOLDS = "-scaffolds";
static const char *STR_CONTIGS = "-contigs";
static const char *STR_GAPS = "-gaps";

/****************************************************************************/

static const char *STR_SCAFFOLD_MARKER = " scaffold ";
static const char *STR_CONTIG_MARKER = " contig ";
static const char *STR_GAP_MARKER = " gap ";
static const char *STR_SPLIT_MARKER = " split ";

static const size_t SCAFFOLD_MARKER_LENGTH = 10u;
static const size_t CONTIG_MARKER_LENGTH = 8u;
static const size_t GAP_MARKER_LENGTH = 5u;
static const size_t SPLIT_MARKER_LENGTH = 7u;

/****************************************************************************/

class GapMerger : public Tool
{
public:
  GapMerger();
  void execute();

  int parseGapIndex(const std::string &) const;
  int parseScaffoldIndex(const std::string &) const;
  int parseContigIndex(const std::string &) const;

  void insertSequence(BankFasta &, const std::string &,
                      const std::string &) const;
};

/*****************************************************************************
Constructor for the tool.
*****************************************************************************/

GapMerger::GapMerger() : Tool("GapMerger")
{
  getParser()->push_front(new OptionOneParam(STR_SCAFFOLDS, "FASTA file of merged scaffolds", true));
  getParser()->push_front(new OptionOneParam(STR_CONTIGS, "FASTA file of contigs", true));
  getParser()->push_front(new OptionOneParam(STR_GAPS, "FASTA file of filled gaps", true));
}

// TODO: Refactor these into a single function
int GapMerger::parseGapIndex(const std::string &comment) const
{
  const size_t gap_pos = comment.find(STR_GAP_MARKER);
  if (gap_pos == std::string::npos) {
    return -1;
  }

  std::string index = "";
  const size_t split_pos = comment.find(STR_SPLIT_MARKER);
  if (split_pos == std::string::npos) {
    index = comment.substr(gap_pos + GAP_MARKER_LENGTH);
  } else {
    index = comment.substr(gap_pos + GAP_MARKER_LENGTH,
                           split_pos - (gap_pos + GAP_MARKER_LENGTH));
  }

  return std::stoi(index);
}

int GapMerger::parseContigIndex(const std::string &comment) const
{
  const size_t contig_pos = comment.find(STR_CONTIG_MARKER);
  if (contig_pos == std::string::npos) {
    return -1;
  }

  std::string index = "";
  const size_t gap_pos = comment.find(STR_GAP_MARKER);
  if (gap_pos == std::string::npos) {
    index = comment.substr(contig_pos + CONTIG_MARKER_LENGTH);
  } else {
    index = comment.substr(contig_pos + CONTIG_MARKER_LENGTH,
                           gap_pos - (contig_pos + CONTIG_MARKER_LENGTH));
  }

  return std::stoi(index);
}

int GapMerger::parseScaffoldIndex(const std::string &comment) const
{
  const size_t scaffold_pos = comment.find(STR_SCAFFOLD_MARKER);
  if (scaffold_pos == std::string::npos) {
    return -1;
  }

  const size_t contig_pos = comment.find(STR_CONTIG_MARKER);
  assert(contig_pos != std::string::npos);

  return std::stoi(comment.substr(scaffold_pos + SCAFFOLD_MARKER_LENGTH,
                                  contig_pos - (scaffold_pos + SCAFFOLD_MARKER_LENGTH)));
}

void GapMerger::insertSequence(BankFasta &bank, const std::string &comment,
                               const std::string &sequence) const
{
  Sequence seq((char *)sequence.c_str());

  // Remove markers from comment, i.e. remove anything after scaffold marker
  const size_t marker_pos = comment.find(STR_SCAFFOLD_MARKER);
  if (marker_pos == std::string::npos) {
    seq._comment = comment;
  } else {
    seq._comment = comment.substr(0, marker_pos);
  }

  bank.insert(seq);
}

void GapMerger::execute()
{
  const std::string scaffoldsFilename = getInput()->getStr(STR_SCAFFOLDS);
  const std::string contigsFilename = getInput()->getStr(STR_CONTIGS);
  const std::string gapsFilename = getInput()->getStr(STR_GAPS);

  std::cout << "Scaffolds file: " << scaffoldsFilename << std::endl;
  std::cout << "Contigs file: " << contigsFilename << std::endl;
  std::cout << "Gaps file: " << gapsFilename << std::endl;

  // Open the input files
  BankFasta contigBank(contigsFilename);
  BankFasta gapBank(gapsFilename);

  BankFasta::Iterator contigIter(contigBank);
  BankFasta::Iterator gapIter(gapBank);

  // Open the output scaffold file
  BankFasta scaffoldBank(scaffoldsFilename);

  int contigs = 0;
  int gaps = 0;
  int scaffoldIndex = 0;

  contigIter.first();
  std::string scaffold = "";
  std::string scaffoldComment = contigIter->getComment();

  for (contigIter.first(); !contigIter.isDone(); contigIter.next()) {
    const std::string comment = contigIter->getComment();
    const std::string contig = contigIter->toString();

    const int contigScaffoldIndex = parseScaffoldIndex(comment);
    const int contigIndex = parseContigIndex(comment);
    const int gapIndex = parseGapIndex(comment);

    assert(contigIndex == contigs);
    contigs++;

    // Scaffold done, write to file
    if (contigScaffoldIndex != scaffoldIndex) {
      assert(contigScaffoldIndex == scaffoldIndex + 1);
      insertSequence(scaffoldBank, scaffoldComment, scaffold);

      scaffold = "";
      scaffoldComment = comment;
      scaffoldIndex = contigScaffoldIndex;
    }

    scaffold += contig;

    // Find the gap corresponding to the contig and append to scaffold
    if (gapIndex != -1) {
      std::string first = "", second = "";
      // TODO: Speed this up, takes quadratic time
      for (gapIter.first(); !gapIter.isDone(); gapIter.next()) {
        if (gapIndex == parseGapIndex(gapIter->getComment())) {
          const size_t markerPos = gapIter->getComment().find(STR_SPLIT_MARKER);
          if (markerPos == std::string::npos) {
            first = gapIter->toString();
            break;
          } else {
            // Combine split gaps
            const int splitNum = std::stoi(gapIter->getComment().substr(markerPos + SPLIT_MARKER_LENGTH, 1));
            if (splitNum == 1) {
              assert(first == "");
              first = gapIter->toString();
            } else {
              assert(second == "");
              // Remove flank from second gap
              // TODO: Fix case where flank is changed by ignored nucleotides
              const int splitLength = std::stoi(gapIter->getComment().substr(markerPos + SPLIT_MARKER_LENGTH + 2));
              second = gapIter->toString().substr(splitLength);
            }

#ifdef DEBUG
            std::cout << "first: " << first << " second: " << second << std::endl;
#endif

            if (first != "" && second != "") {
              break;
            }
          }
        }
      }

      scaffold += first + second;
      assert(gapIndex == gaps);
      gaps++;
    }
  }

  if (scaffold.size() > 0) {
    insertSequence(scaffoldBank, scaffoldComment, scaffold);
    scaffoldIndex++;
  }

  scaffoldBank.flush();

  std::cout << "Merged " << contigs << " contigs and " << gaps << " gaps into " << scaffoldIndex << " scaffolds" << std::endl;
}

/****************************************************************************/

int main(int argc, char *argv[])
{
  try {
    GapMerger().run(argc, argv);
  } catch (Exception &e) {
    std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
