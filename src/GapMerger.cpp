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

/****************************************************************************/

// Constants for command line parameters
static const char* STR_SCAFFOLDS = "-scaffolds";
static const char* STR_CONTIGS = "-contigs";
static const char* STR_GAPS = "-gaps";

/****************************************************************************/

static const char* STR_SCAFFOLD_MARKER = " scaffold ";
static const char* STR_GAP_MARKER = " gap ";
static const char* STR_SPLIT_MARKER = " split ";

static const size_t SCAFFOLD_MARKER_LENGTH = 10u;
static const size_t GAP_MARKER_LENGTH = 5u;
static const size_t SPLIT_MARKER_LENGTH = 7u;

/****************************************************************************/

class GapMerger : public Tool
{
public:

    // Constructor
    GapMerger ();

    // Actual job done by the tool is here
    void execute ();

    int parseGapIndex(const std::string &) const;
    int parseScaffoldIndex(const std::string &) const;

    void insertSequence(BankFasta &, const std::string &,
      const std::string &) const;
};


/*****************************************************************************
Constructor for the tool.
*****************************************************************************/

GapMerger::GapMerger ()  : Tool ("GapMerger")
{
    // We add some custom arguments for command line interface
    getParser()->push_front (new OptionOneParam (STR_SCAFFOLDS, "FASTA file of merged scaffolds",  true));
    getParser()->push_front (new OptionOneParam (STR_CONTIGS, "FASTA file of contigs",  true));
    getParser()->push_front (new OptionOneParam (STR_GAPS, "FASTA file of filled gaps",  true));
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

int GapMerger::parseScaffoldIndex(const std::string &comment) const
{
  const size_t scaffold_pos = comment.find(STR_SCAFFOLD_MARKER);
  if (scaffold_pos == std::string::npos) {
    return -1;
  }

  std::string index = "";
  const size_t gap_pos = comment.find(STR_GAP_MARKER);
  if (gap_pos == std::string::npos) {
    index = comment.substr(scaffold_pos + SCAFFOLD_MARKER_LENGTH);
  } else {
    index = comment.substr(scaffold_pos + SCAFFOLD_MARKER_LENGTH,
      gap_pos - (scaffold_pos + SCAFFOLD_MARKER_LENGTH));
  }

  return std::stoi(index);
}

void GapMerger::insertSequence(BankFasta &bank, const std::string &comment,
  const std::string &sequence) const
{
  Sequence seq((char *) sequence.c_str());
  seq._comment = comment;
  bank.insert(seq);
}

void GapMerger::execute ()
{
  // Get the command line arguments
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
  std::string scaffold = "";
  std::string scaffoldComment = "";

  for (contigIter.first(); !contigIter.isDone(); contigIter.next()) {
    contigs++;

    const std::string comment = contigIter->getComment();
    const std::string contig = contigIter->toString();

    const int contigScaffoldIndex = parseScaffoldIndex(comment);
    const int gapIndex = parseGapIndex(comment);

    // Scaffold done, write to file
    if (contigScaffoldIndex != scaffoldIndex) {
      if (scaffold.size() > 0) {
        insertSequence(scaffoldBank, scaffoldComment, scaffold);
      }

      scaffold = "";
      scaffoldComment = "";
      scaffoldIndex = contigScaffoldIndex;
    }

    // Remove markers from comment, i.e. remove anything after scaffold marker
    if (scaffoldComment.size() == 0) {
      const size_t marker_pos = comment.find(STR_SCAFFOLD_MARKER);

      if (marker_pos == std::string::npos) {
        scaffoldComment = comment;
      } else {
        scaffoldComment = comment.substr(0, marker_pos);
      }
    }

    scaffold += contig;

    if (gapIndex != -1) {
      // Find the gap corresponding to the contig and append to scaffold
      std::string first = "", second = "";
      for (gapIter.first(); !gapIter.isDone(); gapIter.next()) {
        if (gapIndex == parseGapIndex(gapIter->getComment())) {
          gaps++;

          // Combine split gaps
          // TODO: Refactor this
          const size_t marker_pos = gapIter->getComment().find(STR_SPLIT_MARKER);
          if (marker_pos == std::string::npos) {
            first = gapIter->toString();
            break;
          } else {
            const int split_num = std::stoi(gapIter->getComment().substr(marker_pos + SPLIT_MARKER_LENGTH, 1));
            if (split_num == 1) {
              first = gapIter->toString();
            } else {
              // Remove flank from second gap
              // TODO: Fix case where flank is changed by ignored nucleotides
              const int split_length = std::stoi(gapIter->getComment().substr(marker_pos + SPLIT_MARKER_LENGTH + 2));
              second = gapIter->toString().substr(split_length);
            }

            std::cout << "first: " << first << " second: " << second << std::endl;

            if (first != "" && second != "") {
              break;
            }
          }
        }
      }

      scaffold += first + second;
    }
  }

  if (scaffold.size() > 0) {
    insertSequence(scaffoldBank, scaffoldComment, scaffold);
    scaffoldIndex++;
  }

  scaffoldBank.flush();

  std::cout << "Merged " << contigs << " contigs and " << gaps <<
    " gaps into " << scaffoldIndex <<  " scaffolds" << std::endl;
}

/****************************************************************************/

int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.
        GapMerger().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
