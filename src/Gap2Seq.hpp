/*****************************************************************************
 *  Gap2Seq
 *  Copyright (C) Leena Salmela, Kristoffer Sahlin, Veli Mäkinen,
 *  Alexandru Tomescu 2014
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

#ifndef _TOOL_Gap2Seq_HPP_
#define _TOOL_Gap2Seq_HPP_

#include <vector>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

#include <gatb/gatb_core.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS > digraph;
typedef boost::graph_traits<digraph>::vertex_descriptor bnode;
typedef std::vector<bnode> container;

// Check if a file is readable
inline bool is_readable(const std::string & file) {
    std::ifstream f( file.c_str() );
    return !f.fail();
}

// Subgraph statistics
struct subgraph_stats {
  int vertices;
  int edges;
  int nontrivial_components;
  int size_nontrivial_components;
  int vertices_final;
  int edges_final;
};

class Gap2Seq : public Tool
{
public:

    // Constructor
    Gap2Seq ();

    // Actual job done by the tool is here
    void execute ();

    // Extract nodes reachable from a given kmer with given distance
  std::set<Node> extract_reachable_nodes(Graph graph, std::string kmer, int d);

  // Fill a gap
  int fill_gap(Graph graph, std::string kmer_left, std::string kmer_right, int c, int d, int fuz, int k, long long max_mem, int *left_fuz, int *right_fuz, char *fill, bool skip_confident, struct subgraph_stats *substats);

};

/********************************************************************************/

#endif /* _TOOL_Gap2Seq_HPP_ */
