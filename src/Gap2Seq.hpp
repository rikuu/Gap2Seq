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

#ifndef _TOOL_Gap2Seq_HPP_
#define _TOOL_Gap2Seq_HPP_

#include <vector>
#include <string>
#include <unordered_map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>

#include <gatb/gatb_core.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS > digraph;
typedef boost::graph_traits<digraph>::vertex_descriptor bnode;
typedef boost::graph_traits<digraph>::edge_descriptor bedge;
typedef std::vector<bnode> container;

// Check if a file is readable
inline bool is_readable(const std::string & file) {
    std::ifstream f( file.c_str() );
    return !f.fail();
}

#ifndef SINGLE_THREAD
// Synchronizer for access to global variables
extern ISynchronizer *global_lock;
#endif

// Array of memory usage accessed by thread id
extern std::unordered_map<long long, long long> memuse;

// Subgraph statistics
struct subgraph_stats {
  size_t vertices;
  size_t edges;
  size_t nontrivial_components;
  size_t size_nontrivial_components;
  size_t vertices_final;
  size_t edges_final;
};

class Gap2Seq : public Tool
{
public:

    // Constructor
    Gap2Seq ();

    // Actual job done by the tool is here
    void execute ();

    // Extract nodes reachable from a given kmer with given distance
  std::set<Node> extract_reachable_nodes(const Graph &graph, const std::string &kmer, int d);

  // Fill a gap
  int fill_gap(const Graph &graph, const std::string &kmer_left, const std::string &kmer_right, int gap_len,
      int k, int gap_err, int left_max_fuz, int right_max_fuz, int *left_fuz, int *right_fuz,
      long long max_mem, char *fill, bool skip_confident, struct subgraph_stats *substats);
};

/********************************************************************************/

#endif /* _TOOL_Gap2Seq_HPP_ */
