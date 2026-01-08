/*
Copyright 2020 Google LLC.
Copyright 2021 The University of Hong Kong, Department of Computer Science

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <unordered_map>
#include <memory>
#include <vector>
#include <string>
#include<map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

using namespace std;

struct struct_str_arr
{
    int consensus_size;
	char* consensus[500];
};


struct VertexInfo {
  string kmer;
};


struct EdgeInfo {
  int weight;
  bool is_ref;
};


class DeBruijnGraph {
 public:
  using BoostGraph = boost::adjacency_list<
    boost::setS,
    boost::listS,
    boost::bidirectionalS,
    VertexInfo,
    EdgeInfo>;

  using VertexIterator = boost::graph_traits<BoostGraph>::vertex_iterator;
  using EdgeIterator = boost::graph_traits<BoostGraph>::edge_iterator;
  using AdjacencyIterator = boost::graph_traits<BoostGraph>::adjacency_iterator;

 public:
  using Vertex = boost::graph_traits<BoostGraph>::vertex_descriptor;
  using Edge = boost::graph_traits<BoostGraph>::edge_descriptor;
  using Path = vector<Vertex>;

  using RawVertexIndexMap = map<Vertex, int>;
  using VertexIndexMap =
      boost::const_associative_property_map<RawVertexIndexMap>;


 public:

  void RebuildIndexMap();

  VertexIndexMap IndexMap() const;

  Vertex EnsureVertex(string kmer);

  Vertex VertexForKmer(string kmer) const;

  bool HasCycle() const;

  DeBruijnGraph(
      const string& ref,
      const vector<string>& reads,
      vector<set<int> >& base_quality,
      int k);

  Edge AddEdge(Vertex from_vertex, Vertex to_vertex, bool is_ref);


  void AddKmersAndEdges(string bases, int start, int end,
                        bool is_ref);

  void AddEdgesForReference(string ref);

  void AddEdgesForRead(const string& read, set<int>& base_quality_set);

  vector<Path> CandidatePaths() const;

  string HaplotypeForPath(const Path& path) const;

  void Prune();

 public:

  static vector<string> Build(
       const string& ref,

       const vector<string>& reads,
       vector<set<int> >& base_quality);

  vector<string> CandidateHaplotypes() const;

  string GraphViz() const;

  int KmerSize() const { return k_; }

 public:
  BoostGraph g_;
  int k_;
  Vertex source_;
  Vertex sink_;

  unordered_map<string, Vertex> kmer_to_vertex_;
  RawVertexIndexMap vertex_index_map_;
};

