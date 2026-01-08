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

#include "debruijn_graph.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <unordered_set>
#include <boost/algorithm/string.hpp>
#include <set>

using namespace std;

using Vertex = DeBruijnGraph::Vertex;
using VertexIndexMap = DeBruijnGraph::VertexIndexMap;
using Edge = DeBruijnGraph::Edge;
using Path = DeBruijnGraph::Path;


class CycleDetector : public boost::dfs_visitor<> {
public:
    explicit CycleDetector(bool* has_cycle) : has_cycle(has_cycle) {}

template <class Edge, class Graph>
void back_edge(Edge, const Graph&) {
    *has_cycle = true;
}

private:
    bool* has_cycle;
};

template <class BoostGraph>
class EdgeLabelWriter {
public:
    explicit EdgeLabelWriter(const BoostGraph& g) : g_(g) {}

void operator()(ostream& out, const Edge e) const {
     EdgeInfo ei = g_[e];
     out << "[label=" << ei.weight << (ei.is_ref ? " color=red" : "") << "]";
}

private:
    const BoostGraph& g_;
};

class ReachableVertexVisitor : public boost::dfs_visitor<> {
public:
    explicit ReachableVertexVisitor(set<Vertex>* reachable_vertices)
       : reachable_vertices(reachable_vertices) {}

template <class Edge, class Graph>
void tree_edge(Edge e, const Graph& g) {
 Vertex from = boost::source(e, g);
 if (reachable_vertices->find(from) != reachable_vertices->end()) {
   Vertex to = boost::target(e, g);
   reachable_vertices->insert(to);
 }
}

private:
    set<Vertex>* reachable_vertices;
};

template <class BoostGraphT, class VertexIndexMapT>
set<Vertex> VerticesReachableFrom(
 Vertex v, const BoostGraphT& g, const VertexIndexMapT& vertex_index_map) {
    set<Vertex> reachable_vertices{v};
    ReachableVertexVisitor vis(&reachable_vertices);
    boost::depth_first_search(
   g, boost::visitor(vis).root_vertex(v).vertex_index_map(vertex_index_map));
    return reachable_vertices;
}


Vertex DeBruijnGraph::EnsureVertex(string kmer) {
    Vertex v;
    auto vertex_find = kmer_to_vertex_.find(kmer);
    if (vertex_find != kmer_to_vertex_.end()) {
     v = (*vertex_find).second;
    } else {
     string kmer_copy(kmer);
     v = boost::add_vertex(VertexInfo{kmer_copy}, g_);
     kmer_to_vertex_[string(g_[v].kmer)] = v;
    }

    return v;
}

Vertex DeBruijnGraph::VertexForKmer(string kmer) const {
    return kmer_to_vertex_.at(kmer);
}

void DeBruijnGraph::RebuildIndexMap() {
    map<Vertex, int> table;
    VertexIterator vi, vend;
    tie(vi, vend) = boost::vertices(g_);
    int index = 0;
    for (; vi != vend; ++vi) {
     table[*vi] = index;
     ++index;
    }
    vertex_index_map_ = table;
}

VertexIndexMap DeBruijnGraph::IndexMap() const {
    boost::const_associative_property_map<RawVertexIndexMap> vmap(
       vertex_index_map_);
    return vmap;
}

bool DeBruijnGraph::HasCycle() const {
    bool has_cycle = false;
    CycleDetector cycle_detector(&has_cycle);
    boost::depth_first_search(
       g_, boost::visitor(cycle_detector).vertex_index_map(IndexMap()));
    return has_cycle;
}

DeBruijnGraph::DeBruijnGraph(
 const string& ref,
 const vector<string>& reads,
 vector<set<int> >& base_quality,
 int k)
 : k_(k) {

    AddEdgesForReference(ref);
    source_ = VertexForKmer(ref.substr(0, k_));
    sink_ = VertexForKmer(ref.substr(ref.size() - k_, k_));
    for (int i = 0; i < reads.size(); i++) {
      AddEdgesForRead(reads[i], base_quality[i]);
    //  }
    }
    RebuildIndexMap();
}


constexpr int kBoundsNoWorkingK = -1;
struct KBounds {
    int min_k;
    int max_k;
};

KBounds KMinMaxFromReference(const string ref) {
KBounds bounds;
bounds.min_k = kBoundsNoWorkingK;
bounds.max_k  = min(101, static_cast<int>(ref.size()) - 1);

for (int k = 10; k <= bounds.max_k; k++) {
     bool has_cycle = false;
     set<string> kmers;

     for (int i = 0; i < ref.size() - k + 1; i++) {
       string kmer = ref.substr(i, k);
       if (kmers.insert(kmer).second == false) {
         has_cycle = true;
         break;
       }
     }

     if (!has_cycle) {
       bounds.min_k = k;
       break;
     }
    }

    return bounds;
}


vector<string> DeBruijnGraph::Build(
  const string& ref,
  const vector<string>& reads,
  vector<set<int> >& base_quality) {
vector<string> haplotypes;
KBounds bounds = KMinMaxFromReference(ref);
  if (bounds.min_k == kBoundsNoWorkingK) return haplotypes;

for (int k = bounds.min_k; k <= bounds.max_k; k++) {
 shared_ptr<DeBruijnGraph> graph = shared_ptr<DeBruijnGraph>(
     new DeBruijnGraph(ref, reads, base_quality, k));
 if (graph->HasCycle()) {
   continue;
 } else {
   graph->Prune();

    for (const Path& path : graph->CandidatePaths()) {
     haplotypes.push_back(graph->HaplotypeForPath(path));
    }
    sort(haplotypes.begin(), haplotypes.end());
    return haplotypes;
     }
}
  return haplotypes;
}

Edge DeBruijnGraph::AddEdge(Vertex from_vertex, Vertex to_vertex, bool is_ref) {
    bool was_present;
    Edge edge;
    tie(edge, was_present) = boost::edge(from_vertex, to_vertex, g_);
    if (!was_present) {
     tie(edge, ignore) = boost::add_edge(from_vertex, to_vertex,
                                                   EdgeInfo{0, false}, g_);
    }
    EdgeInfo& ei = g_[edge];
    ei.weight++;
    ei.is_ref |= is_ref;
    return edge;
}

void DeBruijnGraph::AddKmersAndEdges(string bases, int start, int end,
                                  bool is_ref) {
    if (end > 0) {
     Vertex vertex_prev = EnsureVertex(bases.substr(start, k_));
     for (int i = start + 1; i <= end; ++i) {
       Vertex vertex_cur = EnsureVertex(bases.substr(i, k_));
       AddEdge(vertex_prev, vertex_cur, is_ref);
       vertex_prev = vertex_cur;
     }
    }
}

void DeBruijnGraph::AddEdgesForReference(string ref) {
    AddKmersAndEdges(ref, 0, ref.size() - k_, true);
}


void DeBruijnGraph::AddEdgesForRead(const string& read, set<int>& base_quality_set) {
    const string bases = read;

    auto NextBadPosition = [&read, &bases, &base_quality_set, this](int start) -> int {
      string ACGT = "ACGT";

      for (int i = start; i < bases.size(); ++i) {
         if (ACGT.find(bases[i]) == string::npos || base_quality_set.find(i) != base_quality_set.end()) {
          return i;
      }
      }
     return bases.size();
};

const string bases_view(bases);
    const int stop = bases.size() - k_;
    int i = 0;
    while (i < stop) {
     int next_bad_position = NextBadPosition(i);
     AddKmersAndEdges(bases_view, i, next_bad_position - k_, false /* is_ref */);
     i = next_bad_position + 1;
    }
}

vector<Path> DeBruijnGraph::CandidatePaths() const {
vector<Path> terminated_paths;
queue<Path> extendable_paths;

extendable_paths.push({source_});

while (!extendable_paths.empty()) {
     int n_total_paths = terminated_paths.size() + extendable_paths.size();
     if (n_total_paths > 256) {
       return {};
     }

     Path path = extendable_paths.front();
     extendable_paths.pop();
     Vertex last_v = path.back();
     AdjacencyIterator vi, vend;
     tie(vi, vend) = boost::adjacent_vertices(last_v, g_);
     for (; vi != vend; ++vi) {
       Path extended_path(path);
       extended_path.push_back(*vi);
       if (*vi == sink_ || boost::out_degree(*vi, g_) == 0) {
         terminated_paths.push_back(extended_path);
       } else {
         extendable_paths.push(extended_path);
       }
     }
    }
    return terminated_paths;
}

string DeBruijnGraph::HaplotypeForPath(const Path& path) const {
    stringstream haplotype;
    for (Vertex v : path) {
     haplotype << g_[v].kmer[0];
    }
    if (!path.empty()) {
     haplotype << g_[path.back()].kmer.substr(1, k_ - 1);
    }
    return haplotype.str();
}

vector<string> DeBruijnGraph::CandidateHaplotypes() const {
    vector<string> haplotypes;
    for (const Path& path : CandidatePaths()) {
     haplotypes.push_back(HaplotypeForPath(path));
    }
    sort(haplotypes.begin(), haplotypes.end());
    return haplotypes;
}

string DeBruijnGraph::GraphViz() const {
    stringstream graphviz;
    auto vertex_label_writer = boost::make_label_writer(
       boost::get(&VertexInfo::kmer, g_));
    boost::write_graphviz(
       graphviz,
       g_,
       vertex_label_writer,
       EdgeLabelWriter<BoostGraph>(g_),
       boost::default_writer(),
       IndexMap());
    return graphviz.str();
}

void DeBruijnGraph::Prune() {
    boost::remove_edge_if(
       [this](const Edge& e) {
         return !g_[e].is_ref && g_[e].weight < 2;
       },
       g_);

    // Remove vertices not reachable forward from src or backward from sink.
    VertexIterator vbegin, vend;
    tie(vbegin, vend) = boost::vertices(g_);
    set<Vertex> all_vertices(vbegin, vend);

    set<Vertex> fwd_reachable_vertices, rev_reachable_vertices;
    fwd_reachable_vertices = VerticesReachableFrom(
       source_, g_, IndexMap());
    rev_reachable_vertices = VerticesReachableFrom(
       sink_, boost::make_reverse_graph(g_), IndexMap());

    set<Vertex> reachable_vertices;
    set_intersection(
       fwd_reachable_vertices.begin(), fwd_reachable_vertices.end(),
       rev_reachable_vertices.begin(), rev_reachable_vertices.end(),
       inserter(reachable_vertices, reachable_vertices.end()));
    for (Vertex v : all_vertices) {
     if (reachable_vertices.find(v) == reachable_vertices.end()) {
       kmer_to_vertex_.erase(g_[v].kmer);
       boost::clear_vertex(v, g_);
       boost::remove_vertex(v, g_);
     }
    }
    RebuildIndexMap();
}


extern "C" {
    struct_str_arr* get_consensus(char* reference, char* c_reads, char* c_base_quality, int read_size) {
    const string ref = reference;
    vector<string> reads;
    string r = c_reads;
    string r2 = c_base_quality;
    string item;
    boost::split(reads, c_reads, boost::is_any_of(","));
    vector<set<int> > base_quality_set;
    vector<string> base_quality_array; //( c_base_quality, c_base_quality + read_size);


   string delimiter = " ";
   boost::split(base_quality_array, c_base_quality, boost::is_any_of(","));

   for (int i = 0; i < base_quality_array.size(); i++) {
        string s = base_quality_array[i];

        size_t pos = 0;
        set<int> bq_pos;
        string token;
        stringstream ss(s);
        int temp;
        while (ss >> temp) {
            bq_pos.insert(temp);
            }

        base_quality_set.push_back(bq_pos);

   }
    vector<string> output = DeBruijnGraph::Build(ref, reads, base_quality_set);
    struct_str_arr* str_arr_ptr = new struct_str_arr();

    str_arr_ptr->consensus_size = output.size();
    for (int i=0; i < output.size(); i++) {
        str_arr_ptr->consensus[i] = new char[output[i].size() + 1];
        strcpy(str_arr_ptr->consensus[i], output[i].c_str());
    }
    return str_arr_ptr;

}
}

extern "C" {
    void free_memory(struct_str_arr* pointer, int size) {
        for (int i=0; i< size; i++) {
            delete [] pointer->consensus[i];
        }
        delete pointer;
        pointer = NULL;
    }
}
