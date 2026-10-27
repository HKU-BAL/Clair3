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

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

#[derive(Clone, Copy, Default)]
struct Edge {
    weight: usize,
    is_reference: bool,
}

struct Graph {
    kmers: Vec<Vec<u8>>,
    forward: Vec<BTreeMap<usize, Edge>>,
    reverse: Vec<BTreeSet<usize>>,
    active: Vec<bool>,
    by_kmer: HashMap<Vec<u8>, usize>,
    source: usize,
    sink: usize,
    k: usize,
}

impl Graph {
    fn new(reference: &[u8], reads: &[&[u8]], low_quality: &[BTreeSet<usize>], k: usize) -> Self {
        let mut graph = Self {
            kmers: Vec::new(),
            forward: Vec::new(),
            reverse: Vec::new(),
            active: Vec::new(),
            by_kmer: HashMap::new(),
            source: 0,
            sink: 0,
            k,
        };
        graph.add_sequence(reference, 0, reference.len() - k, true);
        graph.source = graph.vertex(&reference[..k]);
        graph.sink = graph.vertex(&reference[reference.len() - k..]);
        for (read_index, read) in reads.iter().enumerate() {
            let empty = BTreeSet::new();
            let low_quality = low_quality.get(read_index).unwrap_or(&empty);
            let stop = read.len() as isize - k as isize;
            let mut start = 0_usize;
            while (start as isize) < stop {
                let mut next_bad = read.len();
                for (index, base) in read.iter().enumerate().skip(start) {
                    if !matches!(base, b'A' | b'C' | b'G' | b'T') || low_quality.contains(&index) {
                        next_bad = index;
                        break;
                    }
                }
                if next_bad >= k {
                    graph.add_sequence(read, start, next_bad - k, false);
                }
                start = next_bad + 1;
            }
        }
        graph
    }

    fn vertex(&mut self, kmer: &[u8]) -> usize {
        if let Some(index) = self.by_kmer.get(kmer) {
            return *index;
        }
        let index = self.kmers.len();
        let kmer = kmer.to_vec();
        self.kmers.push(kmer.clone());
        self.forward.push(BTreeMap::new());
        self.reverse.push(BTreeSet::new());
        self.active.push(true);
        self.by_kmer.insert(kmer, index);
        index
    }

    fn add_sequence(&mut self, sequence: &[u8], start: usize, end: usize, is_reference: bool) {
        if end == 0 {
            return;
        }
        let mut previous = self.vertex(&sequence[start..start + self.k]);
        for index in start + 1..=end {
            let current = self.vertex(&sequence[index..index + self.k]);
            let edge = self.forward[previous].entry(current).or_default();
            edge.weight += 1;
            edge.is_reference |= is_reference;
            self.reverse[current].insert(previous);
            previous = current;
        }
    }

    fn has_cycle(&self) -> bool {
        fn visit(vertex: usize, graph: &Graph, colors: &mut [u8]) -> bool {
            if colors[vertex] == 1 {
                return true;
            }
            if colors[vertex] == 2 {
                return false;
            }
            colors[vertex] = 1;
            for next in graph.forward[vertex].keys() {
                if visit(*next, graph, colors) {
                    return true;
                }
            }
            colors[vertex] = 2;
            false
        }
        let mut colors = vec![0; self.kmers.len()];
        (0..self.kmers.len()).any(|vertex| visit(vertex, self, &mut colors))
    }

    fn prune(&mut self) {
        for edges in &mut self.forward {
            edges.retain(|_, edge| edge.is_reference || edge.weight >= 2);
        }
        for incoming in &mut self.reverse {
            incoming.clear();
        }
        for (source, edges) in self.forward.iter().enumerate() {
            for target in edges.keys() {
                self.reverse[*target].insert(source);
            }
        }
        let forward = self.reachable(self.source, false);
        let reverse = self.reachable(self.sink, true);
        for vertex in 0..self.kmers.len() {
            self.active[vertex] = forward[vertex] && reverse[vertex];
        }
        for vertex in 0..self.kmers.len() {
            if self.active[vertex] {
                self.forward[vertex].retain(|next, _| self.active[*next]);
                self.reverse[vertex].retain(|previous| self.active[*previous]);
            } else {
                self.forward[vertex].clear();
                self.reverse[vertex].clear();
            }
        }
    }

    fn reachable(&self, start: usize, reverse: bool) -> Vec<bool> {
        let mut reachable = vec![false; self.kmers.len()];
        let mut pending = vec![start];
        reachable[start] = true;
        while let Some(vertex) = pending.pop() {
            let next: Vec<usize> = if reverse {
                self.reverse[vertex].iter().copied().collect()
            } else {
                self.forward[vertex].keys().copied().collect()
            };
            for next in next {
                if !reachable[next] {
                    reachable[next] = true;
                    pending.push(next);
                }
            }
        }
        reachable
    }

    fn candidate_haplotypes(&self) -> Vec<Vec<u8>> {
        let mut paths = Vec::new();
        let mut pending = VecDeque::from([vec![self.source]]);
        while let Some(path) = pending.pop_front() {
            if paths.len() + pending.len() > 256 {
                return Vec::new();
            }
            let last = *path.last().unwrap();
            for next in self.forward[last].keys() {
                let mut extended = path.clone();
                extended.push(*next);
                if *next == self.sink || self.forward[*next].is_empty() {
                    paths.push(extended);
                } else {
                    pending.push_back(extended);
                }
            }
        }
        let mut haplotypes = Vec::with_capacity(paths.len());
        for path in paths {
            let mut haplotype = Vec::new();
            for vertex in &path {
                haplotype.push(self.kmers[*vertex][0]);
            }
            let last = *path.last().unwrap();
            haplotype.extend_from_slice(&self.kmers[last][1..]);
            haplotypes.push(haplotype);
        }
        haplotypes.sort();
        haplotypes
    }
}

pub fn consensus(
    reference: &[u8],
    reads: &[&[u8]],
    low_quality: &[BTreeSet<usize>],
) -> Vec<Vec<u8>> {
    let max_k = 101.min(reference.len().saturating_sub(1));
    let min_k = (10..=max_k).find(|k| {
        let mut kmers = BTreeSet::new();
        (0..=reference.len() - k).all(|index| kmers.insert(&reference[index..index + k]))
    });
    let Some(min_k) = min_k else {
        return Vec::new();
    };
    for k in min_k..=max_k {
        let mut graph = Graph::new(reference, reads, low_quality, k);
        if !graph.has_cycle() {
            graph.prune();
            return graph.candidate_haplotypes();
        }
    }
    Vec::new()
}
