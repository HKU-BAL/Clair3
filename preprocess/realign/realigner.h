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
#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include<string>
#include <unordered_map>
#include "ssw_cpp.h"


using namespace std;

struct struct_str_arr
{
	int position[1000];
	char* cigar_string[1000];
};

enum Operation {
    OPERATION_UNSPECIFIED = 0,
    ALIGNMENT_MATCH = 1,
    INSERT = 2,
    DELETE = 3,
    SKIP = 4,
    CLIP_SOFT = 5,
    CLIP_HARD = 6,
};

struct Kmer {
  string sequence;
};

struct ReadId {
  ReadId() : is_set(false), id(0) {}
  explicit ReadId(int id) : is_set(true), id(id) {}

  explicit operator int64_t() const { return id; }
  explicit operator uint64_t() const { return id; }

  bool operator<(const ReadId& that) const { return id < that.id; }
  bool operator==(const ReadId& that) const { return id == that.id; }

  bool is_set;
  int id;
};


struct KmerOffset {
  KmerOffset() : is_set(false), pos(0) {}
  explicit KmerOffset(int pos) : is_set(true), pos(pos) {}

  bool operator==(const KmerOffset& that) const {
    return pos == that.pos && is_set == that.is_set;
  }

  bool is_set;
  int pos;
};

struct KmerOccurrence {
  KmerOccurrence() {}
  KmerOccurrence(ReadId read_id, KmerOffset pos)
      : read_id(read_id), read_pos(pos) {}

  bool operator==(const KmerOccurrence& that) const {
    return read_id == that.read_id && read_pos == that.read_pos;
  }

  ReadId read_id;
  KmerOffset read_pos;
};


struct ReadAlignment {
  static const int kNotAligned = -1;
      ReadAlignment() : position(kNotAligned), cigar(""), score(0) {}

  ReadAlignment(int position_param, const string& cigar_param,
                int score_param)
      : position(position_param), cigar(cigar_param), score(score_param) {}

  bool operator==(const ReadAlignment& that) const {
    return score == that.score && position == that.position &&
           cigar == that.cigar;
  }

  void reset() {
    score = 0;
    position = kNotAligned;
    cigar = "";
  }

  int position;
  string cigar;
  int score;
};

struct CigarOp {
  CigarOp()
      : operation(OPERATION_UNSPECIFIED),
        length(0) {}
  CigarOp(Operation op, int len) : operation(op), length(len) {}

  bool operator==(const CigarOp& that) const {
    return operation == that.operation && length == that.length;
  }

  Operation operation;
  int length;
};


struct Read {
    int position;
    list<CigarOp> cigar;
    string cigar_string;
    string seq;

    Read(int position, string cigar_string, string seq)
	{
		this->position = position;
		this->cigar_string = cigar_string;
		this->seq = seq;
	};

};


struct HaplotypeReadsAlignment {
  HaplotypeReadsAlignment() : haplotype_index(0), haplotype_score(0) {}
  HaplotypeReadsAlignment(
      int haplotype_index,
      int score,
      const vector<ReadAlignment>& read_alignment_scores)
        : haplotype_index(haplotype_index), haplotype_score(score) {
    this->read_alignment_scores.assign(read_alignment_scores.begin(),
                                       read_alignment_scores.end());
  }

  bool operator==(const HaplotypeReadsAlignment& that) const {
    return haplotype_index == that.haplotype_index &&
           haplotype_score == that.haplotype_score &&
           read_alignment_scores == that.read_alignment_scores &&
           cigar == that.cigar && cigar_ops == that.cigar_ops &&
           is_reference == that.is_reference &&
           hap_to_ref_positions_map == that.hap_to_ref_positions_map;
  }

  bool operator<(const HaplotypeReadsAlignment& that) const {
    return haplotype_score < that.haplotype_score;
  }

  int haplotype_index;

  int haplotype_score;

  vector<ReadAlignment> read_alignment_scores;

  string cigar;


  list<CigarOp> cigar_ops;

  int ref_pos;


  vector<int> hap_to_ref_positions_map;


  bool is_reference;
};


void SetPositionsMap(int haplotype_size,
                     HaplotypeReadsAlignment* hyplotype_alignment);

void MergeCigarOp(const CigarOp& op, int read_len, list<CigarOp>* cigar);

using KmerIndexType =
    unordered_map<string, vector<KmerOccurrence> >;

class ReAligner {
 public:

    struct_str_arr* realign_reads(char* seqs[], int* positions, char* cigars[], char* reference, char* haplotypes,  int ref_start, int ref_prefix, int ref_suffix, int read_size);

  void set_reference(const string& reference);
  void set_reads(const vector<string>& reads);
  vector<string> get_reads() const { return reads_; }
  void set_ref_start(uint64_t position);
  void set_haplotypes(const vector<string>& haplotypes);
  uint8_t get_match_score() const { return match_score_; }
  uint8_t get_mismatch_penalty() const { return mismatch_penalty_; }
  void set_options();
  void set_ref_prefix_len(int ref_prefix_len) {
    ref_prefix_len_ = ref_prefix_len;
  }
  void set_ref_suffix_len(int ref_suffix_len) {
    ref_suffix_len_ = ref_suffix_len;
  }
  int get_ssw_alignment_score_threshold() const {
    return ssw_alignment_score_threshold_;
  }

vector<struct Read> AlignReads(
    const vector<struct Read>& reads_param);

  void BuildIndex();
//
  KmerIndexType GetKmerIndex() const { return kmer_index_; }

  void FastAlignReadsToHaplotype(
      const string& haplotype, int* haplotype_score,
      vector<ReadAlignment>* haplotype_read_alignment_scores);
  void SswAlignReadsToHaplotypes(int score_threshold);

  void InitSswLib();

  void SswSetReference(const string& reference);
//
  StripedSmithWaterman::Alignment SswAlign(const string& target) const;
//
  void AlignHaplotypesToReference();
//
  const vector<HaplotypeReadsAlignment>& GetReadToHaplotypeAlignments()
      const {
    return read_to_haplotype_alignments_;
  }

  bool GetBestReadAlignment(int readId, int* best_hap_index) const;
  void CalculateReadToRefAlignment(
      int read_index,
      const ReadAlignment& read_to_haplotype_alignment,
      const list<CigarOp>& haplotype_to_ref_cigar_ops_input,
      list<CigarOp>* read_to_ref_cigar_ops) const;

  void RealignReadsToReference(
      const vector<struct Read>& reads,
    vector<struct Read>& realigned_reads);

  void CalculateSswAlignmentScoreThreshold();

// private:

  string reference_;

  string region_chromosome_ ;

  int region_position_in_chr_;

  vector<string> haplotypes_;

  vector<HaplotypeReadsAlignment> read_to_haplotype_alignments_;

  KmerIndexType kmer_index_;

  vector<string> reads_;

  int kmer_size_ = 32;

  int read_size_ = 250;

  int max_num_of_mismatches_ = 2;

  int ssw_alignment_score_threshold_ = 1;

  int match_score_ = 4;
  int mismatch_penalty_ = 6;
  int gap_opening_penalty_ = 8;
  int gap_extending_penalty_ = 1;

  double similarity_threshold_ = 0.85;

//  unique_ptr<Aligner> ssw_aligner_;
  StripedSmithWaterman::Aligner ssw_aligner_;
//  StripedSmithWaterman::Filter filter;
//  StripedSmithWaterman::Alignment alignment;

  int ref_prefix_len_;
  int ref_suffix_len_;

  void FastAlignReadsToHaplotypes();

  void AddReadToIndex(const string& read, ReadId read_id);

  void AddKmerToIndex(string kmer, ReadId read_id,
                      KmerOffset pos);

  int FastAlignStrings(string s1, string s2,
                       int max_mismatches, int* num_of_mismatches) const;

  void UpdateBestHaplotypes(
      int haplotype_index, int haplotype_score,
      const vector<ReadAlignment>& current_read_scores);

  void CalculatePositionMaps();
};
