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

#include "realigner.h"
#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <string>
#include "ssw_cpp.h"
#include <string.h>

#include<regex>
#include <alloca.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
using namespace std;

void ReAligner::set_reference(const string& reference) {
  this->reference_ = reference;
}

void ReAligner::set_reads(const vector<string>& reads) {
  this->reads_ = reads;
}

void ReAligner::set_ref_start(uint64_t position) {
  this->region_position_in_chr_ = position;
}

void ReAligner::set_haplotypes(const vector<string>& haplotypes) {
  this->haplotypes_ = haplotypes;
}

void ReAligner::set_options() {

    this->kmer_size_ = 32;
    this->read_size_ = 250;
    this->max_num_of_mismatches_ = 2;
    this->similarity_threshold_ = 0.16934;
    this->match_score_ = 4;
    this->mismatch_penalty_ = 6;
    this->gap_opening_penalty_ = 8;
    this->gap_extending_penalty_ = 2;
}
//
void ReAligner::CalculateSswAlignmentScoreThreshold() {
  ssw_alignment_score_threshold_ = match_score_
      * read_size_
      * similarity_threshold_
          - mismatch_penalty_
      * read_size_
      * (1 - similarity_threshold_);
  if (ssw_alignment_score_threshold_ < 0) {
    ssw_alignment_score_threshold_ = 1;
  }
}

vector<struct Read> ReAligner::AlignReads(
    const vector<struct Read>& reads_param) {

  for (const auto& read : reads_param) {
    reads_.push_back(read.seq);
  }
  CalculateSswAlignmentScoreThreshold();

  BuildIndex();

  FastAlignReadsToHaplotypes();

  InitSswLib();

  AlignHaplotypesToReference();

  CalculatePositionMaps();

  SswAlignReadsToHaplotypes(ssw_alignment_score_threshold_);

  sort(read_to_haplotype_alignments_.begin(),
            read_to_haplotype_alignments_.end());

  vector<struct Read> realigned_reads;

//  vector<struct Read>* realigned_reads(new vector<struct Read>());
  RealignReadsToReference(reads_param, realigned_reads);

  return realigned_reads;
}

void ReAligner::InitSswLib() {
    StripedSmithWaterman::Aligner ssw_aligner_(match_score_,
      mismatch_penalty_,
      gap_opening_penalty_,
      gap_extending_penalty_);
    StripedSmithWaterman::Filter filter;
}

void ReAligner::SswSetReference(const string& reference) {

  ssw_aligner_.SetReferenceSequence(reference.c_str(), reference.length());
}

StripedSmithWaterman::Alignment ReAligner::SswAlign(const string& target) const {

  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;

  int32_t maskLen = strlen(target.c_str())/2;
  maskLen = maskLen < 15 ? 15 : maskLen;
  if (ssw_aligner_.Align(target.c_str(), filter, &alignment)) {
    return alignment;
  } else {
    return StripedSmithWaterman::Alignment();
  }
}
//

void ReAligner::FastAlignReadsToHaplotypes() {
  vector<ReadAlignment> read_alignment_scores(reads_.size());
  for (int i = 0; i < haplotypes_.size(); i++) {
    const auto& haplotype = haplotypes_[i];
    int haplotype_score = 0;
    for (auto& readAlignment : read_alignment_scores) {
      readAlignment.reset();
    }
    FastAlignReadsToHaplotype(haplotype,
                              &haplotype_score,
                              &read_alignment_scores);

    if (haplotype_score == 0) {
      for (auto& readAlignment : read_alignment_scores) {
        readAlignment.reset();
      }
    }

    read_to_haplotype_alignments_.push_back(
        HaplotypeReadsAlignment(i, haplotype_score, read_alignment_scores));
  }
}

void ReAligner::FastAlignReadsToHaplotype(
    const string& haplotype, int* haplotype_score,
    vector<ReadAlignment>* haplotype_read_alignment_scores) {

  string bases_view = haplotype;

  bool is_ref = (haplotype == reference_);
  vector<int> coverage(haplotype.size(), 0);
  const auto& lastPos = haplotype.length() - kmer_size_;
  for (int i = 0; i <= lastPos; i++) {
    auto index_it = kmer_index_.find(bases_view.substr(i, kmer_size_));
    if (index_it == kmer_index_.end()) {
      continue;
    }
    for (const auto& it : index_it->second) {
      uint64_t read_id_index = static_cast<uint64_t>(it.read_id);
      int target_start_pos = max(
          static_cast<int64_t>(0),
          static_cast<int64_t>(i) - static_cast<int64_t>(it.read_pos.pos));
      int cur_read_size = reads_[read_id_index].size();
      int span = cur_read_size;
      if (target_start_pos + cur_read_size > haplotype.length()) {
        continue;
      }
      auto& read_alignment =
          (*haplotype_read_alignment_scores)[read_id_index];

      if (read_alignment.position != ReadAlignment::kNotAligned &&
          read_alignment.position == target_start_pos) {
        continue;
      }
      int num_of_mismatches = 0;
      int new_read_alignment_score = FastAlignStrings(
          bases_view.substr(target_start_pos, span),
          reads_[read_id_index],
          max_num_of_mismatches_ + 1, &num_of_mismatches);

      if (num_of_mismatches <= max_num_of_mismatches_) {
        int oldScore = read_alignment.score;
        for (auto pos = target_start_pos; pos < target_start_pos + span;
             pos++) {
          coverage[pos]++;
        }

        if (oldScore < new_read_alignment_score) {
          read_alignment.score = new_read_alignment_score;
          *haplotype_score -= oldScore;
          *haplotype_score += read_alignment.score;
          read_alignment.position = target_start_pos;
          read_alignment.cigar = to_string(cur_read_size) + "=";
        }
      }
    }

    if (coverage[i] == 0 && i >= ref_prefix_len_ &&
        i < haplotype.size() - ref_suffix_len_ && !is_ref) {
      *haplotype_score = 0;
      return;
    }
  }
}

int ReAligner::FastAlignStrings(string s1,
                                      string s2,
                                      int max_mismatches,
                                      int* num_of_mismatches) const {
  int num_of_matches = 0;
  *num_of_mismatches = 0;
  for (int i = 0; i < s1.size(); i++) {
    const auto& c1 = s1[i];
    const auto& c2 = s2[i];
    if (c1 != c2 && (c1 != 'N' && c2 != 'N')) {
      if (c1 != c2) {
        (*num_of_mismatches)++;
      }
      if (*num_of_mismatches == max_mismatches) {
        return 0;
      }
    } else {
      num_of_matches++;
    }
  }
  return num_of_matches * match_score_ - *num_of_mismatches * mismatch_penalty_;
}

Operation CigarOperationFromChar(char op) {
  switch (op) {
    case '=':
    case 'X':
      return ALIGNMENT_MATCH;
    case 'S':
      return CLIP_SOFT;
    case 'D':
      return DELETE;
    case 'I':
      return INSERT;
    default:
      return OPERATION_UNSPECIFIED;
  }
}


list<CigarOp> CigarStringToVector(const string& cigar) {
  list<CigarOp> cigarOps;
  string str = cigar;
  smatch result;
  string regex_str("(\\d+)([XIDS=])");
  regex pattern1(regex_str,regex::icase);

    string::const_iterator iter = str.begin();
    string::const_iterator iterEnd= str.end();
    string temp;
    while (regex_search(iter,iterEnd,result,pattern1)) {
        temp=result[0];
        int op_len = atoi(temp.substr(0, temp.length()-1).c_str());
        char op_char = temp[temp.length()-1];
        Operation op = CigarOperationFromChar(op_char);
        cigarOps.push_back(CigarOp(op, op_len));
        iter = result[0].second;
    }
  return cigarOps;
}


string CigarVectorToString(const list<CigarOp>& cigar) {
  string cigar_string = "";
  for (auto& op : cigar) {
      int len = op.length;
      int op_index = static_cast<int>(op.operation);
      string op_string = "";
      switch (op_index) {
        case 1:
          op_string = "X";
          break;
        case 2:
          op_string = "I";
          break;
        case 3:
          op_string = "D";
          break;
        case 5:
          op_string = "S";
          break;
        }
      cigar_string += to_string(len)+ op_string;
    }
  return cigar_string;
}



inline bool AlignmentIsRef(const string& cigar, int target_len) {
  return cigar == to_string(target_len) + "=";
}

void ReAligner::AlignHaplotypesToReference() {
  SswSetReference(reference_);

  if (read_to_haplotype_alignments_.empty()) {
    for (int i = 0; i < haplotypes_.size(); i++) {
      read_to_haplotype_alignments_.push_back(HaplotypeReadsAlignment(
          i, -1, vector<ReadAlignment>(reads_.size())));
    }
  }

  for (auto& haplotype_alignment : read_to_haplotype_alignments_) {
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment =
        SswAlign(haplotypes_[haplotype_alignment.haplotype_index]);
    auto hap_len = haplotypes_[haplotype_alignment.haplotype_index].size();
    if (alignment.sw_score > 0) {
      haplotype_alignment.is_reference =
          AlignmentIsRef(alignment.cigar_string, hap_len);
      haplotype_alignment.cigar = alignment.cigar_string;
      haplotype_alignment.cigar_ops =
          CigarStringToVector(haplotype_alignment.cigar);
      haplotype_alignment.ref_pos = alignment.ref_begin;
    }
  }
}

void ReAligner::SswAlignReadsToHaplotypes(int score_threshold) {

  for (int i = 0; i < reads_.size(); i++) {
    bool has_at_least_one_alignment = false;
    for (const auto& hap_alignment : read_to_haplotype_alignments_) {
      if (hap_alignment.read_alignment_scores[i].score > 0) {
        has_at_least_one_alignment = true;
        break;
      }
    }
    if (!has_at_least_one_alignment) {
      for (auto& hap_alignment : read_to_haplotype_alignments_) {

        if (hap_alignment.haplotype_score == 0) {
          continue;
        }
        SswSetReference(haplotypes_[hap_alignment.haplotype_index]);
        StripedSmithWaterman::Alignment alignment = SswAlign(reads_[i]);
        if (alignment.sw_score > 0) {
          if (alignment.sw_score >= score_threshold) {
            if (hap_alignment.read_alignment_scores[i].score <
                alignment.sw_score) {
              hap_alignment.read_alignment_scores[i].score = alignment.sw_score;
              hap_alignment.read_alignment_scores[i].cigar =
                  alignment.cigar_string;
              hap_alignment.read_alignment_scores[i].position =
                  alignment.ref_begin;
            }
          }
        }
      }
    }
  }  // for all reads
}

void ReAligner::RealignReadsToReference(
    const vector<struct Read>& reads,
    vector<struct Read>& realigned_reads) {

  for (int read_index = 0; read_index < reads.size(); read_index++) {
    const struct Read& read = reads[read_index];
    struct Read realigned_read = read;

    int best_hap_index = -1;

    if (GetBestReadAlignment(read_index, &best_hap_index)) {
      const HaplotypeReadsAlignment& bestHaplotypeAlignments =
          read_to_haplotype_alignments_[best_hap_index];

      int new_position;
      auto read_to_hap_pos = bestHaplotypeAlignments
          .read_alignment_scores[read_index]
          .position;

      new_position = region_position_in_chr_
              + bestHaplotypeAlignments.ref_pos
              + read_to_hap_pos
              + bestHaplotypeAlignments
                  .hap_to_ref_positions_map[read_to_hap_pos];
      list<CigarOp> readToRefCigarOps;

      CalculateReadToRefAlignment(
          read_index, bestHaplotypeAlignments.read_alignment_scores[read_index],
          bestHaplotypeAlignments.cigar_ops, &readToRefCigarOps);

      if (readToRefCigarOps.size() > 0) {
        realigned_read.cigar.clear();
        realigned_read.cigar_string = CigarVectorToString(readToRefCigarOps);
        realigned_read.position = new_position;
      }
      realigned_reads.push_back(realigned_read);

    } else {
      realigned_reads.push_back(realigned_read);
    }
  }  // for
}

void ReAligner::AddKmerToIndex(string kmer,
                                     ReadId read_id, KmerOffset pos) {
  kmer_index_[kmer].push_back(KmerOccurrence(read_id, pos));
}

void ReAligner::AddReadToIndex(const string& read, ReadId read_id) {

  if (read.length() <= kmer_size_) {
    return;
  }
  auto last_pos = read.length() - kmer_size_;
  string bases_view = read;
  for (int i = 0; i <= last_pos; i++) {
    AddKmerToIndex(bases_view.substr(i, kmer_size_), read_id, KmerOffset(i));
  }
}

void ReAligner::BuildIndex() {
  int read_id = 0;
  for (const auto& read : reads_) {
    AddReadToIndex(read, ReadId(read_id++));
  }
}

void SetPositionsMap(int haplotype_size,
                     HaplotypeReadsAlignment* hyplotype_alignment) {
  vector<int>& positions_map =
      hyplotype_alignment->hap_to_ref_positions_map;
  positions_map.resize(haplotype_size);
  string str = hyplotype_alignment->cigar;
  int cur_shift = 0;
  int haplotype_pos = 0;
  int last_pos = 0;
  int operation_len;
  string operation_type;

  smatch result;
  string regex_str("(\\d+)([XIDS=])");
  regex pattern1(regex_str,regex::icase);

    string::const_iterator iter = str.begin();
    string::const_iterator iterEnd = str.end();
    string temp;
    while (regex_search(iter,iterEnd,result,pattern1)) {
        temp=result[0];
        int operation_len = atoi(temp.substr(0, temp.length()-1).c_str());
        char op = temp[temp.length()-1];
        switch (op) {
          case '=':
          case 'X':
            last_pos = haplotype_pos + operation_len;
            while (haplotype_pos != last_pos) {
              positions_map[haplotype_pos] = cur_shift;
              haplotype_pos++;
            }
            break;
          case 'S':
            last_pos = haplotype_pos + operation_len;
            cur_shift -= operation_len;
            while (haplotype_pos != last_pos) {
              positions_map[haplotype_pos] = cur_shift;
              haplotype_pos++;
            }
            break;
          case 'D':
            cur_shift += operation_len;
            break;
          case 'I':
            last_pos = haplotype_pos + operation_len;
            while (haplotype_pos != last_pos) {
              positions_map[haplotype_pos] = cur_shift;
              cur_shift--;
              haplotype_pos++;
            }
            break;
    }
        iter = result[0].second;
    }
}
//
void ReAligner::CalculatePositionMaps() {
  for (auto& hyplotype_alignment : read_to_haplotype_alignments_) {
    SetPositionsMap(haplotypes_[hyplotype_alignment.haplotype_index].size(),
                    &hyplotype_alignment);
  }
}
//
bool ReAligner::GetBestReadAlignment(
    int readId,
    int* best_hap_index) const {
  int best_score = 0;
  bool best_haplotype_found = false;
  for (int hap_index = 0; hap_index < haplotypes_.size(); hap_index++) {
    if (read_to_haplotype_alignments_[hap_index]
                .read_alignment_scores[readId]
                .score > best_score
        || (best_score > 0 &&
            read_to_haplotype_alignments_[hap_index]
                    .read_alignment_scores[readId]
                    .score == best_score &&
            !read_to_haplotype_alignments_[hap_index].is_reference)) {
      best_score = read_to_haplotype_alignments_[hap_index]
                       .read_alignment_scores[readId]
                       .score;
      *best_hap_index = hap_index;
      best_haplotype_found = true;
    }
  }
  return best_haplotype_found;
}

int AlignedLength(const list<CigarOp>& cigar) {
  int len = 0;
  for (auto& op : cigar) {
    if (op.operation != DELETE) {
      len += op.length;
    }
  }
  return len;
}


void MergeCigarOp(const CigarOp& op, int read_len, list<CigarOp>* cigar) {
  const auto& last_cigar_op =
      cigar->empty() ? OPERATION_UNSPECIFIED
                     : cigar->back().operation;
  int aligned_length_before_merge = AlignedLength(*cigar);
  int new_op_length = 0;
  if (op.operation != DELETE) {
    new_op_length = min(op.length, read_len - aligned_length_before_merge);
  } else {
    new_op_length = op.length;
  }


  if (new_op_length <= 0 || aligned_length_before_merge == read_len) {
    return;
  }

  if (op.operation == last_cigar_op) {
    cigar->back().length += new_op_length;

  }  else {
    cigar->push_back(CigarOp(op.operation, new_op_length));
  }
}



list<CigarOp> LeftTrimHaplotypeToRefAlignment(
    const list<CigarOp>& haplotype_to_ref_cigar_ops_input,
    int read_to_haplotype_pos) {
  int cur_pos = 0;
  list<CigarOp> haplotype_to_ref_cigar_ops(
      haplotype_to_ref_cigar_ops_input);
  while (cur_pos != read_to_haplotype_pos) {
    CigarOp cur_hap_op = haplotype_to_ref_cigar_ops.front();
    haplotype_to_ref_cigar_ops.pop_front();
    if (cur_hap_op.operation ==
            ALIGNMENT_MATCH ||
        cur_hap_op.operation == CLIP_HARD ||
        cur_hap_op.operation == CLIP_SOFT ||
        cur_hap_op.operation == INSERT) {
      if (cur_hap_op.length + cur_pos > read_to_haplotype_pos) {
        haplotype_to_ref_cigar_ops.push_front(
            CigarOp(cur_hap_op.operation,
                    cur_hap_op.length - (read_to_haplotype_pos - cur_pos)));
      }
      cur_pos = min(cur_hap_op.length + cur_pos, read_to_haplotype_pos);
    }
  }

  if (haplotype_to_ref_cigar_ops.front().operation ==
      DELETE) {
    haplotype_to_ref_cigar_ops.pop_front();
  }

  return haplotype_to_ref_cigar_ops;
}

inline bool BothOpsAreMatch(const CigarOp& op1, const CigarOp& op2) {
  return (op1.operation == ALIGNMENT_MATCH ||
          op1.operation == CLIP_SOFT) &&
         (op2.operation == ALIGNMENT_MATCH ||
          op2.operation == CLIP_SOFT);
}

inline bool OneOfOpsIsSoftClip(const CigarOp& op1, const CigarOp& op2) {
  return op1.operation == CLIP_SOFT ||
         op2.operation == CLIP_SOFT;
}

inline bool DelAndMatch(const CigarOp& op1, const CigarOp& op2) {
  return op1.operation == DELETE &&
         (op2.operation == ALIGNMENT_MATCH ||
          op2.operation == CLIP_SOFT);
}

inline bool BothOpsAreDel(const CigarOp& op1, const CigarOp& op2) {
  return op1.operation == DELETE &&
         op2.operation == DELETE;
}

inline bool InsAndMatch(const CigarOp& op1, const CigarOp& op2) {
  return op1.operation == INSERT &&
         (op2.operation == ALIGNMENT_MATCH ||
          op2.operation == CLIP_SOFT);
}

inline bool BothOpsAreIns(const CigarOp& op1, const CigarOp& op2) {
  return (op1.operation == INSERT &&
          op2.operation == INSERT);
}

inline void PushFrontIfNotEmpty(const CigarOp& op, list<CigarOp>* cigar) {
  if (cigar == nullptr) {
    return;
  }
  if (op.length > 0) {
    cigar->push_front(op);
  }
}


void ReAligner::CalculateReadToRefAlignment(
    int read_index,
    const ReadAlignment& read_to_haplotype_alignment,
    const list<CigarOp>& haplotype_to_ref_cigar_ops_input,
    list<CigarOp>* read_to_ref_cigar_ops) const {
  int read_len = reads_[read_index].length();
  int read_to_haplotype_pos = read_to_haplotype_alignment.position;
  list<CigarOp> read_to_haplotype_cigar_ops =
      CigarStringToVector(read_to_haplotype_alignment.cigar);


  list<CigarOp> haplotype_to_ref_cigar_ops =
      LeftTrimHaplotypeToRefAlignment(haplotype_to_ref_cigar_ops_input,
                                      read_to_haplotype_pos);


  if (!read_to_haplotype_cigar_ops.empty() &&
      read_to_haplotype_cigar_ops.front().operation ==
          CLIP_SOFT) {
    MergeCigarOp(CigarOp(CLIP_SOFT,
                         read_to_haplotype_cigar_ops.front().length),
                 read_len, read_to_ref_cigar_ops);
    read_to_haplotype_cigar_ops.pop_front();
  }

  while ((!read_to_haplotype_cigar_ops.empty() ||
          !haplotype_to_ref_cigar_ops.empty()) &&
         AlignedLength(*read_to_ref_cigar_ops) < read_len) {

    if (!read_to_haplotype_cigar_ops.empty() &&
        haplotype_to_ref_cigar_ops.empty()) {
      MergeCigarOp(read_to_haplotype_cigar_ops.front(), read_len,
                   read_to_ref_cigar_ops);
      read_to_haplotype_cigar_ops.pop_front();
      continue;
    }

    if (read_to_haplotype_cigar_ops.empty() &&
        !haplotype_to_ref_cigar_ops.empty()) {
      break;
    }

    CigarOp cur_read_to_hap_op = read_to_haplotype_cigar_ops.front();
    read_to_haplotype_cigar_ops.pop_front();
    CigarOp cur_hap_to_ref_op = haplotype_to_ref_cigar_ops.front();
    haplotype_to_ref_cigar_ops.pop_front();


    // cur_read_to_hap_op, cur_hap_to_ref_op = M|S, M|S
    if (BothOpsAreMatch(cur_read_to_hap_op, cur_hap_to_ref_op)) {
      int new_op_len =
          min(cur_read_to_hap_op.length, cur_hap_to_ref_op.length);
      if (OneOfOpsIsSoftClip(cur_read_to_hap_op, cur_hap_to_ref_op)) {
        MergeCigarOp(
            CigarOp(CLIP_SOFT, new_op_len),
            read_len, read_to_ref_cigar_ops);
      } else {
        MergeCigarOp(CigarOp(ALIGNMENT_MATCH,
                             new_op_len),
                     read_len, read_to_ref_cigar_ops);
      }
      cur_read_to_hap_op.length -= new_op_len;
      PushFrontIfNotEmpty(cur_read_to_hap_op, &read_to_haplotype_cigar_ops);
      cur_hap_to_ref_op.length -= new_op_len;
      PushFrontIfNotEmpty(cur_hap_to_ref_op, &haplotype_to_ref_cigar_ops);

    // cur_read_to_hap_op, cur_hap_to_ref_op = D, M
    } else if (DelAndMatch(cur_read_to_hap_op, cur_hap_to_ref_op)) {
      MergeCigarOp(CigarOp(DELETE,
                           cur_read_to_hap_op.length),
                   read_len, read_to_ref_cigar_ops);
      cur_hap_to_ref_op.length -= cur_read_to_hap_op.length;
      PushFrontIfNotEmpty(cur_hap_to_ref_op, &haplotype_to_ref_cigar_ops);

    // cur_read_to_hap_op, cur_hap_to_ref_op = M, D
    } else if (DelAndMatch(cur_hap_to_ref_op, cur_read_to_hap_op)) {
      MergeCigarOp(CigarOp(DELETE,
                           cur_hap_to_ref_op.length),
                   read_len, read_to_ref_cigar_ops);
      PushFrontIfNotEmpty(cur_read_to_hap_op, &read_to_haplotype_cigar_ops);

    // cur_read_to_hap_op, cur_hap_to_ref_op = D, D
    } else if (BothOpsAreDel(cur_read_to_hap_op, cur_hap_to_ref_op)) {
      MergeCigarOp(
          CigarOp(DELETE,
                  cur_hap_to_ref_op.length + cur_read_to_hap_op.length),
          read_len, read_to_ref_cigar_ops);

      // cur_read_to_hap_op, cur_hap_to_ref_op = I, M
    } else if (InsAndMatch(cur_read_to_hap_op, cur_hap_to_ref_op)) {
      cur_read_to_hap_op.length =
          min(read_len - AlignedLength(*read_to_ref_cigar_ops),
                   cur_read_to_hap_op.length);
      MergeCigarOp(CigarOp(INSERT,
                           cur_read_to_hap_op.length),
                   read_len, read_to_ref_cigar_ops);
      PushFrontIfNotEmpty(cur_hap_to_ref_op, &haplotype_to_ref_cigar_ops);

    // cur_read_to_hap_op, cur_hap_to_ref_op = M, I
    } else if (InsAndMatch(cur_hap_to_ref_op, cur_read_to_hap_op)) {
      cur_hap_to_ref_op.length =
          min(read_len - AlignedLength(*read_to_ref_cigar_ops),
                   cur_hap_to_ref_op.length);
      MergeCigarOp(CigarOp(INSERT,
                           cur_hap_to_ref_op.length),
                   read_len, read_to_ref_cigar_ops);
      // We need to decrease the length of cur_read_to_hap_op by INS length
      cur_read_to_hap_op.length =
          max(0, cur_read_to_hap_op.length - cur_hap_to_ref_op.length);
      PushFrontIfNotEmpty(cur_read_to_hap_op, &read_to_haplotype_cigar_ops);

    // cur_read_to_hap_op, cur_hap_to_ref_op = I, I
    } else if (BothOpsAreIns(cur_hap_to_ref_op, cur_read_to_hap_op)) {
      cur_hap_to_ref_op.length =
          cur_hap_to_ref_op.length + cur_read_to_hap_op.length;
      MergeCigarOp(CigarOp(INSERT,
                           cur_hap_to_ref_op.length),
                   read_len, read_to_ref_cigar_ops);
    } else {

      read_to_ref_cigar_ops->clear();
      return;
    }
  }
}




struct_str_arr* ReAligner::realign_reads(char* seqs[], int* positions, char* cigars[], char* reference, char* haplotypes,  int ref_start, int ref_prefix, int ref_suffix, int read_size) {

    string reference_string = reference;
    string str = haplotypes;

    istringstream in_hal(str);
    vector<string> haplotypes_string;
    string t;
    while (in_hal >> t) {
        haplotypes_string.push_back(t);
    }

    vector<string> seq_list;
    for (int i=0; i < read_size; i++) {
        t = seqs[i];
        seq_list.push_back(t);
    }

    vector<int> position_list;
    for (int i=0; i < read_size; i++) {
        int pos = positions[i];
        position_list.push_back(pos);
    }


    vector<string> cigar_list;
    for (int i=0; i < read_size; i++) {
        t = cigars[i];
        cigar_list.push_back(t);
    }

    ReAligner aligner;
    aligner.set_options();
    aligner.set_reference(reference_string);
    aligner.set_haplotypes(haplotypes_string);
    aligner.set_ref_start(static_cast<uint64_t>(ref_start));
    aligner.set_ref_prefix_len(ref_prefix);
    aligner.set_ref_suffix_len(ref_suffix);
    vector<struct Read> Reads;
    for (int i = 0; i < read_size; i++) {
        string seq = seq_list[i];
        int position = position_list[i];
        string cigar_string = cigar_list[i];
        struct Read tmp_read(position, cigar_string, seq);

        list<CigarOp> cigar_vector = CigarStringToVector(cigar_string);
        for (auto& op : cigar_vector) {
            tmp_read.cigar.push_back(op);
        }
        Reads.push_back(tmp_read);
    }
    vector<struct Read> result = aligner.AlignReads(Reads);

    vector<string> result_cigar;
    vector<int> result_position;
    for (auto& read : result) {
        result_cigar.push_back(read.cigar_string);
        result_position.push_back(read.position);
    }
    struct_str_arr* str_arr_ptr = new struct_str_arr();

    for (int i=0; i < read_size; i++) {
        str_arr_ptr->cigar_string[i] = new char[result_cigar[i].size() + 1];
        strcpy(str_arr_ptr->cigar_string[i], result_cigar[i].c_str());
        str_arr_ptr->position[i] = result_position[i];
    }

    return str_arr_ptr;
    }



extern "C" {
    ReAligner obj;
    struct_str_arr* realign_reads(char * seqs[], int* positions, char* cigars[], char* reference, char* haplotypes,  int ref_start, int ref_prefix, int ref_suffix, int read_size) {
        return obj.realign_reads(seqs, positions, cigars, reference, haplotypes, ref_start, ref_prefix, ref_suffix, read_size);
    }
}

extern "C" {
    void free_memory(struct_str_arr* pointer, int size) {
        for (int i=0; i< size; i++) {
            delete [] pointer->cigar_string[i];
        }
        delete pointer;
        pointer = NULL;
    }
}


