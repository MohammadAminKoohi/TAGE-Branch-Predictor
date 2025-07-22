#ifndef BRANCH_TAGE_H
#define BRANCH_TAGE_H

#include <array>
#include <vector>

#include "modules.h"
#include "msl/fwcounter.h"

using namespace std;

typedef struct {
  unsigned long long tag;
  champsim::msl::fwcounter<3> ctr;
  champsim::msl::fwcounter<2> u;
} ENTRY;

struct TABLE {
  size_t num_entries;
  size_t history_length;
  size_t tag_length;

  vector<ENTRY> rows;

  TABLE() : num_entries(0), history_length(0), tag_length(0) {}

  TABLE(size_t num_entries, size_t history_length, size_t tag_length)
  {
    this->num_entries = num_entries;
    this->history_length = history_length;
    this->tag_length = tag_length;
    rows.resize(num_entries);
  }
};

class tage : champsim::modules::branch_predictor
{
  constexpr static size_t COMPONENTS = 8;
  constexpr static size_t BIMODAL_TABLE_SIZE = 16384;
  constexpr static double HISTORY_ALPHA = 1.6;
  constexpr static int MIN_HISTORY_LENGTH = 4;
  constexpr static int RESET_USEFUL_INTERVAL = 512000;

  const uint8_t INDEX_BITS[12] = {10, 10, 11, 11, 11, 11, 10, 10, 10, 10, 9, 9};
  const uint8_t TAG_BITS[12] = {7, 7, 8, 8, 9, 10, 11, 12, 12, 13, 14, 15};

  int RESET_COUNTER = 0;
  int num_branches = 0;
  bitset<1024> global_history;
  bitset<32> path_history;
  int provider_table = -1;
  int alternate_table = -1;
  uint8_t use_alt = 8;

  uint8_t tage_pred, pred, alt_pred;
  int pred_comp, alt_comp;
  int STRONG;

  std::vector<champsim::msl::fwcounter<2>> bimodal_table;
  array<TABLE, COMPONENTS> tables;

  uint64_t compute_index(champsim::address ip, int table_id);
  uint64_t compute_tag(champsim::address ip, int table_id);
  uint64_t compress_history(int in_size, int out_size);
  uint64_t get_path_hash(int component);
  uint16_t get_bimodal_index(uint64_t ip);
  int get_match_below_n(uint64_t ip, int component);
  void ctr_update(uint8_t& ctr, bool cond, int low, int high);
  uint8_t get_prediction(uint64_t ip, int comp);
  void allocate_new_entries(uint64_t pc);
  void reset_usefulness_counters();

public:
  using branch_predictor::branch_predictor;

  void initialize_branch_predictor();
  bool predict_branch(champsim::address ip);
  void last_branch_result(champsim::address ip, champsim::address branch_target, bool taken, uint8_t branch_type);
  void print_stats();
};

#endif