#include "tage.h"

#include <bitset>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{
constexpr uint8_t BIMODAL_INIT_VALUE = 1;
constexpr uint8_t TAGE_INIT_CTR = 4;
constexpr uint8_t TAGE_INIT_USEFUL = 0;
constexpr uint8_t BIMODAL_THRESHOLD = 2;
constexpr uint8_t TAGE_THRESHOLD = 4;
constexpr uint8_t USE_ALT_INIT = 8;
constexpr uint8_t USE_ALT_MAX = 15;
constexpr uint8_t BIMODAL_CTR_RANGE = 3;
constexpr uint8_t TAGE_CTR_RANGE = 7;
constexpr int MAX_PATH_HISTORY_BITS = 16;
} // namespace

void tage::initialize_branch_predictor()
{
  global_history.reset();
  path_history.reset();
  use_alt = USE_ALT_INIT;
  num_branches = 0;

  bimodal_table.resize(BIMODAL_TABLE_SIZE);
  for (auto& entry : bimodal_table) {
    entry = champsim::msl::fwcounter<2>(BIMODAL_INIT_VALUE);
  }

  double power = 1.0;
  for (int i = 0; i < COMPONENTS; i++) {
    int history_len = static_cast<int>(MIN_HISTORY_LENGTH * power + 0.5);
    power *= HISTORY_ALPHA;

    tables[i] = TABLE(1 << INDEX_BITS[i], history_len, TAG_BITS[i]);
    for (auto& entry : tables[i].rows) {
      entry.ctr = champsim::msl::fwcounter<3>(TAGE_INIT_CTR);
      entry.u = champsim::msl::fwcounter<2>(TAGE_INIT_USEFUL);
      entry.tag = 0;
    }
  }
}

uint64_t tage::compress_history(int in_size, int out_size)
{
  uint64_t compressed_history = 0;
  uint64_t temporary_history = 0;

  for (int i = 0; i < in_size && i < global_history.size(); i++) {
    if (i % out_size == 0) {
      compressed_history ^= temporary_history;
      temporary_history = 0;
    }
    temporary_history = (temporary_history << 1) | global_history[i];
  }
  compressed_history ^= temporary_history;
  return compressed_history;
}

uint16_t tage::get_bimodal_index(uint64_t ip) { return ip % BIMODAL_TABLE_SIZE; }

uint8_t tage::get_prediction(uint64_t ip, int comp)
{
  if (comp == 0) {
    uint16_t index = get_bimodal_index(ip);
    return bimodal_table[index].value() >= BIMODAL_THRESHOLD;
  }

  uint16_t index = compute_index(champsim::address{ip}, comp);
  return tables[comp - 1].rows[index].ctr.value() >= TAGE_THRESHOLD;
}

int tage::get_match_below_n(uint64_t ip, int component)
{
  for (int i = component - 1; i >= 1; i--) {
    uint16_t index = compute_index(champsim::address{ip}, i);
    uint16_t tag = compute_tag(champsim::address{ip}, i);

    if (tables[i - 1].rows[index].tag == tag) {
      return i;
    }
  }
  return 0;
}

void tage::ctr_update(uint8_t& ctr, bool cond, int low, int high)
{
  if (cond && ctr < high)
    ctr++;
  else if (!cond && ctr > low)
    ctr--;
}

uint64_t tage::get_path_hash(int component)
{
  uint64_t path_hash = 0;
  int effective_size = std::min(static_cast<int>(tables[component - 1].history_length), MAX_PATH_HISTORY_BITS);

  for (int i = 0; i < path_history.size() && i < effective_size; i++) {
    path_hash = (path_hash << 1) | path_history[i];
  }

  int index_width = INDEX_BITS[component - 1];
  uint64_t index_mask = (1ULL << index_width) - 1;

  uint64_t path1 = path_hash & index_mask;
  uint64_t path2 = path_hash >> index_width;

  path2 = ((path2 << component) & index_mask) + (path2 >> abs(index_width - component));
  path_hash = path1 ^ path2;
  path_hash = ((path_hash << component) & index_mask) + (path_hash >> abs(index_width - component));

  return path_hash;
}

uint64_t tage::compute_index(champsim::address ip, int table_id)
{
  if (table_id == 0) {
    return ip.to<uint64_t>() % tables[0].num_entries;
  }

  uint64_t pc = ip.to<uint64_t>();

  int index_width = INDEX_BITS[table_id - 1];

  uint64_t global_history_hash = compress_history(tables[table_id - 1].history_length, index_width);

  uint64_t path_history_hash = get_path_hash(table_id);

  uint64_t index = global_history_hash ^ pc ^ (pc >> (abs(index_width - table_id) + 1)) ^ path_history_hash;

  return index & ((1 << index_width) - 1);
}

uint64_t tage::compute_tag(champsim::address ip, int table_id)
{
  if (table_id == 0)
    return 0;

  uint64_t pc = ip.to<uint64_t>();
  int tag_width = tables[table_id - 1].tag_length;

  uint64_t global_history_hash = compress_history(tables[table_id - 1].history_length, tag_width);
  global_history_hash ^= compress_history(tables[table_id - 1].history_length, tag_width - 1);

  return (global_history_hash ^ pc) & ((1ULL << tag_width) - 1);
}

bool tage::predict_branch(champsim::address ip)
{
  uint64_t pc = ip.to<uint64_t>();

  pred_comp = get_match_below_n(pc, COMPONENTS + 1);
  alt_comp = get_match_below_n(pc, pred_comp);

  pred = get_prediction(pc, pred_comp);
  alt_pred = get_prediction(pc, alt_comp);

  if (pred_comp == 0) {
    tage_pred = pred;
  } else {
    uint16_t index = compute_index(champsim::address{pc}, pred_comp);
    STRONG = tables[pred_comp - 1].rows[index].ctr.value() > TAGE_THRESHOLD;

    tage_pred = (STRONG || use_alt < USE_ALT_INIT) ? pred : alt_pred;
  }

  return tage_pred;
}

void tage::last_branch_result(champsim::address ip, champsim::address branch_target, bool taken, uint8_t branch_type)
{
  uint64_t pc = ip.to<uint64_t>();

  if (pred_comp > 0) {
    uint16_t pred_index = compute_index(champsim::address{pc}, pred_comp);
    auto& pred_entry = tables[pred_comp - 1].rows[pred_index];
    uint8_t useful = pred_entry.u.value();

    if (!STRONG && pred != alt_pred) {
      if (pred != taken && use_alt < USE_ALT_MAX)
        use_alt++;
      else if (pred == taken && use_alt > 0)
        use_alt--;
    }

    if (alt_comp > 0) {
      uint16_t alt_index = compute_index(champsim::address{pc}, alt_comp);
      auto& alt_entry = tables[alt_comp - 1].rows[alt_index];
      if (useful == 0) {
        uint8_t temp_ctr = alt_entry.ctr.value();
        ctr_update(temp_ctr, taken, 0, TAGE_CTR_RANGE);
        alt_entry.ctr = champsim::msl::fwcounter<3>(temp_ctr);
      }
    } else if (useful == 0) {
      uint16_t bimodal_index = get_bimodal_index(pc);
      uint8_t temp_ctr = bimodal_table[bimodal_index].value();
      ctr_update(temp_ctr, taken, 0, BIMODAL_CTR_RANGE);
      bimodal_table[bimodal_index] = champsim::msl::fwcounter<2>(temp_ctr);
    }

    if (pred != alt_pred) {
      pred_entry.u += (pred == taken) ? 1 : -1;
    }

    uint8_t temp_ctr = pred_entry.ctr.value();
    ctr_update(temp_ctr, taken, 0, TAGE_CTR_RANGE);
    pred_entry.ctr = champsim::msl::fwcounter<3>(temp_ctr);
  } else {
    uint16_t bimodal_index = get_bimodal_index(pc);
    uint8_t temp_ctr = bimodal_table[bimodal_index].value();
    ctr_update(temp_ctr, taken, 0, BIMODAL_CTR_RANGE);
    bimodal_table[bimodal_index] = champsim::msl::fwcounter<2>(temp_ctr);
  }

  if (tage_pred != taken) {
    allocate_new_entries(pc);
  }

  global_history <<= 1;
  global_history[0] = taken;

  path_history <<= 1;
  path_history[0] = pc & 1;

  num_branches++;
  if (num_branches % RESET_USEFUL_INTERVAL == 0) {
    reset_usefulness_counters();
  }
}

void tage::allocate_new_entries(uint64_t pc)
{
  long rand_val = random() & ((1 << (COMPONENTS - pred_comp - 1)) - 1);
  int start_component = pred_comp + 1;

  if (rand_val & 1) {
    start_component++;
    if (rand_val & 2)
      start_component++;
  }

  bool is_free = false;
  for (int i = pred_comp + 1; i <= COMPONENTS; i++) {
    uint16_t index = compute_index(champsim::address{pc}, i);
    if (tables[i - 1].rows[index].u.value() == 0) {
      is_free = true;
      break;
    }
  }

  if (!is_free && start_component <= COMPONENTS) {
    uint16_t index = compute_index(champsim::address{pc}, start_component);
    tables[start_component - 1].rows[index].u = champsim::msl::fwcounter<2>(0);
  }

  for (int i = start_component; i <= COMPONENTS; i++) {
    uint16_t index = compute_index(champsim::address{pc}, i);
    auto& entry = tables[i - 1].rows[index];
    if (entry.u.value() == 0) {
      entry.tag = compute_tag(champsim::address{pc}, i);
      entry.ctr = champsim::msl::fwcounter<3>(TAGE_INIT_CTR);
      break;
    }
  }
}

void tage::reset_usefulness_counters()
{
  num_branches = 0;
  for (int i = 0; i < COMPONENTS; i++) {
    for (auto& entry : tables[i].rows) {
      entry.u = champsim::msl::fwcounter<2>(entry.u.value() >> 1);
    }
  }
}
