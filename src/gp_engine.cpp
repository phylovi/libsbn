// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t pcss_count)
    : site_pattern_(std::move(site_pattern)), pcss_count_(pcss_count) {
  plvs_ = std::vector<NucleotidePLV>(
      pcss_count_, NucleotidePLV::Zero(site_pattern_.PatternCount(), 4));
  branch_lengths_.resize(pcss_count_);
  likelihoods_.resize(pcss_count_);
  q_.resize(pcss_count_);

  InitializePLVsWithSitePatterns();
}

void GPEngine::InitializePLVsWithSitePatterns() {
  size_t taxon_idx = 0;
  for (const auto &pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative pattern!");
      if (symbol < 4) {
        plvs_[taxon_idx](site_idx, symbol) = 1.;
      }
      site_idx++;
    }
    taxon_idx++;
  }
}

