// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_PHYLO_MODEL_HPP_
#define SRC_PHYLO_MODEL_HPP_

#include <memory>

#include "clock_model.hpp"
#include "site_model.hpp"
#include "substitution_model.hpp"

struct PhyloModelSpecification {
  int substitution_model_specification_;
  int site_model_specification_;
  int clock_model_specification_;
};

class PhyloModel {
 public:
  PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
             std::unique_ptr<SiteModel> site_model,
             std::unique_ptr<ClockModel> clock_model);

  SubstitutionModel* GetSubstitutionModel() const {
    return substitution_model_.get();
  }
  SiteModel* GetSiteModel() const { return site_model_.get(); }
  ClockModel* GetClockModel() const { return clock_model_.get(); }

  static std::unique_ptr<PhyloModel> OfSpecification(
      const PhyloModelSpecification& specification) {
    return std::make_unique<PhyloModel>(
        std::make_unique<JCModel>(), std::make_unique<ConstantSiteModel>(),
        std::make_unique<StrictClockModel>(1.0));
  }

 private:
  std::unique_ptr<SubstitutionModel> substitution_model_;
  std::unique_ptr<SiteModel> site_model_;
  std::unique_ptr<ClockModel> clock_model_;
};

#endif  // SRC_PHYLO_MODEL_HPP_
