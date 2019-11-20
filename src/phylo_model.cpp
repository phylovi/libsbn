// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "phylo_model.hpp"

PhyloModel::PhyloModel(std::unique_ptr<SubstitutionModel> substitution_model,
                       std::unique_ptr<SiteModel> site_model,
                       std::unique_ptr<ClockModel> clock_model)
    : BlockModel({}),
      substitution_model_(std::move(substitution_model)),
      site_model_(std::move(site_model)),
      clock_model_(std::move(clock_model)) {
  size_t next_available_idx = 0;
  Incorporate(next_available_idx, substitution_entire_key_,
              substitution_model_->GetBlockSpecification());
  Incorporate(next_available_idx, site_entire_key_,
              site_model_->GetBlockSpecification());
  Incorporate(next_available_idx, clock_entire_key_,
              clock_model_->GetBlockSpecification());
  InsertEntire({0, next_available_idx});
}

std::unique_ptr<PhyloModel> PhyloModel::OfSpecification(
    const PhyloModelSpecification& specification) {
  return std::make_unique<PhyloModel>(
      SubstitutionModel::OfSpecification(specification.substitution_),
      SiteModel::OfSpecification(specification.site_),
      ClockModel::OfSpecification(specification.clock_));
}

void PhyloModel::SetParameters(const EigenVectorXdRef parameterization) {
  substitution_model_->SetParameters(
      ExtractSegment(parameterization, substitution_entire_key_));
  site_model_->SetParameters(
      ExtractSegment(parameterization, site_entire_key_));
  clock_model_->SetParameters(
      ExtractSegment(parameterization, clock_entire_key_));
}

