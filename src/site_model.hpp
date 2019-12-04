// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SITE_MODEL_HPP_
#define SRC_SITE_MODEL_HPP_

#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include "block_model.hpp"

class SiteModel : public BlockModel {
 public:
  SiteModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~SiteModel() = default;

  virtual size_t GetCategoryCount() = 0;
  virtual const Eigen::VectorXd& GetCategoryRates() = 0;
  virtual const Eigen::VectorXd& GetCategoryProportions() = 0;

  static std::unique_ptr<SiteModel> OfSpecification(
      const std::string& specification);
};

class ConstantSiteModel : public SiteModel {
 public:
  ConstantSiteModel() : SiteModel({}) {
    one_.resize(1);
    one_[0] = 1.0;
  }

  size_t GetCategoryCount() override { return 1; }

  const Eigen::VectorXd& GetCategoryRates() override { return one_; }

  const Eigen::VectorXd& GetCategoryProportions() override { return one_; }

  void SetParameters(const EigenVectorXdRef param_vector) override{};

 private:
  Eigen::VectorXd one_;
};

class WeibullSiteModel : public SiteModel {
 public:
  explicit WeibullSiteModel(size_t category_count, double shape = 1)
      // Issue #147: Why are these commented out?
      : SiteModel({// {rates_key_, category_count},
                   // {proportions_key_, category_count},
                   {shape_key_, 1}}),
        category_count_(category_count),
        shape_(shape) {
    category_rates_.resize(category_count);
    category_proportions_.resize(category_count);
    for (int i = 0; i < category_count; i++) {
      category_proportions_[i] = 1.0 / category_count;
    }
    UpdateRates();
  }

  size_t GetCategoryCount() override;
  const Eigen::VectorXd& GetCategoryRates() override;
  const Eigen::VectorXd& GetCategoryProportions() override;

  void SetParameters(const EigenVectorXdRef param_vector) override;

  inline const static std::string rates_key_ = "Weibull category rates";
  inline const static std::string proportions_key_ =
      "Weibull category proportions";
  inline const static std::string shape_key_ = "Weibull shape";

 private:
  void UpdateRates();

  size_t category_count_;
  double shape_;  // shape of the Weibull distribution
  Eigen::VectorXd category_rates_;
  Eigen::VectorXd category_proportions_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SiteModel") {
  // Test 1: First we test using the "built in" default values.
  auto weibull_model = std::make_unique<WeibullSiteModel>(4, 1.0);
  const Eigen::VectorXd& rates = weibull_model->GetCategoryRates();
  Eigen::VectorXd rates_r(4);
  rates_r << 0.1457844, 0.5131316, 1.0708310, 2.2702530;
  CheckVectorXdEquality(rates, rates_r, 0.0001);

  // Test 2: Now set param_vector using SetParameters.
  weibull_model = std::make_unique<WeibullSiteModel>(4);
  Eigen::VectorXd param_vector(1);
  param_vector << 0.1;
  weibull_model->SetParameters(param_vector);
  rates_r << 4.766392e-12, 1.391131e-06, 2.179165e-03, 3.997819e+00;
  const Eigen::VectorXd& rates2 = weibull_model->GetCategoryRates();
  CheckVectorXdEquality(rates2, rates_r, 0.0001);

  // Test 3: Check proportions.
  const Eigen::VectorXd& proportions = weibull_model->GetCategoryProportions();
  for (size_t i = 0; i < proportions.size(); i++) {
    CHECK_LT(fabs(proportions[i] - 0.25), 0.0001);
  }

  // Test 4: Check sum rates[i]*proportions[i]==1.
  CHECK_LT(fabs(rates.dot(proportions) - 1.), 0.0001);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SITE_MODEL_HPP_
