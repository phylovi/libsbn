// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_SUBSTITUTION_MODEL_HPP_
#define SRC_SUBSTITUTION_MODEL_HPP_
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>
#include "block_model.hpp"
#include "sugar.hpp"

class SubstitutionModel : public BlockModel {
 public:
  SubstitutionModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~SubstitutionModel() = default;

  size_t GetStateCount() const { return frequencies_.size(); }

  const EigenMatrixXd& GetQMatrix() { return Q_; }

  const Eigen::VectorXd& GetFrequencies() { return frequencies_; }

  // We follow BEAGLE in terminology. "Inverse Eigenvectors" means the inverse
  // of the matrix containing the eigenvectors.
  const EigenMatrixXd& GetEigenvectors() { return eigenvectors_; }
  const EigenMatrixXd& GetInverseEigenvectors() {
    return inverse_eigenvectors_;
  }
  const Eigen::VectorXd& GetEigenvalues() { return eigenvalues_; }

  virtual void SetParameters(const EigenVectorXdRef param_vector) = 0;

  // This is the factory method that will be the typical way of buiding
  // substitution models.
  static std::unique_ptr<SubstitutionModel> OfSpecification(
      const std::string& specification);

 protected:
  Eigen::VectorXd frequencies_;
  EigenMatrixXd eigenvectors_;
  EigenMatrixXd inverse_eigenvectors_;
  Eigen::VectorXd eigenvalues_;
  EigenMatrixXd Q_;
};

class DNAModel : public SubstitutionModel {
 public:
  DNAModel(const BlockSpecification::ParamCounts& param_counts)
      : SubstitutionModel(param_counts) {
    frequencies_.resize(4);
    eigenvectors_.resize(4, 4);
    inverse_eigenvectors_.resize(4, 4);
    eigenvalues_.resize(4);
    Q_.resize(4, 4);
  }
};

class JC69Model : public DNAModel {
 public:
  JC69Model() : DNAModel({}) {
    frequencies_ << 0.25, 0.25, 0.25, 0.25;
    eigenvectors_ << 1.0, 2.0, 0.0, 0.5, 1.0, -2.0, 0.5, 0.0, 1.0, 2.0, 0.0,
        -0.5, 1.0, -2.0, -0.5, 0.0;
    inverse_eigenvectors_ << 0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125,
        -0.125, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0;
    eigenvalues_ << 0.0, -1.3333333333333333, -1.3333333333333333,
        -1.3333333333333333;
    Q_ << 1.0 / 3.0, -1.0, -1.0, -1.0, -1.0, 1.0 / 3.0, -1.0, -1.0, -1.0, -1.0,
        1.0 / 3.0, -1.0, -1.0, -1.0, -1.0, 1.0 / 3.0;
  }

  // No parameters to set for JC!
  void SetParameters(const EigenVectorXdRef param_vector){};
};

class GTRModel : public DNAModel {
 public:
  GTRModel() : DNAModel({{rates_key_, 6}, {frequencies_key_, 4}}) {
    rates_.resize(6);
    rates_ << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    frequencies_ << 0.25, 0.25, 0.25, 0.25;
    Update();
  }

  void SetParameters(const EigenVectorXdRef param_vector) override;

  inline const static std::string rates_key_ = "GTR rates";
  inline const static std::string frequencies_key_ = "frequencies";

 protected:
  void UpdateQMatrix();
  // Update the Q matrix _and_ the eigendecomposition.
  void Update();

 private:
  Eigen::VectorXd rates_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SubstitutionModel") {
  auto CheckEigenvalueEquality = [](Eigen::VectorXd eval1,
                                    Eigen::VectorXd eval2) {
    std::sort(eval1.begin(), eval1.end());
    std::sort(eval2.begin(), eval2.end());
    for (size_t i = 0; i < eval1.size(); i++) {
      CHECK_LT(fabs(eval1[i] - eval2[i]), 0.0001);
    }
  };
  auto gtr_model = std::make_unique<GTRModel>();
  auto jc_model = std::make_unique<JC69Model>();
  // Test 1: First we test using the "built in" default values.
  CheckEigenvalueEquality(jc_model->GetEigenvalues(),
                          gtr_model->GetEigenvalues());
  Eigen::VectorXd param_vector(10);
  // Test 2: Now set param_vector using SetParameters.
  gtr_model = std::make_unique<GTRModel>();
  param_vector << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25;
  gtr_model->SetParameters(param_vector);
  CheckEigenvalueEquality(jc_model->GetEigenvalues(),
                          gtr_model->GetEigenvalues());
  // Test 3: Now try out ParameterSegmentMapOf.
  gtr_model = std::make_unique<GTRModel>();
  // First zero out our param_vector.
  param_vector.setZero();
  // We can use ParameterSegmentMapOf to get two "views" into our parameter
  // vector.
  auto parameter_map =
      gtr_model->GetBlockSpecification().ParameterSegmentMapOf(param_vector);
  auto frequencies = parameter_map.at(GTRModel::frequencies_key_);
  auto rates = parameter_map.at(GTRModel::rates_key_);
  // When we modify the contents of these views, that changes param_vector.
  frequencies.setConstant(0.25);
  rates.setOnes();
  // We can then set param_vector and go forward as before.
  gtr_model->SetParameters(param_vector);
  CheckEigenvalueEquality(jc_model->GetEigenvalues(),
                          gtr_model->GetEigenvalues());
  // TODO @M: would it be easy for you to add an eigenvector/value test with
  // different coefficients given that you do the same thing in physher?
  // frequencies << 0.2, 0.55, 0.1, 0.15;
  // rates << 1., 0.5, 1., 2., 1., 1.;
  // gtr_model->SetParameters(param_vector);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSTITUTION_MODEL_HPP_
