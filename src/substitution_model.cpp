// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "substitution_model.hpp"

void GTRModel::SetParameters(Parameterization parameterization) {
  frequencies_ = GetFromParameterization(parameterization, "frequencies", 4);
  rates_ = GetFromParameterization(parameterization, "rates", 6);
  UpdateEigenDecomposition();
}

void GTRModel::UpdateQMatrix() {
  int rate_index = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      double rate = rates_[rate_index];
      rate_index++;
      Q_(i, j) = rate * frequencies_[j];
      Q_(j, i) = rate * frequencies_[i];
    }
  }
  // Set the diagonal entries so the rows sum to one.
  double total_substitution_rate = 0;
  for (int i = 0; i < 4; i++) {
    double row_sum = 0;
    for (int j = 0; j < 4; j++) {
      if (i != j) {
        row_sum += Q_(i, j);
      }
    }
    Q_(i, i) = -row_sum;
    total_substitution_rate += row_sum * frequencies_[i];
  }
  // Rescale matrix for unit substitution rate.
  Q_ /= total_substitution_rate;
}

void GTRModel::UpdateEigenDecomposition() {
  Eigen::Map<const Eigen::Array4d> tmp(&frequencies_[0]);
  EigenMatrixXd sqrt_frequencies =
      EigenMatrixXd(tmp.sqrt().matrix().asDiagonal());
  EigenMatrixXd sqrt_frequencies_inv =
      EigenMatrixXd(sqrt_frequencies.inverse());

  UpdateQMatrix();
  EigenMatrixXd S = EigenMatrixXd(sqrt_frequencies * Q_ * sqrt_frequencies_inv);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(S);

  evec_ = sqrt_frequencies_inv * solver.eigenvectors();
  ivec_ = solver.eigenvectors().transpose() * sqrt_frequencies;
  eval_ = solver.eigenvalues();
}