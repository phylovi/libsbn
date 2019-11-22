// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "fat_beagle.hpp"
#include <numeric>
#include <utility>
#include <vector>

FatBeagle::FatBeagle(const PhyloModelSpecification &specification,
                     const SitePattern &site_pattern)
    : phylo_model_(PhyloModel::OfSpecification(specification)),
      rescaling_(false),
      // TODO @M do we have to store pattern_count_? It seems like the answer is
      // no here, but I seem to remember you might need it in the future.
      pattern_count_(static_cast<int>(site_pattern.PatternCount())) {
  beagle_instance_ = CreateInstance(site_pattern);
  SetTipStates(site_pattern);
  UpdatePhyloModelInBeagle();
};

FatBeagle::~FatBeagle() {
  auto finalize_result = beagleFinalizeInstance(beagle_instance_);
  if (finalize_result != 0) {
    std::cout << "beagleFinalizeInstance gave nonzero return value!";
    std::terminate();
  }
}

const BlockSpecification &FatBeagle::GetBlockSpecification() const {
  return phylo_model_->GetBlockSpecification();
}

void FatBeagle::SetParameters(const EigenVectorXdRef param_vector) {
  phylo_model_->SetParameters(param_vector);
  UpdatePhyloModelInBeagle();
}

FatBeagle::BeagleInstance FatBeagle::CreateInstance(
    const SitePattern &site_pattern) {
  int tip_count = static_cast<int>(site_pattern.SequenceCount());
  // Number of partial buffers to create (input):
  // tip_count - 1 for lower partials (internal nodes only)
  // 2*tip_count - 2 for upper partials (every node except the root)
  int partials_buffer_count = 3 * tip_count - 3;
  // Number of compact state representation buffers to create -- for use with
  // setTipStates (input) */
  int compact_buffer_count = tip_count;
  // DNA assumption here.
  int state_count =
      static_cast<int>(phylo_model_->GetSubstitutionModel()->GetStateCount());
  // Number of site patterns to be handled by the instance (input) -- not
  // compressed in this case
  int pattern_count = pattern_count_;
  // Number of eigen-decomposition buffers to allocate (input)
  int eigen_buffer_count = 1;
  // Number of transition matrix buffers (input) -- two per edge
  int matrix_buffer_count = 2 * (2 * tip_count - 1);
  // Number of rate categories
  int category_count =
      static_cast<int>(phylo_model_->GetSiteModel()->GetCategoryCount());
  // Number of scaling buffers -- 1 buffer per partial buffer and 1 more
  // for accumulating scale factors in position 0.
  int scale_buffer_count = partials_buffer_count + 1;
  // List of potential resources on which this instance is allowed (input,
  // NULL implies no restriction
  int *allowed_resources = nullptr;
  // Length of resourceList list (input) -- not needed to use the default
  // hardware config
  int resource_count = 0;
  // Bit-flags indicating preferred implementation charactertistics, see
  // BeagleFlags (input)
  int64_t preference_flags = 0;
  // Bit-flags indicating required implementation characteristics, see
  // BeagleFlags (input)
  int requirement_flags = BEAGLE_FLAG_SCALING_MANUAL;

  BeagleInstanceDetails return_info;
  auto beagle_instance = beagleCreateInstance(
      tip_count, partials_buffer_count, compact_buffer_count, state_count,
      pattern_count, eigen_buffer_count, matrix_buffer_count, category_count,
      scale_buffer_count, allowed_resources, resource_count, preference_flags,
      requirement_flags, &return_info);
  if (return_info.flags &
      (BEAGLE_FLAG_PROCESSOR_CPU | BEAGLE_FLAG_PROCESSOR_GPU)) {
    return beagle_instance;
  }  // else
  Failwith("Couldn't get a CPU or a GPU from BEAGLE.");
}

void FatBeagle::SetTipStates(const SitePattern &site_pattern) {
  int taxon_number = 0;
  for (const auto &pattern : site_pattern.GetPatterns()) {
    beagleSetTipStates(beagle_instance_, taxon_number++, pattern.data());
  }
  beagleSetPatternWeights(beagle_instance_, site_pattern.GetWeights().data());
}

void FatBeagle::UpdateSiteModelInBeagle() {
  const std::vector<double> &weights =
      phylo_model_->GetSiteModel()->GetCategoryProportions();
  const std::vector<double> &rates =
      phylo_model_->GetSiteModel()->GetCategoryRates();
  beagleSetCategoryWeights(beagle_instance_, 0, weights.data());
  beagleSetCategoryRates(beagle_instance_, rates.data());
}

void FatBeagle::UpdateSubstitutionModelInBeagle() {
  const auto &substitution_model = phylo_model_->GetSubstitutionModel();
  const EigenMatrixXd &eigen_vectors = substitution_model->GetEigenvectors();
  const EigenMatrixXd &inverse_eigen_vectors =
      substitution_model->GetInverseEigenvectors();
  const Eigen::VectorXd &eigen_values = substitution_model->GetEigenvalues();
  const Eigen::VectorXd &frequencies = substitution_model->GetFrequencies();

  beagleSetStateFrequencies(beagle_instance_, 0, frequencies.data());
  beagleSetEigenDecomposition(beagle_instance_,
                              0,  // eigenIndex
                              &eigen_vectors.data()[0],
                              &inverse_eigen_vectors.data()[0],
                              &eigen_values.data()[0]);
}

void FatBeagle::UpdatePhyloModelInBeagle() {
  // TODO Update clock model.
  UpdateSiteModelInBeagle();
  UpdateSubstitutionModelInBeagle();
}

Tree PrepareTreeForLikelihood(const Tree &tree) {
  if (tree.Children().size() == 3) {
    return tree.Detrifurcate();
  }  // else
  if (tree.Children().size() == 2) {
    return tree;
  }
  // else
  Failwith(
      "Tree likelihood calculations should be done on a tree with a "
      "bifurcation or a trifurcation at the root.");
}

double FatBeagle::LogLikelihood(const Tree &in_tree) const {
  beagleResetScaleFactors(beagle_instance_, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  BeagleAccessories ba(beagle_instance_, rescaling_, tree);
  std::vector<BeagleOperation> operations;
  tree.Topology()->PostOrder(
      [&operations, &ba](const Node *node) {
        if (!node->IsLeaf()) {
          Assert(node->Children().size() == 2,
                 "Tree isn't bifurcating in LogLikelihood.");
          int dest = static_cast<int>(node->Id());
          int child0_id = static_cast<int>(node->Children()[0]->Id());
          int child1_id = static_cast<int>(node->Children()[1]->Id());
          BeagleOperation op = {
              dest,                            // dest
              BEAGLE_OP_NONE, BEAGLE_OP_NONE,  // src and dest scaling
              child0_id,      child0_id,       // src1 and matrix1
              child1_id,      child1_id        // src2 and matrix2
          };
          if (ba.rescaling_) {
            // We don't need scaling buffers for the leaves.
            // Index 0 is reserved for accumulating the sum of log scalers.
            // Thus the scaling buffers are indexed by the edge number minus the
            // number of leaves + 1.
            op.destinationScaleWrite = dest - ba.taxon_count_ + 1;
          }
          operations.push_back(op);
        }
      });
  beagleUpdateTransitionMatrices(beagle_instance_,
                                 0,  // eigenIndex
                                 ba.node_indices_.data(),
                                 nullptr,  // firstDerivativeIndices
                                 nullptr,  // secondDervativeIndices
                                 tree.BranchLengths().data(),
                                 ba.node_count_ - 1);

  beagleUpdatePartials(beagle_instance_,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       ba.cumulative_scale_index_[0]);

  double log_like = 0;
  std::vector<int> root_id = {static_cast<int>(tree.Id())};
  beagleCalculateRootLogLikelihoods(
      beagle_instance_, root_id.data(), ba.category_weight_index_.data(),
      ba.state_frequency_index_.data(), ba.cumulative_scale_index_.data(),
      ba.mysterious_count_, &log_like);
  return log_like;
}

std::pair<double, std::vector<double>> FatBeagle::BranchGradient(
    const Tree &in_tree) const {
  beagleResetScaleFactors(beagle_instance_, 0);
  auto tree = PrepareTreeForLikelihood(in_tree);
  tree.SlideRootPosition();

  BeagleAccessories ba(beagle_instance_, rescaling_, tree);
  std::vector<BeagleOperation> operations;
  std::vector<int> gradient_indices(ba.internal_count_);
  std::iota(gradient_indices.begin(), gradient_indices.end(), ba.node_count_);

  // Calculate lower partials.
  tree.Topology()->BinaryIdPostOrder(
      [&operations, &ba](int node_id, int child0_id, int child1_id) {
        BeagleOperation op = {
            node_id,                         // dest
            BEAGLE_OP_NONE, BEAGLE_OP_NONE,  // src and dest scaling
            child0_id,      child0_id,       // src1 and matrix1
            child1_id,      child1_id        // src2 and matrix2
        };
        if (ba.rescaling_) {
          op.destinationScaleWrite = node_id - ba.taxon_count_ + 1;
        }
        operations.push_back(op);
      });

  // Calculate upper partials.
  tree.Topology()->TripleIdPreOrderBifurcating(
      [&operations, &ba](int node_id, int sister_id, int parent_id) {
        if (node_id != ba.root_child_id_ && node_id != ba.fixed_node_id_) {
          int upper_partial_index;
          int upper_matrix_index;
          if (parent_id == ba.root_child_id_) {
            upper_matrix_index = static_cast<int>(ba.root_child_id_);
            upper_partial_index = static_cast<int>(ba.fixed_node_id_);
          } else if (parent_id == ba.fixed_node_id_) {
            upper_matrix_index = static_cast<int>(ba.root_child_id_);
            upper_partial_index = static_cast<int>(ba.root_child_id_);
          } else {
            upper_partial_index = parent_id + ba.node_count_;
            upper_matrix_index = parent_id;
          }
          BeagleOperation op = {node_id + ba.node_count_,
                                BEAGLE_OP_NONE,
                                BEAGLE_OP_NONE,
                                upper_partial_index,
                                upper_matrix_index,
                                sister_id,
                                sister_id};
          // Scalers are indexed differently for the upper conditional
          // likelihood. They start at the number of internal nodes + 1 because
          // of the lower conditional likelihoods. Also, in this case the leaves
          // have scalers.
          if (ba.rescaling_) {
            // Scaling factors are recomputed every time so we don't read them
            // using destinationScaleRead.
            op.destinationScaleWrite = node_id + 1 + ba.internal_count_;
          }
          operations.push_back(op);
        }
      });

  beagleUpdateTransitionMatrices(
      beagle_instance_,
      0,  // eigenIndex
      ba.node_indices_.data(),
      gradient_indices.data(),  // firstDerivativeIndices
      nullptr,                  // secondDervativeIndices
      tree.BranchLengths().data(), ba.node_count_ - 1);

  beagleUpdatePartials(beagle_instance_,
                       operations.data(),  // eigenIndex
                       static_cast<int>(operations.size()),
                       BEAGLE_OP_NONE);  // cumulative scale index

  std::vector<double> gradient(ba.node_count_, 0);
  double log_like = 0;

  SizeVectorVector indices_above = tree.Topology()->IdsAbove();
  for (auto &indices : indices_above) {
    // Reverse vector so we have indices from node to root.
    std::reverse(indices.begin(), indices.end());
    // Remove the root scalers.
    indices.pop_back();
  }

  // Actually compute gradient.
  tree.Topology()->TripleIdPreOrderBifurcating([&ba, &gradient, &indices_above,
                                                &log_like](int node_id,
                                                           int sister_id, int) {
    if (node_id == ba.fixed_node_id_) {
      return;
    }
    double dlogLp;
    ba.upper_partials_index_[0] = node_id + ba.node_count_;
    ba.node_partial_indices_[0] = node_id;
    ba.node_mat_indices_[0] = node_id;
    ba.node_deriv_index_[0] = node_id + ba.node_count_;

    if (node_id == ba.root_child_id_) {
      ba.upper_partials_index_[0] = sister_id;
    }
    // parent partial Buffers cannot be a taxon in
    // beagleCalculateEdgeLogLikelihoods
    if (ba.node_partial_indices_[0] > ba.upper_partials_index_[0]) {
      std::swap(ba.node_partial_indices_, ba.upper_partials_index_);
    }

    if (ba.rescaling_) {
      beagleResetScaleFactors(ba.beagle_instance_,
                              ba.cumulative_scale_index_[0]);
      std::vector<int> scaler_indices(
          static_cast<size_t>(ba.internal_count_ - 1));
      std::iota(scaler_indices.begin(), scaler_indices.end(), 1);
      // Replace lower scaler index with upper scaler index for nodes between
      // node_index and root.
      int child = node_id;
      for (size_t upper : indices_above[static_cast<size_t>(node_id)]) {
        int int_upper = static_cast<int>(upper);
        int scaler_indices_index = int_upper - ba.taxon_count_;
        Assert(scaler_indices_index >= 0, "int_upper must be >= taxon count.");
        scaler_indices[static_cast<size_t>(scaler_indices_index)] =
            child + ba.internal_count_ + 1;
        child = int_upper;
      }
      beagleAccumulateScaleFactors(ba.beagle_instance_, scaler_indices.data(),
                                   static_cast<int>(scaler_indices.size()),
                                   ba.cumulative_scale_index_[0]);
    }

    beagleCalculateEdgeLogLikelihoods(
        ba.beagle_instance_, ba.upper_partials_index_.data(),
        ba.node_partial_indices_.data(), ba.node_mat_indices_.data(),
        ba.node_deriv_index_.data(),
        nullptr,  // second derivative matrices
        ba.category_weight_index_.data(), ba.state_frequency_index_.data(),
        ba.cumulative_scale_index_.data(), ba.mysterious_count_, &log_like,
        &dlogLp,
        nullptr);  // destination for second derivative
    gradient[static_cast<size_t>(node_id)] = dlogLp;
  });

  return {log_like, gradient};
}

double FatBeagle::StaticLogLikelihood(FatBeagle *fat_beagle,
                                      const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "NULL FatBeagle pointer!");
  return fat_beagle->LogLikelihood(in_tree);
}

std::pair<double, std::vector<double>> FatBeagle::StaticBranchGradient(
    FatBeagle *fat_beagle, const Tree &in_tree) {
  Assert(fat_beagle != nullptr, "NULL FatBeagle pointer!");
  return fat_beagle->BranchGradient(in_tree);
}
