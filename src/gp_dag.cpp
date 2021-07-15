// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag.hpp"
#include "numerical_utils.hpp"

using namespace GPOperations;  // NOLINT

GPDAG::PLVType RPLVType(bool rotated) {
  return rotated ? GPDAG::PLVType::R_TILDE : GPDAG::PLVType::R;
}

// The PLVTypes are documented where the enum is defined in the header file.
size_t GPDAG::GetPLVIndexStatic(PLVType plv_type, size_t node_count, size_t src_idx) {
  switch (plv_type) {
    case PLVType::P:
      return src_idx;
    case PLVType::P_HAT:
      return node_count + src_idx;
    case PLVType::P_HAT_TILDE:
      return 2 * node_count + src_idx;
    case PLVType::R_HAT:
      return 3 * node_count + src_idx;
    case PLVType::R:
      return 4 * node_count + src_idx;
    case PLVType::R_TILDE:
      return 5 * node_count + src_idx;
    default:
      Failwith("Invalid PLV index requested.");
  }
}

size_t GPDAG::GetPLVIndex(PLVType plv_type, size_t src_idx) const {
  return GetPLVIndexStatic(plv_type, dag_nodes_.size(), src_idx);
}

// The R PLV update that corresponds to our rotation status.
GPOperation GPDAG::RUpdateOfRotated(size_t node_id, bool rotated) const {
  return rotated ? Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                            GetPLVIndex(PLVType::R_HAT, node_id),
                            GetPLVIndex(PLVType::P_HAT, node_id)}
                 : Multiply{GetPLVIndex(PLVType::R, node_id),
                            GetPLVIndex(PLVType::R_HAT, node_id),
                            GetPLVIndex(PLVType::P_HAT_TILDE, node_id)};
}

// After this traversal, we will have optimized branch lengths, but we cannot assume
// that all of the PLVs are in a valid state.
//
// Update the terminology in this function as part of #288.
GPOperationVector GPDAG::ApproximateBranchLengthOptimization() const {
  GPOperationVector operations;
  SubsplitDAG::DepthFirstWithAction(SubsplitDAGTraversalAction(
      // BeforeNode
      [this, &operations](size_t node_id) {
        if (!GetDAGNode(node_id)->IsRootsplit()) {
          // Update R-hat if we're not at the root.
          UpdateRHat(node_id, operations);
        }
      },
      // AfterNode
      [this, &operations](size_t node_id) {
        // Make P the elementwise product ("o") of the two P PLVs for the node-clades.
        operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                      GetPLVIndex(PLVType::P_HAT, node_id),
                                      GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
      },
      // BeforeNodeClade
      [this, &operations](size_t node_id, bool rotated) {
        const PLVType p_hat_plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
        // Update the R PLV corresponding to our rotation status.
        operations.push_back(RUpdateOfRotated(node_id, rotated));
        // Zero out the node-clade PLV so we can fill it as part of VisitEdge.
        operations.push_back(ZeroPLV{GetPLVIndex(p_hat_plv_type, node_id)});
      },
      // VisitEdge
      [this, &operations](size_t node_id, size_t child_id, bool rotated) {
        // #310 this is temporary:
        // We do a full PLV population and then marginal likelihood calculation.
        // GPOperations::AppendGPOperations(operations, PopulatePLVs());
        // GPOperations::AppendGPOperations(operations, MarginalLikelihood());

        // Optimize each branch for a given node-clade and accumulate the resulting
        // P-hat PLVs in the parent node.
        OptimizeBranchLengthUpdatePHat(node_id, child_id, rotated, operations);
      }));
  return operations;
}

// After this traversal, we will have optimized branch lengths, but we cannot assume
// that all of the PLVs are in a valid state.
//
// Update the terminology in this function as part of #288.
GPOperationVector GPDAG::BranchLengthOptimization() {
  GPOperationVector operations;
  DepthFirstWithTidyAction(TidySubsplitDAGTraversalAction(
      // BeforeNode
      [this, &operations](size_t node_id) {
        if (!GetDAGNode(node_id)->IsRootsplit()) {
          // Update R-hat if we're not at the root.
          UpdateRHat(node_id, operations);
        }
      },
      // AfterNode
      [this, &operations](size_t node_id) {
        // Make P the elementwise product ("o") of the two P PLVs for the node-clades.
        operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                      GetPLVIndex(PLVType::P_HAT, node_id),
                                      GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
      },
      // BeforeNodeClade
      [this, &operations](size_t node_id, bool rotated) {
        const PLVType p_hat_plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
        // Update the R PLV corresponding to our rotation status.
        operations.push_back(RUpdateOfRotated(node_id, rotated));
        // Zero out the node-clade PLV so we can fill it as part of VisitEdge.
        operations.push_back(ZeroPLV{GetPLVIndex(p_hat_plv_type, node_id)});
      },
      // ModifyEdge
      [this, &operations](size_t node_id, size_t child_id, bool rotated) {
        // Optimize each branch for a given node-clade and accumulate the resulting
        // P-hat PLVs in the parent node.
        OptimizeBranchLengthUpdatePHat(node_id, child_id, rotated, operations);
      },
      // UpdateEdge
      [this, &operations](size_t node_id, size_t child_id, bool rotated) {
        // Accumulate all P-hat PLVs in the parent node without optimization.
        // #321 I don't think we need this Likelihood call... just the update PHat.
        UpdatePHatComputeLikelihood(node_id, child_id, rotated, operations);
      }));
  return operations;
}

GPOperationVector GPDAG::ComputeLikelihoods() const {
  GPOperationVector operations;
  IterateOverRealNodes([this, &operations](const SubsplitDAGNode *node) {
    IterateOverLeafwardEdges(
        node, [this, node, &operations](const bool rotated,
                                        const SubsplitDAGNode *child_node) {
          const auto gpcsp_idx = GPCSPIndexOfIds(node->Id(), child_node->Id());
          operations.push_back(Likelihood{gpcsp_idx,
                                          GetPLVIndex(RPLVType(rotated), node->Id()),
                                          GetPLVIndex(PLVType::P, child_node->Id())});
        });
  });

  const auto marginal_likelihood_operations = MarginalLikelihood();
  operations.insert(operations.end(), marginal_likelihood_operations.begin(),
                    marginal_likelihood_operations.end());
  return operations;
}

GPOperationVector GPDAG::LeafwardPass() const {
  return LeafwardPass(LeafwardPassTraversal());
}

GPOperationVector GPDAG::MarginalLikelihood() const {
  GPOperationVector operations = {GPOperations::ResetMarginalLikelihood{}};
  for (const auto &rootsplit_id : root_node_->GetLeafwardSorted()) {
    operations.push_back(GPOperations::IncrementMarginalLikelihood{
        GetPLVIndex(PLVType::R_HAT, rootsplit_id), RootsplitIndexOfId(rootsplit_id),
        GetPLVIndex(PLVType::P, rootsplit_id)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass() const {
  return RootwardPass(RootwardPassTraversal());
}

GPOperationVector GPDAG::OptimizeSBNParameters() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (size_t &id : LeafwardPassTraversal()) {
    const auto node = GetDAGNode(id);
    OptimizeSBNParametersForASubsplit(node->GetBitset(), operations);
    OptimizeSBNParametersForASubsplit(node->GetBitset().RotateSubsplit(), operations);
  }
  operations.push_back(UpdateSBNProbabilities{0, RootsplitCount()});
  return operations;
}

GPOperationVector GPDAG::SetLeafwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = 0; i < node_count; i++) {
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::R_HAT, i)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::R, i)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::R_TILDE, i)});
  }
  return operations;
}

GPOperationVector GPDAG::SetRhatToStationary() const {
  GPOperationVector operations;
  IterateOverRootsplitIds([this, &operations](size_t rootsplit_id) {
    size_t root_gpcsp_idx = RootsplitIndexOfId(rootsplit_id);
    operations.push_back(SetToStationaryDistribution{
        GetPLVIndex(PLVType::R_HAT, rootsplit_id), root_gpcsp_idx});
  });
  return operations;
}

GPOperationVector GPDAG::SetRootwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = taxon_count_; i < node_count; i++) {
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::P, i)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::P_HAT, i)});
    operations.push_back(ZeroPLV{GetPLVIndex(PLVType::P_HAT_TILDE, i)});
  }
  return operations;
}

GPOperationVector GPDAG::LeafwardPass(const SizeVector &visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDAGNode(node_id);

    // Build rhat(s) via rhat(s) += \sum_t q(s|t) P'(s|t) r(t)
    AddRhatOperations(node, operations);
    // Multiply to get r(s) = rhat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
    // Multiply to get r(s_tilde) = rhat(s) \circ phat(s).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass(const SizeVector &visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDAGNode(node_id);
    if (node->IsLeaf()) {
      continue;
    }
    // Build phat(s).
    AddPhatOperations(node, false, operations);
    // Build phat(s_tilde).
    AddPhatOperations(node, true, operations);
    // Multiply to get p(s) = phat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{node_id, GetPLVIndex(PLVType::P_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::PopulatePLVs() const {
  GPOperationVector operations;
  GPOperations::AppendGPOperations(operations, SetRootwardZero());
  GPOperations::AppendGPOperations(operations, SetLeafwardZero());
  GPOperations::AppendGPOperations(operations, SetRhatToStationary());
  GPOperations::AppendGPOperations(operations, RootwardPass());
  GPOperations::AppendGPOperations(operations, LeafwardPass());
  return operations;
}

// Take in some new operations, determine an appropriate PrepForMarginalization for
// them, then append the PrepForMarginalization and the new operations to `operations`
// (in that order).
void AppendOperationsAfterPrepForMarginalization(
    GPOperationVector &operations, const GPOperationVector &new_operations) {
  if (!new_operations.empty()) {
    operations.push_back(PrepForMarginalizationOfOperations(new_operations));
    operations.insert(operations.end(), new_operations.begin(), new_operations.end());
  }
}

void GPDAG::AddPhatOperations(const SubsplitDAGNode *node, bool rotated,
                              GPOperationVector &operations) const {
  PLVType plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  const auto parent_idx = node->Id();
  const size_t dest_idx = GetPLVIndex(plv_type, node->Id());
  GPOperationVector new_operations;
  for (const size_t &child_idx : node->GetLeafward(rotated)) {
    const auto gpcsp_idx = GPCSPIndexOfIds(parent_idx, child_idx);
    new_operations.push_back(IncrementWithWeightedEvolvedPLV{
        dest_idx, gpcsp_idx, GetPLVIndex(PLVType::P, child_idx)});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::AddRhatOperations(const SubsplitDAGNode *node,
                              GPOperationVector &operations) const {
  GPOperationVector new_operations;
  IterateOverRootwardEdges(
      node, [this, node, &new_operations](const bool rotated,
                                          const SubsplitDAGNode *parent_node) {
        new_operations.push_back(IncrementWithWeightedEvolvedPLV{
            GetPLVIndex(PLVType::R_HAT, node->Id()),
            GPCSPIndexOfIds(parent_node->Id(), node->Id()),
            GetPLVIndex(RPLVType(rotated), parent_node->Id())});
      });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                              GPOperationVector &operations) const {
  if (parent_to_range_.count(subsplit) > 0) {
    const auto param_range = parent_to_range_.at(subsplit);
    if (param_range.second - param_range.first > 1) {
      operations.push_back(
          UpdateSBNProbabilities{param_range.first, param_range.second});
    }
  }
}

void GPDAG::UpdateRHat(size_t node_id, GPOperationVector &operations) const {
  operations.push_back(ZeroPLV{GetPLVIndex(PLVType::R_HAT, node_id)});
  GPOperationVector new_operations;
  const auto node = GetDAGNode(node_id);
  for (const bool rotated : {false, true}) {
    PLVType src_plv_type = rotated ? PLVType::R_TILDE : PLVType::R;
    for (size_t parent_id : node->GetRootward(rotated)) {
      new_operations.push_back(IncrementWithWeightedEvolvedPLV{
          GetPLVIndex(PLVType::R_HAT, node_id), GPCSPIndexOfIds(parent_id, node_id),
          GetPLVIndex(src_plv_type, parent_id)});
    }
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

// #311 there's some work to be done here.
// There's a lot of common code between this function and the next.
// Also, the prep for marginalization isn't actually working correctly: we need to
// gather more operations first.
void GPDAG::UpdatePHatComputeLikelihood(size_t node_id, size_t child_node_id,
                                        bool rotated,
                                        GPOperationVector &operations) const {
  const auto gpcsp_idx = GPCSPIndexOfIds(node_id, child_node_id);
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  new_operations.push_back(Likelihood{gpcsp_idx,
                                      GetPLVIndex(RPLVType(rotated), node_id),
                                      GetPLVIndex(PLVType::P, child_node_id)});
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeBranchLengthUpdatePHat(size_t node_id, size_t child_node_id,
                                           bool rotated,
                                           GPOperationVector &operations) const {
  size_t gpcsp_idx = GPCSPIndexOfIds(node_id, child_node_id);
  operations.push_back(OptimizeBranchLength{GetPLVIndex(PLVType::P, child_node_id),
                                            GetPLVIndex(RPLVType(rotated), node_id),
                                            gpcsp_idx});
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

QuartetHybridRequest GPDAG::QuartetHybridRequestOf(size_t parent_id, bool rotated,
                                                   size_t child_id) const {
  QuartetTipVector rootward_tips;
  IterateOverRootwardEdgesAndParents(
      GetDAGNode(parent_id),
      [this, &rootward_tips](const size_t gpcsp_idx, const bool rootward_rotated,
                             const size_t grandparent_id) {
        rootward_tips.emplace_back(
            grandparent_id, GetPLVIndex(RPLVType(rootward_rotated), grandparent_id),
            gpcsp_idx);
      });

  QuartetTipVector sister_tips;
  const auto &parent_node = GetDAGNode(parent_id);
  const bool is_edge_to_sister_rotated = !rotated;
  IterateOverLeafwardEdges(
      parent_node, is_edge_to_sister_rotated,
      [this, &parent_node, &sister_tips,
       &is_edge_to_sister_rotated](const SubsplitDAGNode *sister_node) {
        const auto sister_id = sister_node->Id();
        sister_tips.emplace_back(
            sister_id, GetPLVIndex(PLVType::P, sister_id),
            GetGPCSPIndex(parent_node->GetBitset(is_edge_to_sister_rotated),
                          sister_node->GetBitset()));
      });

  QuartetTipVector rotated_tips;
  QuartetTipVector sorted_tips;
  IterateOverLeafwardEdgesAndChildren(
      GetDAGNode(child_id), [this, &rotated_tips, &sorted_tips](
                                const size_t gpcsp_idx, const bool leafward_rotated,
                                const size_t grandchild_id) {
        if (leafward_rotated) {
          rotated_tips.emplace_back(grandchild_id,
                                    GetPLVIndex(PLVType::P, grandchild_id), gpcsp_idx);
        } else {
          sorted_tips.emplace_back(grandchild_id,
                                   GetPLVIndex(PLVType::P, grandchild_id), gpcsp_idx);
        }
      });
  return QuartetHybridRequest(GPCSPIndexOfIds(parent_id, child_id),
                              std::move(rootward_tips), std::move(sister_tips),
                              std::move(rotated_tips), std::move(sorted_tips));
}
