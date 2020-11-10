// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"

#include <chrono>
#include <iomanip>
#include <string>

#include "csv.hpp"
#include "driver.hpp"
#include "gp_operation.hpp"
#include "numerical_utils.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_probability.hpp"

using namespace GPOperations;  // NOLINT

void GPInstance::PrintStatus() {
  const auto tree_count = tree_collection_.TreeCount();
  const auto taxon_count = tree_collection_.TaxonCount();
  if (tree_count > 0) {
    std::cout << tree_count << " trees loaded on " << taxon_count << " leaves.\n";
  } else {
    std::cout << "No trees loaded.\n";
  }
  std::cout << alignment_.Data().size() << " sequences loaded.\n";
  std::cout << dag_.NodeCount() << " DAG nodes representing " << dag_.TopologyCount()
            << " trees.\n";
  std::cout << dag_.GPCSPCountWithFakeSubsplits() << " continuous parameters.\n";
  if (HasEngine()) {
    std::cout << "Engine available using " << GetEngine()->PLVByteCount() / 1e9
              << "G virtual memory.\n";
  } else {
    std::cout << "Engine has not been made.\n";
  }
}

void GPInstance::ReadFastaFile(const std::string &fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void GPInstance::ReadNewickFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void GPInstance::ReadNexusFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

void GPInstance::CheckSequencesAndTreesLoaded() const {
  if (alignment_.SequenceCount() == 0) {
    Failwith(
        "Load an alignment into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
  if (tree_collection_.TreeCount() == 0) {
    Failwith(
        "Load some trees into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
}

void GPInstance::MakeEngine(double rescaling_threshold) {
  CheckSequencesAndTreesLoaded();
  ProcessLoadedTrees();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());

  dag_ = GPDAG(tree_collection_);
  engine_ = std::make_unique<GPEngine>(
      std::move(site_pattern), plv_count_per_node_ * dag_.NodeCount(),
      dag_.GPCSPCountWithFakeSubsplits(), mmap_file_path_, rescaling_threshold);
  InitializeGPEngine();
}

GPEngine *GPInstance::GetEngine() const {
  if (engine_ != nullptr) {
    return engine_.get();
  }
  // else
  Failwith(
      "Engine not available. Call MakeEngine to make an engine for phylogenetic "
      "likelihood computation.");
}

bool GPInstance::HasEngine() const { return engine_ != nullptr; }

void GPInstance::ProcessOperations(const GPOperationVector &operations) {
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::ClearTreeCollectionAssociatedState() {
  sbn_parameters_.resize(0);
  dag_ = GPDAG();
}

void GPInstance::ProcessLoadedTrees() {
  sbn_parameters_.resize(dag_.GPCSPCount());
  sbn_parameters_.setOnes();
}

void GPInstance::HotStartBranchLengths() {
  if (HasEngine()) {
    GetEngine()->HotStartBranchLengths(tree_collection_, dag_.GetGPCSPIndexer());
  } else {
    Failwith(
        "Please load and process some trees before calling HotStartBranchLengths.");
  }
}

void GPInstance::PrintDAG() { dag_.Print(); }

void GPInstance::PrintGPCSPIndexer() {
  std::cout << "Vector of taxon names: " << tree_collection_.TaxonNames() << std::endl;
  dag_.PrintGPCSPIndexer();
}

void GPInstance::InitializeGPEngine() {
  GetEngine()->SetSBNParameters(dag_.BuildUniformPrior());
}

void GPInstance::ResetMarginalLikelihoodAndPopulatePLVs() {
  GetEngine()->ResetLogMarginalLikelihood();
  ProcessOperations(dag_.SetRootwardZero());
  ProcessOperations(dag_.SetLeafwardZero());
  ProcessOperations(dag_.SetRhatToStationary());
  ProcessOperations(dag_.RootwardPass());
  ProcessOperations(dag_.LeafwardPass());
}

void GPInstance::ComputeLikelihoods() { ProcessOperations(dag_.ComputeLikelihoods()); }

void GPInstance::EstimateBranchLengths(double tol, size_t max_iter) {
  auto now = std::chrono::high_resolution_clock::now;
  auto t_start = now();
  std::cout << "Begin branch optimization\n";
  GPOperationVector branch_optimization_operations = dag_.BranchLengthOptimization();
  GPOperationVector marginal_lik_operations = dag_.MarginalLikelihood();

  std::cout << "Populating PLVs\n";
  ResetMarginalLikelihoodAndPopulatePLVs();
  std::chrono::duration<double> warmup_duration = now() - t_start;
  t_start = now();
  std::cout << "Computing initial likelihood\n";
  ProcessOperations(marginal_lik_operations);
  double current_marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
  std::chrono::duration<double> initial_likelihood_duration = now() - t_start;
  t_start = now();

  for (size_t i = 0; i < max_iter; i++) {
    std::cout << "Iteration: " << (i + 1) << std::endl;
    ProcessOperations(branch_optimization_operations);
    GetEngine()->ResetLogMarginalLikelihood();
    ProcessOperations(marginal_lik_operations);
    double marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
    std::cout << "Current marginal log likelihood: ";
    std::cout << std::setprecision(9) << current_marginal_log_lik << std::endl;
    std::cout << "New marginal log likelihood: ";
    std::cout << std::setprecision(9) << marginal_log_lik << std::endl;
    if (marginal_log_lik < current_marginal_log_lik) {
      std::cout << "Marginal log likelihood decreased.\n";
    }
    if (abs(current_marginal_log_lik - marginal_log_lik) < tol) {
      std::cout << "Converged.\n";
      break;
    }
    current_marginal_log_lik = marginal_log_lik;
  }
  std::chrono::duration<double> optimization_duration = now() - t_start;
  std::cout << "\n# Timing Report\n";
  std::cout << "warmup: " << warmup_duration.count() << "s\n";
  std::cout << "initial likelihood: " << initial_likelihood_duration.count() << "s\n";
  std::cout << "optimization: " << optimization_duration.count() << "s or "
            << optimization_duration.count() / 60 << "m\n";
  std::ofstream branch_length_file("_output/branch_lengths.txt");
  branch_length_file << GetEngine()->GetBranchLengths();
}

void GPInstance::EstimateSBNParameters() {
  std::cout << "Begin SBN parameter optimization\n";
  GPOperationVector sbn_param_optimization_operations = dag_.OptimizeSBNParameters();

  ResetMarginalLikelihoodAndPopulatePLVs();
  ComputeLikelihoods();

  ProcessOperations(sbn_param_optimization_operations);
  sbn_parameters_ = engine_->GetSBNParameters();
}

size_t GPInstance::GetGPCSPIndexForLeafNode(const Bitset &parent_subsplit,
                                            const Node *leaf_node) {
  Assert(leaf_node->IsLeaf(), "Only leaf node is permitted.");
  return dag_.GetGPCSPIndex(parent_subsplit, Bitset::FakeSubsplit(leaf_node->Leaves()));
}

RootedTreeCollection GPInstance::GenerateCompleteRootedTreeCollection() {
  RootedTree::RootedTreeVector tree_vector;
  Node::NodePtrVec topologies = dag_.GenerateAllGPNodeIndexedTopologies();
  const EigenVectorXd gpcsp_indexed_branch_lengths = engine_->GetBranchLengths();

  // Construct Node pointer to parent subsplit encoding.
  // We will use this indexer to look up the GPCSP index.
  std::unordered_map<const Node *, Bitset> node_to_subsplit_indexer;
  for (const auto &topology : topologies) {
    topology->PreOrder([this, &node_to_subsplit_indexer](const Node *node) {
      GPDAGNode *dag_node = dag_.GetDagNode(node->Id());
      if (node_to_subsplit_indexer.count(node) == 0) {
        SafeInsert(node_to_subsplit_indexer, node, dag_node->GetBitset());
      }
    });
  }

  for (auto &root_node : topologies) {
    // Polish will re-assign the node Ids.
    root_node->Polish();

    size_t node_count = 2 * root_node->LeafCount() - 1;
    std::vector<double> branch_lengths(node_count);

    root_node->RootedPCSPPreOrder([this, &branch_lengths,
                                   &gpcsp_indexed_branch_lengths](
                                      const Node *sister, const Node *focal,
                                      const Node *child0, const Node *child1) {
      Bitset parent_subsplit = sister->Leaves() + focal->Leaves();
      Bitset child_subsplit = child0->Leaves() + child1->Leaves();
      size_t gpcsp_idx = dag_.GetGPCSPIndex(parent_subsplit, child_subsplit);
      branch_lengths[focal->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];

      if (sister->IsLeaf()) {
        gpcsp_idx = GetGPCSPIndexForLeafNode(parent_subsplit.RotateSubsplit(), sister);
        branch_lengths[sister->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
      }
      if (child0->IsLeaf()) {
        gpcsp_idx = GetGPCSPIndexForLeafNode(child_subsplit.RotateSubsplit(), child0);
        branch_lengths[child0->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
      }
      if (child1->IsLeaf()) {
        gpcsp_idx = GetGPCSPIndexForLeafNode(child_subsplit, child1);
        branch_lengths[child1->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
      }
    });

    tree_vector.emplace_back(root_node, std::move(branch_lengths));
  }

  RootedTreeCollection rooted_tree_collection(tree_vector,
                                              tree_collection_.TagTaxonMap());
  return rooted_tree_collection;
}

StringVector GPInstance::PrettyIndexer() const {
  StringVector pretty_representation(dag_.GetGPCSPIndexer().size());
  for (const auto &[key, idx] : dag_.GetGPCSPIndexer()) {
    if (idx < dag_.RootsplitCount()) {
      // We have decided to keep the "expanded" rootsplit representation for the time
      // being (see #273). Here we convert it to the representation used in the rest of
      // libsbn.
      auto classic_rootsplit_representation =
          std::min(key.SubsplitChunk(0), key.SubsplitChunk(1));
      pretty_representation[idx] = classic_rootsplit_representation.ToString();
    } else {
      pretty_representation[idx] = key.PCSPToString();
    }
  }
  return pretty_representation;
}

void GPInstance::ProbabilityNormalizeSBNParametersInLog(
    EigenVectorXdRef sbn_parameters) const {
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters, dag_.RootsplitCount(),
                                                  dag_.GetParentToRange());
}

EigenVectorXd GPInstance::NormalizedSBNParametersIncludingFakeSubsplits() const {
  // Extend SBN parameters to include the fake subsplits (which are indexed by the
  // largest set of indices after running
  // AddFakeSubsplitsToGPCSPIndexerAndParentToRange).
  // They don't have any splitting so we set them to zero in log space.
  EigenVectorXd sbn_parameters_result =
      EigenVectorXd::Zero(dag_.GPCSPCountWithFakeSubsplits());
  sbn_parameters_result.segment(0, sbn_parameters_.size()) = sbn_parameters_;
  ProbabilityNormalizeSBNParametersInLog(sbn_parameters_result);
  NumericalUtils::Exponentiate(sbn_parameters_result);
  return sbn_parameters_result;
}

StringDoubleVector GPInstance::PrettyIndexedSBNParameters() {
  StringDoubleVector result;
  auto sbn_parameters = NormalizedSBNParametersIncludingFakeSubsplits();
  result.reserve(sbn_parameters.size());
  const auto pretty_indexer = PrettyIndexer();
  for (size_t i = 0; i < sbn_parameters.size(); i++) {
    result.push_back({pretty_indexer.at(i), sbn_parameters(i)});
  }
  return result;
}

void GPInstance::SBNParametersToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(PrettyIndexedSBNParameters(), file_path);
}
