// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag.hpp"
#include "combinatorics.hpp"
#include "numerical_utils.hpp"
#include "sbn_probability.hpp"

SubsplitDAG::SubsplitDAG()
    : taxon_count_(0), gpcsp_count_without_fake_subsplits_(0), topology_count_(0.) {}

SubsplitDAG::SubsplitDAG(size_t taxon_count,
                         const Node::TopologyCounter &topology_counter)
    : taxon_count_(taxon_count) {
  Assert(topology_counter.size() > 0, "Empty topology counter given to SubsplitDAG.");
  Assert(topology_counter.begin()->first->LeafCount() == taxon_count,
         "Taxon count mismatch in SubsplitDAG constructor.");
  auto [gpcsp_indexer, index_to_child] = ProcessTopologyCounter(topology_counter);
  BuildNodes(index_to_child);
  BuildEdges(index_to_child);
  BuildDAGEdgesFromGPCSPIndexer(gpcsp_indexer);
  AddFakeSubsplitsToDAGEdgesAndParentToRange();
  CountTopologies();
}

SubsplitDAG::SubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection.TaxonCount(), tree_collection.TopologyCounter()) {}

void SubsplitDAG::CountTopologies() {
  topology_count_below_ = EigenVectorXd::Ones(NodeCount());
  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDAGNode(node_id);
    if (!node->IsLeaf()) {
      for (const bool rotated : {false, true}) {
        double per_rotated_count = 0.;
        // Sum options across the possible children.
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id];
        }
        // Take the product across the number of options for the left and right branches
        // of the tree.
        topology_count_below_[node_id] *= per_rotated_count;
      }
    }
  }
  topology_count_ = 0;
  IterateOverRootsplitIds(
      [this](size_t root_id) { topology_count_ += topology_count_below_[root_id]; });
}

size_t SubsplitDAG::NodeCount() const {
  return dag_nodes_.size() - 1;
}  // replace -1 #273

double SubsplitDAG::TopologyCount() const { return topology_count_; }

size_t SubsplitDAG::RootsplitCount() const { return rootsplits_.size(); }

size_t SubsplitDAG::GPCSPCount() const { return gpcsp_count_without_fake_subsplits_; }

size_t SubsplitDAG::GPCSPCountWithFakeSubsplits() const { return dag_edges_.size(); }

void SubsplitDAG::Print() const {
  for (const auto &dag_node : dag_nodes_) {
    std::cout << dag_node->ToString() << std::endl;
  }
}

void SubsplitDAG::PrintGPCSPIndexer() const {
  for (const auto &[pcsp, idx] : BuildGPCSPIndexer()) {
    // Since pcsp may be a rootsubsplit, we need to check and get subsplit.
    std::string gpcsp_string = pcsp.size() == (taxon_count_ * 3)
                                   ? pcsp.PCSPToString()
                                   : pcsp.SubsplitToString();
    std::cout << gpcsp_string << ", " << idx << std::endl;
  }
}

std::string SubsplitDAG::ToDot(bool show_index_labels) const {
  std::stringstream string_stream;
  string_stream << "digraph g {\n";
  string_stream << "node [shape=record];\n";
  string_stream << "edge [colorscheme=dark23];\n";
  DepthFirstWithAction(SubsplitDAGTraversalAction(
      // BeforeNode
      [this, &string_stream, &show_index_labels](size_t node_id) {
        auto bs = GetDAGNode(node_id)->GetBitset();
        string_stream << node_id << " [label=\"<f0>"
                      << bs.SubsplitChunk(0).ToIndexSetString() << "|<f1>";
        if (show_index_labels) {
          string_stream << node_id;
        }
        string_stream << "|<f2>" << bs.SubsplitChunk(1).ToIndexSetString() << "\"]\n";
      },
      // AfterNode
      [](size_t node_id) {},
      // BeforeNodeClade
      [](size_t node_id, bool rotated) {},
      // VisitEdge
      [this, &string_stream, &show_index_labels](size_t node_id, size_t child_id,
                                                 bool rotated) {
        if (GetDAGNode(child_id)->IsLeaf()) {
          string_stream << child_id << " [label=\"<f1>" << child_id << "\"]\n";
        }
        string_stream << "\"" << node_id << "\":";
        string_stream << (rotated ? "f0" : "f2");
        string_stream << "->\"";
        string_stream << child_id << "\":f1";
        if (show_index_labels) {
          string_stream << " [label=\"" << GPCSPIndexOfIds(node_id, child_id);
          if (rotated) {
            string_stream << "\", color=1, fontcolor=1]";
          } else {
            string_stream << "\", color=3, fontcolor=3]";
          }
        }
        string_stream << "\n";
      }));
  string_stream << "}";
  return string_stream.str();
}

BitsetSizeMap SubsplitDAG::BuildGPCSPIndexer() const {
  auto gpcsp_indexer = BitsetSizeMap();
  ReversePostorderIndexTraversal([this, &gpcsp_indexer](size_t parent_id, bool rotated,
                                                        size_t child_id,
                                                        size_t gpcsp_idx) {
    const auto parent_subsplit = GetDAGNode(parent_id)->GetBitset(rotated);
    const auto child_subsplit = GetDAGNode(child_id)->GetBitset();
    SafeInsert(gpcsp_indexer, Bitset::PCSPOfPair(parent_subsplit, child_subsplit),
               gpcsp_idx);
  });
  for (const auto &rootsplit : rootsplits_) {
    const auto subsplit = rootsplit + ~rootsplit;
    SafeInsert(gpcsp_indexer, subsplit, GetRootsplitIndex(subsplit));
  }
  return gpcsp_indexer;
}

SubsplitDAGNode *SubsplitDAG::GetDAGNode(const size_t node_id) const {
  return dag_nodes_.at(node_id).get();
}

size_t SubsplitDAG::GetRootsplitIndex(const Bitset &rootsplit) const {
  return RootsplitIndexOfId(subsplit_to_id_.at(rootsplit));
}

size_t SubsplitDAG::RootsplitIndexOfId(size_t root_id) const {
  return dag_edges_.at({root_id, root_id});  // #273
}

size_t SubsplitDAG::GetGPCSPIndex(const Bitset &parent_subsplit,
                                  const Bitset &child_subsplit) const {
  return GPCSPIndexOfIds(subsplit_to_id_.at(parent_subsplit),
                         subsplit_to_id_.at(child_subsplit));
}

size_t SubsplitDAG::GPCSPIndexOfIds(size_t parent_id, size_t child_id) const {
  return dag_edges_.at({parent_id, child_id});
}

EigenVectorXd SubsplitDAG::BuildUniformOnTopologicalSupportPrior() const {
  EigenVectorXd q = EigenVectorXd::Ones(GPCSPCountWithFakeSubsplits());

  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDAGNode(node_id);
    if (!node->IsLeaf()) {
      for (const bool rotated : {false, true}) {
        double per_rotated_count = 0.;
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id];
        }
        for (const auto &child_id : node->GetLeafward(rotated)) {
          size_t gpcsp_idx = GPCSPIndexOfIds(node->Id(), child_id);
          q(gpcsp_idx) = topology_count_below_(child_id) / per_rotated_count;
        }
      }
    }
  }

  IterateOverRootsplitIds([this, &q](size_t root_id) {
    auto gpcsp_idx = RootsplitIndexOfId(root_id);
    q[gpcsp_idx] = topology_count_below_[root_id] / topology_count_;
  });

  return q;
}

Node::NodePtrVec SubsplitDAG::GenerateAllTopologies() const {
  std::vector<Node::NodePtrVec> topology_below(NodeCount());

  auto GetSubtopologies = [&topology_below](const SubsplitDAGNode *node) {
    Node::NodePtrVec rotated_subtopologies, sorted_subtopologies;
    for (const bool rotated : {false, true}) {
      for (const auto &child_id : node->GetLeafward(rotated)) {
        for (const auto &subtopology : topology_below.at(child_id)) {
          rotated ? rotated_subtopologies.push_back(subtopology)
                  : sorted_subtopologies.push_back(subtopology);
        }
      }
    }
    return std::make_pair(rotated_subtopologies, sorted_subtopologies);
  };

  auto MergeTopologies = [](size_t node_id, Node::NodePtrVec &rotated_subtopologies,
                            Node::NodePtrVec &sorted_subtopologies) {
    Node::NodePtrVec topologies;
    for (const auto &rotated_subtopology : rotated_subtopologies) {
      for (const auto &sorted_subtopology : sorted_subtopologies) {
        Node::NodePtr new_topology =
            Node::Join(sorted_subtopology, rotated_subtopology, node_id);
        topologies.push_back(new_topology);
      }
    }
    return topologies;
  };

  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDAGNode(node_id);
    if (node->IsLeaf()) {
      topology_below.at(node_id).push_back(Node::Leaf(node_id));
    } else {
      auto [rotated_topologies, sorted_topologies] = GetSubtopologies(node);
      topology_below[node_id] =
          MergeTopologies(node_id, rotated_topologies, sorted_topologies);
    }
  }

  Node::NodePtrVec topologies;
  IterateOverRootsplitIds([&topologies, &topology_below](size_t root_id) {
    topologies.insert(topologies.end(), topology_below.at(root_id).begin(),
                      topology_below.at(root_id).end());
  });

  Assert(topologies.size() == TopologyCount(),
         "The realized number of topologies does not match the expected count.");

  // We return a deep copy of every Polished topology to avoid loops in the pointer
  // structure. Such loops can create problems when we Polish the topologies one at a
  // time: polishing a second topology can change the numbering of a previous topology.
  // This is checked for in the "GPInstance: GenerateCompleteRootedTreeCollection" test.
  Node::NodePtrVec final_topologies;
  final_topologies.reserve(topologies.size());
  for (auto &topology : topologies) {
    topology->Polish();
    final_topologies.push_back(topology->DeepCopy());
  }

  return final_topologies;
}

EigenVectorXd SubsplitDAG::BuildUniformOnAllTopologiesPrior() const {
  EigenVectorXd result = EigenVectorXd::Zero(GPCSPCountWithFakeSubsplits());
  for (const auto &[parent_child_id, gpcsp_idx] : dag_edges_) {
    const auto &[parent_id, child_id] = parent_child_id;
    size_t child0_taxon_count, child1_taxon_count;
    if (parent_id == child_id) {  // #273
      child0_taxon_count = GetDAGNode(child_id)->GetBitset().SubsplitChunk(0).Count();
      child1_taxon_count = static_cast<size_t>(taxon_count_ - child0_taxon_count);
    } else {
      // "child 0" of a PCSP refers to the clade of the child subsplit whose bitset has
      // a smaller binary representation, and "child 1" is the larger bitset in binary.
      // For a sorted subsplit, child 0 is the second chunk and child 1 is the first
      // chunk. child 0 and child 1 of a PCSP are further described in bitset.hpp.
      child0_taxon_count = GetDAGNode(child_id)->GetBitset().SubsplitChunk(1).Count();
      child1_taxon_count = GetDAGNode(child_id)->GetBitset().SubsplitChunk(0).Count();
    }
    result(gpcsp_idx) = Combinatorics::LogChildSubsplitCountRatio(child0_taxon_count,
                                                                  child1_taxon_count);
  }
  NumericalUtils::Exponentiate(result);
  return result;
}

void SubsplitDAG::IterateOverRealNodes(const NodeLambda &f) const {
  Assert(taxon_count_ < dag_nodes_.size() - 1,
         "No real DAG nodes!");  // replace -1 #273
  for (auto it = dag_nodes_.cbegin() + taxon_count_; it < dag_nodes_.cend() - 1;
       it++) {  // replace -1 #273
    f((*it).get());
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(const SubsplitDAGNode *node, bool rotated,
                                           const NodeLambda &f) const {
  for (const size_t child_id : node->GetLeafward(rotated)) {
    f(GetDAGNode(child_id));
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(const SubsplitDAGNode *node,
                                           const EdgeDestinationLambda &f) const {
  for (bool rotated : {false, true}) {
    for (const size_t child_id : node->GetLeafward(rotated)) {
      f(rotated, GetDAGNode(child_id));
    }
  }
}

void SubsplitDAG::IterateOverLeafwardEdgesAndChildren(
    const SubsplitDAGNode *node, const EdgeAndNodeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *child) {
        f(GPCSPIndexOfIds(node->Id(), child->Id()), rotated, child->Id());
      });
}

void SubsplitDAG::IterateOverRootwardEdges(const SubsplitDAGNode *node,
                                           const EdgeDestinationLambda &f) const {
  if (not node->IsRoot()) {
    for (bool rotated : {false, true}) {
      for (const size_t parent_idx : node->GetRootward(rotated)) {
        f(rotated, GetDAGNode(parent_idx));
      }
    }
  }
}

void SubsplitDAG::IterateOverRootwardEdgesAndParents(const SubsplitDAGNode *node,
                                                     const EdgeAndNodeLambda &f) const {
  IterateOverRootwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *parent) {
        f(GPCSPIndexOfIds(parent->Id(), node->Id()), rotated, parent->Id());
      });
}

void SubsplitDAG::IterateOverRootsplitIds(const std::function<void(size_t)> &f) const {
  for (const auto &rootsplit : rootsplits_) {
    f(subsplit_to_id_.at(rootsplit + ~rootsplit));
  }
}

RootedIndexerRepresentation SubsplitDAG::IndexerRepresentationOf(
    const BitsetSizeMap &indexer, const Node::NodePtr &topology,
    size_t default_index) const {
  return RootedSBNMaps::IndexerRepresentationOf(indexer, topology, default_index);
}

EigenVectorXd SubsplitDAG::UnconditionalNodeProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  EigenVectorXd node_probabilities(NodeCount());
  node_probabilities.setZero();

  IterateOverRootsplitIds(
      [this, &node_probabilities, &normalized_sbn_parameters](size_t rootsplit_id) {
        Assert(node_probabilities[rootsplit_id] == 0.,
               "We have iterated over the same rootsplit multiple times.");
        node_probabilities[rootsplit_id] +=
            normalized_sbn_parameters[RootsplitIndexOfId(rootsplit_id)];
      });

  ReversePostorderIndexTraversal([&node_probabilities, &normalized_sbn_parameters](
                                     const size_t parent_id, const bool,
                                     const size_t child_id, const size_t gpcsp_idx) {
    const double child_probability_given_parent = normalized_sbn_parameters[gpcsp_idx];
    Assert(child_probability_given_parent >= 0. && child_probability_given_parent <= 1.,
           "UnconditionalNodeProbabilities: got an out-of-range probability. Are these "
           "normalized and in linear space?");
    const double parent_probability = node_probabilities[parent_id];
    node_probabilities[child_id] += parent_probability * child_probability_given_parent;
  });

  return node_probabilities;
}

BitsetDoubleMap SubsplitDAG::UnconditionalSubsplitProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  auto node_probabilities = UnconditionalNodeProbabilities(normalized_sbn_parameters);
  BitsetDoubleMap subsplit_probability_map;
  for (size_t node_id = 0; node_id < node_probabilities.size(); node_id++) {
    const auto &subsplit_bitset = GetDAGNode(node_id)->GetBitset();
    if (!subsplit_bitset.SubsplitIsFake()) {
      SafeInsert(subsplit_probability_map, subsplit_bitset,
                 node_probabilities[node_id]);
    }
  }
  return subsplit_probability_map;
}

EigenVectorXd SubsplitDAG::InvertedGPCSPProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters,
    EigenConstVectorXdRef node_probabilities) const {
  EigenVectorXd inverted_probabilities =
      EigenVectorXd(normalized_sbn_parameters.size());
  // The traversal doesn't set the rootsplit probabilities, but those are always 1
  // (there is only one "parent" of a rootsplit).
  inverted_probabilities.setOnes();
  ReversePostorderIndexTraversal(
      [&node_probabilities, &normalized_sbn_parameters, &inverted_probabilities](
          const size_t parent_id, const bool, const size_t child_id,
          const size_t gpcsp_idx) {
        // For a PCSP t -> s:
        inverted_probabilities[gpcsp_idx] =         // P(t|s)
            node_probabilities[parent_id] *         // P(t)
            normalized_sbn_parameters[gpcsp_idx] /  // P(s|t)
            node_probabilities[child_id];           // P(s)
      });
  return inverted_probabilities;
}

std::vector<Bitset> SubsplitDAG::GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                                   const Bitset &subsplit,
                                                   bool include_fake_subsplits) {
  std::vector<Bitset> children_subsplits;
  if (parent_to_range_.count(subsplit) > 0) {
    const auto [start, stop] = parent_to_range_.at(subsplit);
    for (auto idx = start; idx < stop; idx++) {
      children_subsplits.push_back(index_to_child.at(idx));
    }
  } else if (include_fake_subsplits) {
    // This method is designed to be called before calling
    // AddFakeSubsplitsToDAGEdgesAndParentToRange. In that case, if the second chunk
    // of the subsplit is just a single taxon, the subsplit will not map to any value in
    // parent_to_range_.
    //
    // But we still need to create and connect to fake subsplits in the DAG. So, here we
    // make a fake child subsplit.
    children_subsplits.push_back(Bitset::FakeChildSubsplit(subsplit));
  }

  return children_subsplits;
}

std::pair<BitsetSizeMap, SizeBitsetMap> SubsplitDAG::ProcessTopologyCounter(
    const Node::TopologyCounter &topology_counter) {
  BitsetSizeMap gpcsp_indexer;
  SizeBitsetMap index_to_child;
  std::tie(rootsplits_, gpcsp_indexer, index_to_child, parent_to_range_,
           gpcsp_count_without_fake_subsplits_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topology_counter),
                                  RootedSBNMaps::PCSPCounterOf(topology_counter));
  return {gpcsp_indexer, index_to_child};
}

void SubsplitDAG::CreateAndInsertNode(const Bitset &subsplit) {
  size_t id = dag_nodes_.size();
  SafeInsert(subsplit_to_id_, subsplit, id);
  SafeInsert(subsplit_to_id_, subsplit.RotateSubsplit(), id);
  dag_nodes_.push_back(std::make_unique<SubsplitDAGNode>(id, subsplit));
}

void SubsplitDAG::ConnectNodes(const SizeBitsetMap &index_to_child, size_t idx,
                               bool rotated) {
  const auto node = GetDAGNode(idx);
  // Retrieve children subsplits, set edge relation.
  const Bitset subsplit = node->GetBitset(rotated);
  const auto children = GetChildSubsplits(index_to_child, subsplit, true);
  for (const auto &child_subsplit : children) {
    const auto child_node = GetDAGNode(subsplit_to_id_.at(child_subsplit));
    if (rotated) {
      node->AddLeafwardRotated(child_node->Id());
      child_node->AddRootwardRotated(node->Id());
    } else {
      node->AddLeafwardSorted(child_node->Id());
      child_node->AddRootwardSorted(node->Id());
    }
  }
}

void SubsplitDAG::ConnectRootNode() {
  const auto root_node = GetDAGNode(dag_nodes_.size() - 1);  // replace -1 #273
  for (const auto &rootsplit : rootsplits_) {
    const auto child_node = GetDAGNode(subsplit_to_id_.at(rootsplit + ~rootsplit));
    root_node->AddLeafwardSorted(child_node->Id());
    child_node->AddRootwardSorted(root_node->Id());
    child_node->SetRootsplit();
  }
}

void SubsplitDAG::BuildNodesDepthFirst(const SizeBitsetMap &index_to_child,
                                       const Bitset &subsplit,
                                       std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit : GetChildSubsplits(
             index_to_child, PerhapsRotateSubsplit(subsplit, rotated), false)) {
      if (visited_subsplits.count(child_subsplit) == 0) {
        BuildNodesDepthFirst(index_to_child, child_subsplit, visited_subsplits);
      }
    }
  }
  CreateAndInsertNode(subsplit);
}

void SubsplitDAG::BuildNodes(const SizeBitsetMap &index_to_child) {
  std::unordered_set<Bitset> visited_subsplits;

  // We will create fake subsplits and insert to dag_nodes_.
  // These nodes will take IDs in [0, taxon_count_).
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    CreateAndInsertNode(
        Bitset::FakeSubsplit(Bitset::Singleton(taxon_count_, taxon_idx)));
  }
  // We are going to add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (const auto &rootsplit : rootsplits_) {
    const auto subsplit = rootsplit + ~rootsplit;
    BuildNodesDepthFirst(index_to_child, subsplit, visited_subsplits);
  }
  // We will create a root node with bitset 0000|1111.
  Bitset zeros(taxon_count_);
  CreateAndInsertNode(zeros + ~zeros);
}

void SubsplitDAG::BuildEdges(const SizeBitsetMap &index_to_child) {
  for (size_t i = taxon_count_; i < dag_nodes_.size() - 1; i++) {  // replace -1 #273
    ConnectNodes(index_to_child, i, false);
    ConnectNodes(index_to_child, i, true);
  }
  ConnectRootNode();
}

void SubsplitDAG::BuildDAGEdgesFromGPCSPIndexer(BitsetSizeMap &gpcsp_indexer) {
  for (const auto &[gpcsp, index] : gpcsp_indexer) {
    if (gpcsp.size() == 3 * taxon_count_) {
      const auto parent_id = subsplit_to_id_.at(gpcsp.PCSPParent());
      const auto child_id = subsplit_to_id_.at(gpcsp.PCSPChildSubsplit());
      SafeInsert(dag_edges_, {parent_id, child_id}, index);
    } else if (gpcsp.size() == taxon_count_) {
      const auto rootsplit_id = subsplit_to_id_.at(gpcsp + ~gpcsp);
      SafeInsert(dag_edges_, {rootsplit_id, rootsplit_id}, index);  // #273
    } else {
      Failwith("Don't recognize bitset size!");
    }
  }
}

void SubsplitDAG::AddFakeSubsplitsToDAGEdgesAndParentToRange() {
  for (size_t i = 0; i < taxon_count_; i++) {
    const auto current_bitset = dag_nodes_.at(i)->GetBitset();
    IterateOverRootwardEdges(
        GetDAGNode(i),
        [this, current_bitset](const bool rotated, const SubsplitDAGNode *node) {
          SafeInsert(parent_to_range_, node->GetBitset(rotated),
                     {dag_edges_.size(), dag_edges_.size() + 1});
          SafeInsert(dag_edges_, {node->Id(), subsplit_to_id_.at(current_bitset)},
                     dag_edges_.size());
        });
  }
}

void RootwardDepthFirst(size_t id,
                        const std::vector<std::unique_ptr<SubsplitDAGNode>> &dag_nodes,
                        SizeVector &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  std::copy(visited_nodes.begin(), visited_nodes.end(),
            std::ostream_iterator<int>(std::cout, " "));
  std::cout << "inserting " << id << std::endl;
  SafeInsert(visited_nodes, id);
  std::cout << "inserted" << std::endl;
  const auto &node = dag_nodes.at(id);
  if (not node->IsRoot()) {
    for (size_t child_id : node->GetRootwardSorted()) {
      if (visited_nodes.count(child_id) == 0) {
        RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
      }
    }
    for (size_t child_id : node->GetRootwardRotated()) {
      if (visited_nodes.count(child_id) == 0) {
        RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
      }
    }
  }
  visit_order.push_back(id);
}

void LeafwardDepthFirst(size_t id,
                        const std::vector<std::unique_ptr<SubsplitDAGNode>> &dag_nodes,
                        SizeVector &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  SafeInsert(visited_nodes, id);
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardSorted()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardRotated()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

SizeVector SubsplitDAG::LeafwardPassTraversal() const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (size_t leaf_id = 0; leaf_id < taxon_count_; leaf_id++) {
    RootwardDepthFirst(leaf_id, dag_nodes_, visit_order, visited_nodes);
  }
  return visit_order;
}

SizeVector SubsplitDAG::RootwardPassTraversal() const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  IterateOverRootsplitIds([this, &visit_order, &visited_nodes](size_t root_id) {
    LeafwardDepthFirst(root_id, dag_nodes_, visit_order, visited_nodes);
  });
  return visit_order;
}

SizeVector SubsplitDAG::ReversePostorderTraversal() const {
  auto visit_order = RootwardPassTraversal();
  std::reverse(visit_order.begin(), visit_order.end());
  return visit_order;
}

void SubsplitDAG::ReversePostorderIndexTraversal(
    ParentRotationChildEdgeLambda f) const {
  for (const auto node_id : ReversePostorderTraversal()) {
    IterateOverLeafwardEdgesAndChildren(
        GetDAGNode(node_id), [&f, &node_id](const size_t gpcsp_idx, const bool rotated,
                                            const size_t child_id) {
          f(node_id, rotated, child_id, gpcsp_idx);
        });
  }
}

Bitset SubsplitDAG::PerhapsRotateSubsplit(const Bitset &subsplit, bool rotated) {
  return rotated ? subsplit.RotateSubsplit() : subsplit;
}
