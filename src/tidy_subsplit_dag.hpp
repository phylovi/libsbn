// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A "tidy" subsplit DAG has a notion of clean and dirty vectors.

#ifndef SRC_TIDY_SUBSPLIT_DAG_HPP_
#define SRC_TIDY_SUBSPLIT_DAG_HPP_

#include "subsplit_dag.hpp"

class TidySubsplitDAG : public SubsplitDAG {
 public:
  TidySubsplitDAG();
  explicit TidySubsplitDAG(const RootedTreeCollection &tree_collection);
  // This constructor is really just meant for testing.
  explicit TidySubsplitDAG(size_t node_count);

  // What nodes are above or below the specified node? We consider a node to be both
  // above and below itself (this just happens to be handy for the implementation).
  EigenArrayXb BelowNode(size_t node_idx);
  EigenArrayXbRef BelowNode(bool rotated, size_t node_idx);
  EigenArrayXb AboveNode(size_t node_idx);
  EigenArrayXb AboveNode(bool rotated, size_t node_idx);

  // Set the below matrix up to have the node at src_idx below the subsplit-clade
  // described by (dst_rotated, dst_idx).
  void SetBelow(size_t dst_idx, bool dst_rotated, size_t src_idx);

  // (0,(1,(2,3))) and ((0,(2,3)),1)
  // See https://github.com/phylovi/libsbn/issues/307#issuecomment-765901588
  // Update during #288
  static TidySubsplitDAG MotivatingExample();

 private:
  // above_rotated_.(i,j) is true iff i,true is above j.
  EigenMatrixXb above_rotated_;
  // above_sorted.(i,j) is true iff i,false is above j.
  EigenMatrixXb above_sorted_;

  TidySubsplitDAG(size_t taxon_count, const Node::TopologyCounter &topology_counter);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TidySubsplitDAG: slicing") {
  auto dag = TidySubsplitDAG(5);

  // The tree ((0,1)3,2)4:
  dag.SetBelow(3, true, 0);
  dag.SetBelow(3, false, 1);
  dag.SetBelow(4, true, 3);
  dag.SetBelow(4, false, 2);

  CHECK_EQ(GenericToString(dag.AboveNode(0)), "[1, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(dag.AboveNode(1)), "[0, 1, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(dag.AboveNode(2)), "[0, 0, 1, 0, 1]\n");
  CHECK_EQ(GenericToString(dag.AboveNode(3)), "[0, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(dag.AboveNode(4)), "[0, 0, 0, 0, 1]\n");
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TIDY_SUBSPLIT_DAG_HPP_
