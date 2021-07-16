// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A node in a directed acyclic graph representing a collection of subsplits with their
// corresponding parent-child relationships.
//
// Each node represents a sorted osubsplit, which is stored as a bitset in `subsplit_`.
// The leafward edges are divided into two groups based on if they split apart the
// right clade (in which case they are called sorted children) or if they split apart
// the left clade (in which case they are called rotated children).
// Similarly, the rootward edges are divided into two groups based on if the child of
// the edge splits apart the right clade of the parent (in which case they are called
// sorted parents) or if they split apart the left clade of the parent (in which case
// they are called rotated parents).

#ifndef SRC_SUBSPLIT_DAG_NODE_HPP_
#define SRC_SUBSPLIT_DAG_NODE_HPP_

#include "bitset.hpp"
#include "sugar.hpp"

class SubsplitDAGNode {
 public:
  SubsplitDAGNode(size_t id, Bitset subsplit)
      : id_(id), subsplit_(std::move(subsplit)) {}

  size_t Id() const { return id_; }
  const Bitset &GetBitset() const { return subsplit_; }
  const Bitset GetBitset(bool rotated) const {
    return rotated ? subsplit_.RotateSubsplit() : subsplit_;
  }
  bool IsRootNode() const {
    return (rootward_sorted_.empty() && rootward_rotated_.empty());
  };
  bool IsRootsplit() const {
    return subsplit_.Count() == subsplit_.SubsplitChunkSize() && !IsRootNode();
  }
  bool IsLeaf() const { return leafward_rotated_.empty() && leafward_sorted_.empty(); }

  void AddLeafwardRotated(size_t node_id) { leafward_rotated_.push_back(node_id); }
  void AddLeafwardSorted(size_t node_id) { leafward_sorted_.push_back(node_id); }
  void AddRootwardRotated(size_t node_id) { rootward_rotated_.push_back(node_id); }
  void AddRootwardSorted(size_t node_id) { rootward_sorted_.push_back(node_id); }
  const SizeVector &GetLeafwardRotated() const { return leafward_rotated_; }
  const SizeVector &GetLeafwardSorted() const { return leafward_sorted_; }
  const SizeVector &GetLeafward(bool rotated) const {
    return rotated ? GetLeafwardRotated() : GetLeafwardSorted();
  }
  const SizeVector &GetRootwardRotated() const { return rootward_rotated_; }
  const SizeVector &GetRootwardSorted() const { return rootward_sorted_; }
  const SizeVector &GetRootward(bool rotated) const {
    return rotated ? GetRootwardRotated() : GetRootwardSorted();
  }

  std::string ToString() const;

 private:
  size_t id_;
  Bitset subsplit_;

  SizeVector leafward_rotated_;
  SizeVector leafward_sorted_;
  SizeVector rootward_rotated_;
  SizeVector rootward_sorted_;
};

#endif  // SRC_SUBSPLIT_DAG_NODE_HPP_
