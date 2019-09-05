// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

// This is a class implementing the an indexing scheme for the Primary Subsplit
// Pair branch length parameterization.
// See the 2018 ICLR paper for details.

#ifndef SRC_PSP_INDEXER_HPP_
#define SRC_PSP_INDEXER_HPP_

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "sbn_maps.hpp"
#include "sugar.hpp"
#include "tree.hpp"

class PSPIndexer {
 public:
  PSPIndexer() : after_rootsplits_index_(0), first_empty_index_(0) {}
  PSPIndexer(BitsetVector rootsplits, BitsetSizeMap in_indexer);

  size_t AfterRootsplitsIndex() const { return after_rootsplits_index_; }
  size_t FirstEmptyIndex() const { return first_empty_index_; }
  // These are just some things that we may want to know about the indexer.
  StringSizeMap Details() const {
    return {
        {"after_rootsplits_index", after_rootsplits_index_},
        {"first_empty_index", first_empty_index_},
        // This is the "official" definition of a PSP indexer representation of
        // a tree. It's a vector of vectors, where the order of entries of the
        // outer vector is laid out as follows.
        {"rootsplit_position", 0},
        {"subsplit_down_position", 1},
        {"subsplit_up_position", 2},
    };
  }

  // Reverse the indexer to a vector of strings.
  // We add in another extra empty string at the end for "not found."
  StringVector ToStringVector() const;

  // Get the PSP representation of a given topology.
  SizeVectorVector RepresentationOf(const Node::NodePtr& topology) const;
  // Get the string version of the representation.
  // Inefficiently implemented, so for testing only.
  StringVectorVector StringRepresentationOf(
      const Node::NodePtr& topology) const;

 private:
  BitsetSizeMap indexer_;
  size_t after_rootsplits_index_;
  size_t first_empty_index_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("PSPIndexer") {}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_PSP_INDEXER_HPP_