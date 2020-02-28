// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_SBN_INSTANCE_HPP_
#define SRC_ROOTED_SBN_INSTANCE_HPP_

#include "libsbn.hpp"

// TODO refactor so that we have a shared object from which the current SBNInstance and
// the RootedSBNInstance descend.

class RootedSBNInstance : public SBNInstance {
 public:
  explicit RootedSBNInstance(const std::string &name) : SBNInstance(name) {}

  void UseCurrentTreesAsRooted();
  std::vector<double> LogLikelihoods();
  // For each loaded tree, returns a pair of (likelihood, gradient).
  std::vector<std::pair<double, std::vector<double>>> BranchGradients();

 private:
  RootedTreeCollection rooted_tree_collection_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedSBNInstance") {
  // TODO
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_SBN_INSTANCE_HPP_