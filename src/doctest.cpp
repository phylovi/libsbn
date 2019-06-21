#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "driver.hpp"
#include "sbn.hpp"

TEST_CASE("Node") {

Driver driver;

std::vector<std::string> trace;
auto t = driver.ParseString("((((0_1,1_1),(2_1,3_1)),4_1),((5_1,(6_1,7_1)),(8_1,9_1)));");

// preorder:
t->PreOrder([&trace](Node* node) { trace.push_back(node->TagString()); });
CHECK(std::vector<std::string>({"9_10","4_5","3_4","1_2","0_1","1_1","3_2","2_1","3_1","4_1","9_5","7_3","5_1","7_2","6_1","7_1","9_2","8_1","9_1"}) == trace);
trace.clear();

// postorder:
t->PostOrder([&trace](Node* node) { trace.push_back(node->TagString()); });
CHECK(std::vector<std::string>({"0_1","1_1","1_2","2_1","3_1","3_2","3_4","4_1","4_5","5_1","6_1","7_1","7_2","7_3","8_1","9_1","9_2","9_5","9_10"}) == trace);
trace.clear();

// levelorder:
t->LevelOrder([&trace](Node* node) { trace.push_back(node->TagString()); });
CHECK(std::vector<std::string>({"9_10","4_5","9_5","3_4","4_1","7_3","9_2","1_2","3_2","5_1","7_2","8_1","9_1","0_1","1_1","2_1","3_1","6_1","7_1"}) == trace);
trace.clear();

}
