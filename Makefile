default:
	make bison && scons && ./_build/newick_parser _ignore/ex.nwk | paste _ignore/ex.nwk -
	./_build/doctest

bison: src/parser.yy src/scanner.ll
	bison -o src/parser.cpp --defines=src/parser.hpp src/parser.yy
	flex -o src/scanner.cpp src/scanner.ll

format:
	clang-format -i -style=file src/* test/test.c

