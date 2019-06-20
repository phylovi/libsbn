import os

env = Environment(ENV=os.environ, CPPPATH=['include', 'src']) #, CPPFLAGS=['-Wall', '-Wextra', '-Werror', '-Wconversion'])

env.VariantDir('_build', 'src')
env.VariantDir('_build/test', 'test')

env.Library('sbn', ['_build/libsbn.cpp', '_build/sbn.cpp'])

env.Program('_build/newick_parser', ['_build/newick_parser.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp'], LIBS=['fl','sbn'], LIBPATH='.')

doctest = env.Program(['_build/test/doctest.cpp', '_build/sbn.cpp'])
#test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.', LINK='g++')
