project('bayes','cpp','cython', version : '1.0.0',
        default_options : ['optimization=3'])
add_project_arguments(['-std=c++20'], language : 'cpp')

liblognormal = static_library(
    'lognormal', ['src/lognormal.cpp','src/math.cpp'],
    dependencies : [dependency('boost'), dependency('openmp')],
    include_directories : include_directories('include')
)


python = import('python').find_installation()
dep_py = python.dependency()

python.extension_module('bayes',
                        'pybalonor/bayes.pyx',
                        dependencies : [dep_py, dependency('boost')],
                        include_directories : include_directories('include'),
                        override_options : ['cython_language=cpp'],
                        link_with : liblognormal)
