project('bayes','cpp','cython',
        default_options : ['optimization=3'])
add_project_arguments(['-std=c++20'], language : 'cpp')

liblognormal = static_library(
    'lognormal', ['lognormal.cpp'],
    dependencies : [dependency('boost'), dependency('openmp')],
    include_directories : include_directories('.')
)


python = import('python').find_installation()
dep_py = python.dependency()

python.extension_module('bayes',
                        'bayes.pyx',
                        dependencies : [dep_py],
                        include_directories : include_directories('.'),
                        override_options : ['cython_language=cpp'],
                        link_with : liblognormal)
