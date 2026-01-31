dependencies = [
    mfem_dep,
    config_dep,
]

subdir('include/stroid')

stroid_include_files = include_directories('include')
stroid_sources = files(
    'lib/topology/curvilinear.cpp',
    'lib/topology/mapping.cpp',
    'lib/topology/topology.cpp',
    'lib/IO/mesh.cpp',
    'lib/utils/mesh_utils.cpp',
)

stroid_lib = static_library(
    'libstroid',
    stroid_sources,
    include_directories: stroid_include_files,
    dependencies: dependencies,
    install: true
)

stroid_dep = declare_dependency(
    link_with: stroid_lib,
    include_directories: stroid_include_files,
    dependencies: dependencies
)

install_subdir(
    'include/stroid',
    install_dir: get_option('includedir') / 'stroid'
)