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
    'lib/topology/optimize.cpp',
    'lib/IO/mesh.cpp',
    'lib/utils/mesh_utils.cpp',
    'lib/utils/mesh_stats.cpp',
    'lib/refinement/uniform.cpp',
)

if get_option('build_python')
    if host_machine.system() == 'darwin'
        stroid_lib_rpath_args = [
            '-Wl,-rpath,@loader_path',
            '-Wl,-rpath,@loader_path/../../fourdst/lib',
            '-Wl,-rpath,@loader_path/../../fourdst/lib/vendor',
        ]
        stroid_lib_rpath = ''

        stroid_ext_rpath_args = [
            '-Wl,-rpath,@loader_path/lib',
            '-Wl,-rpath,@loader_path/../fourdst/lib',
            '-Wl,-rpath,@loader_path/../fourdst/lib/vendor',
        ]
        stroid_ext_rpath = ''
    else
        stroid_lib_rpath_args = []
        stroid_lib_rpath = '$ORIGIN:' + '$ORIGIN/../../fourdst/lib:' + '$ORIGIN/../../fourdst/lib/vendor'

        stroid_ext_rpath_args = []
        stroid_ext_rpath = '$ORIGIN/lib:' + '$ORIGIN/../fourdst/lib:' + '$ORIGIN/../fourdst/lib/vendor'
    endif

    libstroid = static_library(
        'stroid',
        stroid_sources,
        include_directories: stroid_include_files,
        dependencies: dependencies,
        install_dir: stroid_libdir,
        link_args: stroid_lib_rpath_args,
        build_rpath: stroid_lib_rpath,
        install_rpath: stroid_lib_rpath,
    )
else
    libstroid = static_library(
        'stroid',
        stroid_sources,
        include_directories: stroid_include_files,
        dependencies: dependencies,
        install: true
    )
endif


if get_option('build_python')
    stroid_iface_dep = declare_dependency(
        dependencies: dependencies
    ).partial_dependency(compile_args: true, includes: true)

    stroid_dep = declare_dependency(
        link_with: libstroid,
        include_directories: stroid_include_files,
        dependencies: [stroid_iface_dep]
    )
else
    stroid_dep = declare_dependency(
        link_with: libstroid,
        include_directories: stroid_include_files,
        dependencies: dependencies
    )
endif

meson.override_dependency('stroid', stroid_dep)

message('stroid include dir: ' + stroid_includedir)
install_subdir(
    'include/stroid',
    install_dir: stroid_includedir,
    exclude_files: ['version.h.in']
)