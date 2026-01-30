project('stroid', 'cpp', meson_version : '>= 1.3.0', version : 'v0.1.0', default_options : ['cpp_std=c++23'])

subdir('build-config')
subdir('src')

if get_option('build_tests')
    subdir('tests')
endif

if get_option('build_tools')
    subdir('tools')
endif

if get_option('pkg_config')
    pkg = import('pkgconfig')
    pkg.generate(
        name: 'stroid',
        description: 'Stroid multi-block curvilinear mesh generation library',
        version: meson.project_version(),
        libraries: [
            stroid_lib
        ],
        subdirs: ['stroid'],
        filebase: 'stroid',
        install_dir: join_paths(get_option('libdir'), 'pkgconfig')
    )
endif