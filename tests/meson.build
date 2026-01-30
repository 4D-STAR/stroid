gtest_dep = dependency('gtest', main: true, required : true)
gtest_main = dependency('gtest_main', required: true)
gtest_nomain_dep = dependency('gtest', main: false, required : true)

# Test files for const
test_sources = [
    'stroidTest.cpp'
]

foreach test_file : test_sources
    exe_name = test_file.split('.')[0]

    test_exe = executable(
        exe_name,
        test_file,
        dependencies: [
            stroid_dep,
            gtest_dep,
            gtest_main
        ],
    )

    test(
        exe_name,
        test_exe,
        env: ['MESON_SOURCE_ROOT=' + meson.project_source_root(), 'MESON_BUILD_ROOT=' + meson.project_build_root()])
endforeach
