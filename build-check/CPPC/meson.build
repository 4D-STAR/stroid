cppc = meson.get_compiler('cpp')

if cppc.get_id() == 'clang'
    message('disabling bitwise-instead-of-logical warnings for clang')
    add_project_arguments('-Wno-bitwise-instead-of-logical', language: 'cpp')
endif

if cppc.get_id() == 'gcc'
    message('disabling psabi warnings for gcc')
    add_project_arguments('-Wno-psabi', language: 'cpp')
    
    if (cppc.version().version_compare('<14.0'))
      error('g++ version must be at least 14.0, found ' + cppc.version())
    endif
endif

if not cppc.has_header('print')
    error('C++ standard library header <print> not found. Please ensure your compiler and standard library supports C++23. We have already validated your compiler version so this is likely an issue with your standard library installation.')
endif
if not cppc.has_header('format')
    error('C++ standard library header <format> not found. Please ensure your compiler and standard library supports C++23. We have already validated your compiler version so this is likely an issue with your standard library installation.')
endif
