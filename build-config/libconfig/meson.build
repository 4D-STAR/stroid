config_p = subproject('libconfig',
default_options:[
  'default_library=static',
  'pkg_config=' + get_option('pkg_config').to_string(),
  'build_tests=' + get_option('build_tests').to_string(),
  'build_examples=false'
])
config_dep  = config_p.get_variable('config_dep')
