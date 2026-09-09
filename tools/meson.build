stroid_cli = executable('stroid', 'stroid.cpp', dependencies: [stroid_dep, cli11_dep, magic_enum_dep], install: true)

executable('geometry_quality_experiment', 'geometry_quality_experiment.cpp', dependencies: [stroid_dep, cli11_dep], build_by_default: false, install: false)

executable('vacuum_refinement_experiment', 'vacuum_refinement_experiment.cpp', dependencies: [stroid_dep, cli11_dep], build_by_default: false, install: false)
