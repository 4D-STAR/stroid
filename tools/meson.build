executable('stroid', 'stroid.cpp', dependencies: [stroid_dep, cli11_dep, magic_enum_dep], install: true)

# Opt-in diagnostic driver; deliberately not part of the installed API/tools.
executable('geometry_quality_experiment', 'geometry_quality_experiment.cpp', dependencies: [stroid_dep, cli11_dep], build_by_default: false, install: false)
