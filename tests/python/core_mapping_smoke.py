"""Focused core-mapping smoke checks for a built or staged STROID Python module.

Example:
    python core_mapping_smoke.py --module _stroid --module-dir build/build-python
    python core_mapping_smoke.py --module-dir /path/to/staged/site-packages
"""

import argparse
import importlib
import json
from pathlib import Path
import sys
import tempfile


def element_counts(stroid, mesh):
    result = stroid.stats.ComputeMeshStats(
        mesh, stroid.stats.MeshStatFeatures.ELEMENT_COUNT
    )
    assert not result.errors, result.errors
    assert result.element_counts is not None
    return result.element_counts


def run(stroid):
    default = stroid.config.MeshConfig()
    assert default.core_mapping == "spherified"
    default.core_mapping = "multi_block"
    assert default.core_mapping == "multi_block"
    assert "core_mapping: multi_block" in repr(default)

    summaries = []
    with tempfile.TemporaryDirectory(prefix="stroid-python-smoke-") as output:
        output_path = Path(output)
        config_path = output_path / "multi_block.toml"
        config_path.write_text(
            '[main]\ncore_mapping = "multi_block"\nrefinement_levels = 0\n'
            'order = 3\ninclude_external_domain = false\n'
            'r_core = 0.25\nr_star = 1.0\nr_infinity = 6.0\n'
            'flattening = 0.0\nr_instability = 1e-14\n'
            'core_steepness = 1.0\ncontinuity_order = 2\n'
            'surface_bdr_id = 1\ninf_bdr_id = 2\n'
            'core_id = 1\nenvelope_id = 2\nvacuum_id = 3\n'
            '[main.optimization_methods]\ntmop = false\nsmoothstep = true\n'
        )
        configured_mesh = stroid.GenerateMesh(str(config_path))
        assert configured_mesh.config.core_mapping == "multi_block"
        assert element_counts(stroid, configured_mesh).total == 13

        for mapping in ("spherified", "multi_block"):
            for external in (False, True):
                config = stroid.config.MeshConfig(
                    core_mapping=mapping,
                    refinement_levels=0,
                    order=3,
                    include_external_domain=external,
                    optimization_methods=stroid.config.OptimizationMethods(
                        tmop=False, smoothstep=True
                    ),
                )
                mesh = stroid.GenerateMesh(config)
                assert mesh.has_mesh() and mesh.has_rmesh()
                counts = element_counts(stroid, mesh)
                expected_core = 7 if mapping == "multi_block" else 1
                expected_total = expected_core + 6 + (6 if external else 0)
                assert counts.total == expected_total
                assert counts.core == expected_core
                assert counts.envelope == 6
                assert counts.vacuum == (6 if external else 0)

                path = output_path / f"{mapping}-{external}.smesh"
                stroid.IO.SaveStroidMesh(mesh, str(path), "Python core-mapping smoke")
                loaded = stroid.IO.LoadStroidMesh(str(path))
                assert loaded.config.core_mapping == mapping
                assert element_counts(stroid, loaded).total == expected_total

                stroid.refinement.UniformRefinement(loaded, 1)
                assert loaded.config.core_mapping == mapping
                assert loaded.refinement_levels == 1
                refined_counts = element_counts(stroid, loaded)
                assert refined_counts.total == expected_total * 8
                assert refined_counts.core == expected_core * 8

                if mapping == "spherified":
                    legacy = "\n".join(
                        line for line in path.read_text().splitlines()
                        if not line.startswith("core_mapping:")
                    )
                    legacy_mesh = stroid.IO.ParseStroidMesh(legacy)
                    assert legacy_mesh.config.core_mapping == "spherified"
                    assert element_counts(stroid, legacy_mesh).total == expected_total

                summaries.append({
                    "mapping": mapping,
                    "external": external,
                    "initial_elements": expected_total,
                    "refined_elements": refined_counts.total,
                })

        invalid = stroid.config.MeshConfig(
            core_mapping="unknown", refinement_levels=0
        )
        try:
            stroid.GenerateMesh(invalid)
        except (ValueError, RuntimeError):
            pass
        else:
            raise AssertionError("Unsupported core_mapping was accepted")

    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", default="stroid")
    parser.add_argument("--module-dir", type=Path)
    args = parser.parse_args()
    if args.module_dir is not None:
        sys.path.insert(0, str(args.module_dir.resolve()))
    stroid = importlib.import_module(args.module)
    summaries = run(stroid)
    print(json.dumps({"module": stroid.__file__, "cases": summaries}, indent=2))


if __name__ == "__main__":
    main()
