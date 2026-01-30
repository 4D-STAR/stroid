#include <memory>
#include "mfem.hpp"

#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"
#include "stroid/topology/curvilinear.h"
#include "stroid/topology/topology.h"

#include "stroid/stroid.h"

#include "fourdst/config/config.h"

#include "CLI/CLI.hpp"
#include <optional>

enum INFO_MODE {
    VERSION,
    DEFAULT_CONFIG
};

int main(int argc, char** argv) {
    fourdst::config::Config<stroid::config::MeshConfig> cfg;

    CLI::App app{"stroid - A tool for generating multi-block meshes for stellar modeling"};
    auto* generate = app.add_subcommand("generate", "Generate a multi-block mesh for stellar modeling");
    auto* info = app.add_subcommand("info", "Access information about the program stroid");

    std::optional<std::string> config_filename;
    std::string output_filename;
    bool view_mesh;
    bool no_save;
    std::optional<std::string> glvis_host;
    std::optional<int> glvis_port;

    generate->add_option("-c,--config", config_filename, "Path to configuration file")->check(CLI::ExistingFile);
    generate->add_flag("-v,--view", view_mesh, "View the generated mesh using GLVis");
    generate->add_flag("-n,--nosave", no_save, "Do not save the generated mesh to a file");
    generate->add_option("--glvis-host", glvis_host, "GLVis server host (default: localhost)")->default_val("localhost");
    generate->add_option("--glvis-port", glvis_port, "GLVis server port (default: 19916)")->default_val(19916);
    generate->add_option("-o,--output", output_filename, "Output filename for the generated mesh")->default_val("stroid");


    info->add_flag_callback("-v,--version", []() {
        std::println("Stroid Version {}", stroid::version::toString());
        return 0;
    }, "Display stroid version information");

    info->add_flag_callback("-d,--default", [&cfg]() {
        cfg.save("default.toml");
    }, "Save the default configuration to a file default.toml");


    CLI11_PARSE(app, argc, argv);

    // Ensure that if view is requested, host and port are set
    if (view_mesh) {
        if (!glvis_host.has_value()) {
            glvis_host = "localhost";
        }
        if (!glvis_port.has_value()) {
            glvis_port = 19916;
        }
    }

    // If config filename is provided, load configuration from file
    if (config_filename.has_value()) {
        cfg.load(config_filename.value());
    }


    const std::unique_ptr<mfem::Mesh> mesh = stroid::topology::BuildSkeleton(cfg);
    stroid::topology::Finalize(*mesh, cfg);
    stroid::topology::PromoteToHighOrder(*mesh, cfg);
    stroid::topology::ProjectMesh(*mesh, cfg);


    if (!output_filename.empty() && !no_save) {
        stroid::IO::SaveVTU(*mesh, output_filename);
    }

    if (view_mesh) {
        stroid::IO::ViewMesh(*mesh, "Spheroidal Mesh - Colored by Element ID", stroid::IO::VISUALIZATION_MODE::ELEMENT_ID, glvis_host.value(), glvis_port.value());
    }
}