#include "mfem.hpp"
#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"
#include "stroid/topology/curvilinear.h"

#include <algorithm>
#include <charconv>
#include <cmath>

#include "stroid/version.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdint>
#include <format>
#include <chrono>
#include <string>
#include <string_view>
#include <expected>
#include <stdexcept>
#include <concepts>
#include <limits>
#include <vector>

namespace stroid::IO {

    namespace {
        std::string format_header(const StroidMesh& mesh, const std::string& comment) {
            auto now = std::chrono::system_clock::now();
            version v;

            std::stringstream vs;
            vs << v;

            std::string header = std::format(R"(# STROID MESH
# NOTE: STROID MESH IS A THIN WRAPPER AROUND MFEM's NATIVE MESH FORMAT
# STRUCTURE:
#           - Type                : Serial or Parallel (S for Serial, P for Parallel)
#           - mesh                : the primary computational domain which can be of n order and be h-refined
#           - reference mesh      : a reference, linear order mesh, used to ensure that the primary mesh remains well formed
#           - exterior coordinate : a scalar material coordinate which is zero at the stellar surface and one at infinity
#           - config              : The configuration options initially used to generate the mesh
#           - refinement-levels   : the total number of refinement levels the primary mesh has been subjected too
# NOTE: EACH BLOCK OF DATA IS STORED BETWEEN "BEGIN BLOCK <NAME>\n ... \nEND BLOCK <NAME>
#       PARSING THE UNDERLYING MFEM NATIVE MESH FORMAT CAN BE DONE WITH MFEM'S STREAM READER
#       IF YOU EXTRACT THE RAW CONTENTS BETWEEN THOSE LINES
BEGIN BLOCK HEADER
    MESH_TYPE:{}
    REFINEMENT_LEVELS:{}
    DATE_CREATED:{:%Y-%m-%d}
    COMMENT:{}
    STROID_VERSION:{}
END BLOCK HEADER)",
            mesh.type == MFEM_MESH_TYPE::PARALLEL ? "P" : "S",
            mesh.refinement_levels,
            now,
            comment,
            vs.str(),
            mesh.refinement_levels
            );
        return header;
        }

        std::string format_primary_mesh(const StroidMesh& mesh) {
            std::stringstream ss;
            ss.precision(std::numeric_limits<double>::max_digits10);
            mesh.mesh->Print(ss);

            std::string pmesh = std::format("BEGIN BLOCK PMESH\n{}END BLOCK PMESH", ss.str());

            return pmesh;
        }

        template <typename T>
        std::string format_opt(const std::optional<T> opt, T default_val) {
            if (opt.has_value()) {
                return std::format("{}", opt.value());
            }
            return std::format("{}", default_val);

        }

        std::string format_reference_mesh(const StroidMesh& mesh) {
            std::stringstream ss;
            ss.precision(std::numeric_limits<double>::max_digits10);
            mesh.reference_mesh->Print(ss);

            std::string rmesh = std::format("BEGIN BLOCK RMESH\n{}END BLOCK RMESH", ss.str());
            return rmesh;
        }

        std::string format_config(const StroidMesh& mesh) {
            config::MeshConfig d;

            config::OptimizationMethods d_opt = d.optimization_methods.value_or(config::OptimizationMethods{false, true});
            config::OptimizationMethods m_opt = mesh.config.optimization_methods.value_or(d_opt);

            std::string config_str = std::format(R"(BEGIN BLOCK CONFIG
# refiniment_levels: Initial number of levels of refinmenet, note the value in the header may be more up to date
# std::optional<int>
# default: 4
refinement_levels:{}

# order: Polynomial / geometric order to use when constructing the mesh
# std::optional<int>
# default: 3
order:{}

# include_external_domain: Whether or not to include the external domain in the mesh generally used for applying boundary conditions at infinity
# std::optional<bool>
# default: true
include_external_domain:{}

# r_core: the radius of the stellar core region (in reference space)
# std::optional<double>
# default: 0.25
r_core:{}

# r_star: the radius of the stellar surface (in reference space)
# std::optional<double>
# default: 1.0
r_star:{}

# flattening: the flattening of the star (in reference space) where 0 is spherical and >0 is oblate. Note that this parameter is not equivalent to solving for the structure of a rotating model
# std::optional<float>
# default: 0.0
flattening:{}

# r_infinity: the radius of the outer boundary of the mesh (in reference space)
# std::optional<double>
# default: 6.0
r_infinity:{}

# r_instability: the radius inside which computations of geometry are skipped to avoid a core singularity
# std::optional<double>
# default: 1e-14
r_instability:{}

# core_steepness: Controls the rate of transition of the core-to-envelope transition
# std::optional<double>
# default: 1.0
core_steepness:{}

# continuity_order: order of continuity to force from teh core-envelope transition (0 = discontinuous, 1=C1 continuity, etc...)
# std::optional<double>
# default: 2
continuity_order:{}

# surface_bdr_id: the boundary id to tag the stellar surface boundary elements as
# std::optional<size_t>
# default: 1
surface_bdr_id:{}

# inf_bdr_id: the boundary id to tag the outer boundary elements as
# std::optional<size_t>
# default: 2
inf_bdr_id:{}

# core_id: the material attribute to tag elements in the core region as
# std::optional<size_t>
# default 1
core_id:{}

# envelope_id: the material attribute to tag elements in the envelope as
# std::optional<size_t>
# default 2
envelope_id:{}

# vacuum_id: the material attribute to tag elements in the vacuum region as
# std::optional<size_t>
# default 3
vacuum_id:{}

# optimization_method: struct for storing which optimization methods are being used
# includes tmop and smoothstep booleans
optimization_methods-tmop:{}
optimization_methods-smoothstep:{}

# core_mapping: Core mapping strategy, either spherified or multi_block
# std::optional<std::string>
# default: spherified
core_mapping:{}
END BLOCK CONFIG)",
            format_opt(mesh.config.refinement_levels, d.refinement_levels.value()),
            format_opt(mesh.config.order, d.order.value()),
            format_opt(mesh.config.include_external_domain, d.include_external_domain.value()),
            format_opt(mesh.config.r_core, d.r_core.value()),
            format_opt(mesh.config.r_star, d.r_star.value()),
            format_opt(mesh.config.flattening, d.flattening.value()),
            format_opt(mesh.config.r_infinity, d.r_infinity.value()),
            format_opt(mesh.config.r_instability, d.r_instability.value()),
            format_opt(mesh.config.core_steepness, d.core_steepness.value()),
            format_opt(mesh.config.continuity_order, d.continuity_order.value()),
            format_opt(mesh.config.surface_bdr_id, d.surface_bdr_id.value()),
            format_opt(mesh.config.inf_bdr_id, d.inf_bdr_id.value()),
            format_opt(mesh.config.core_id, d.core_id.value()),
            format_opt(mesh.config.envelope_id, d.envelope_id.value()),
            format_opt(mesh.config.vacuum_id, d.vacuum_id.value()),
            m_opt.tmop.value_or(false),
            m_opt.smoothstep.value_or(true),
            format_opt(mesh.config.core_mapping, d.core_mapping.value()));

            return config_str;
        }

        std::string format_exterior_coordinate(const StroidMesh& mesh) {
            const bool include_external_domain = mesh.config.include_external_domain.value_or(true);

            if (!include_external_domain) {
                if (mesh.exterior_coordinate) throw std::runtime_error("A mesh without an external domain cannot contain an exterior-coordinate field.");
                return "BEGIN BLOCK EXTERIOR_COORDINATE\nPRESENT:false\nEND BLOCK EXTERIOR_COORDINATE";
            }

            if (!mesh.exterior_coordinate || !mesh.exterior_coordinate->space || !mesh.exterior_coordinate->values) {
                throw std::runtime_error("A mesh with an external domain must contain a complete exterior-coordinate field before it can be saved.");
            }
            if (mesh.exterior_coordinate->space->GetMesh() != mesh.mesh.get()) {
                throw std::runtime_error("The exterior-coordinate finite-element space is attached to the wrong mesh.");
            }
            if (mesh.exterior_coordinate->values->FESpace() != mesh.exterior_coordinate->space.get()) {
                throw std::runtime_error("The exterior-coordinate grid function is attached to the wrong finite-element space.");
            }

            const int scalar_dofs = mesh.exterior_coordinate->space->GetNDofs();
            if (mesh.exterior_coordinate->values->Size() != scalar_dofs) {
                throw std::runtime_error("The exterior-coordinate grid function has an invalid size.");
            }

            std::stringstream ss;
            ss << std::setprecision(std::numeric_limits<double>::max_digits10);
            ss << "BEGIN BLOCK EXTERIOR_COORDINATE\n";
            ss << "PRESENT:true\n";
            ss << "NDOFS:" << scalar_dofs << '\n';
            ss << "VALUES:\n";

            for (int dof = 0; dof < scalar_dofs; ++dof) {
                const double coordinate = (*mesh.exterior_coordinate->values)(dof);
                if (!std::isfinite(coordinate) || coordinate < 0.0 || coordinate > 1.0) {
                    throw std::runtime_error(std::format("Exterior-coordinate DOF {} has invalid value {}.", dof, coordinate));
                }
                ss << coordinate << '\n';
            }

            ss << "END BLOCK EXTERIOR_COORDINATE";
            return ss.str();
        }
    }

    namespace {
        constexpr std::string_view BEGIN_PREFIX = "BEGIN BLOCK ";
        constexpr std::string_view END_PREFIX   = "END BLOCK ";

        std::string_view trim(std::string_view s) {
            const auto b = s.find_first_not_of(" \t\r\n");
            if (b == std::string_view::npos) return {};
            const auto e = s.find_last_not_of(" \t\r\n");
            return s.substr(b, e - b + 1);
        }


        std::expected<bool, std::string> parse_bool(std::string_view v) {
            std::string s(trim(v));
            std::ranges::transform(s, s.begin(),
                                   [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (s == "true"  || s == "1") return true;
            if (s == "false" || s == "0") return false;
            return std::unexpected(std::format("invalid bool value '{}'", v));
        }

        template <std::integral T>
        std::expected<T, std::string> parse_int(std::string_view v) {
            const std::string_view s = trim(v);

            std::string temp(s);

            try {
                size_t pos = 0;

                if constexpr (std::is_signed_v<T>) {
                    long long val = std::stoll(temp, &pos);

                    if (pos != temp.size() || val < std::numeric_limits<T>::min() || val > std::numeric_limits<T>::max()) {
                        return std::unexpected(std::format("invalid integer value '{}'", v));
                    }
                    return static_cast<T>(val);

                } else {
                    unsigned long long val = std::stoull(temp, &pos);

                    if (pos != temp.size() || val > std::numeric_limits<T>::max()) {
                        return std::unexpected(std::format("invalid integer value '{}'", v));
                    }
                    return static_cast<T>(val);
                }
            } catch (const std::exception&) {
                return std::unexpected(std::format("invalid integer value '{}'", v));
            }
        }

        std::expected<double, std::string> parse_double(std::string_view v) {
            const std::string_view s = trim(v);
            std::string temp(s);

            try {
                size_t pos = 0;
                double out = std::stod(temp, &pos);

                if (pos != temp.size()) {
                    return std::unexpected(std::format("invalid floating-point value '{}'", v));
                }
                return out;

            } catch (const std::exception&) {
                return std::unexpected(std::format("invalid floating-point value '{}'", v));
            }
        }

        std::expected<std::map<std::string, std::string>, std::string> extract_blocks(std::istream& is) {
            std::map<std::string, std::string> blocks;
            std::string line;
            std::string current;
            std::string buffer;
            bool in_block = false;

            while (std::getline(is, line)) {
                const std::string_view t = trim(line);

                if (!in_block) {
                    if (t.starts_with(BEGIN_PREFIX)) {
                        current = std::string(trim(t.substr(BEGIN_PREFIX.size())));
                        if (current.empty())
                            return std::unexpected("found 'BEGIN BLOCK' with no block name");
                        if (blocks.contains(current))
                            return std::unexpected(std::format("duplicate block '{}'", current));
                        buffer.clear();
                        in_block = true;
                    }
                } else {
                    if (t.starts_with(END_PREFIX)) {
                        if (const std::string end_name(trim(t.substr(END_PREFIX.size()))); end_name != current)
                            return std::unexpected(std::format(
                                "mismatched block markers: opened '{}' but closed '{}'",
                                current, end_name));
                        blocks.emplace(std::move(current), std::move(buffer));
                        current.clear();
                        buffer.clear();
                        in_block = false;
                    } else {
                        std::string_view raw = line;
                        if (!raw.empty() && raw.back() == '\r') raw.remove_suffix(1);
                        buffer.append(raw);
                        buffer.push_back('\n');
                    }
                }
            }

            if (in_block)
                return std::unexpected(std::format("unterminated block '{}' (missing END BLOCK)", current));
            return blocks;
        }

        std::expected<void, std::string> parse_header(const std::string& content, StroidMesh& out) {
            std::istringstream iss(content);
            std::string line;
            std::optional<MFEM_MESH_TYPE> type;
            std::optional<size_t> ref_levels;

            while (std::getline(iss, line)) {
                const std::string_view t = trim(line);
                if (t.empty() || t.starts_with('#')) continue;

                const auto colon = t.find(':');
                if (colon == std::string_view::npos) continue;

                const std::string_view key = trim(t.substr(0, colon));
                const std::string_view val = trim(t.substr(colon + 1));

                if (key == "MESH_TYPE") {
                    if      (val == "P") type = MFEM_MESH_TYPE::PARALLEL;
                    else if (val == "S") type = MFEM_MESH_TYPE::SERIAL;
                    else return std::unexpected(std::format("unknown MESH_TYPE '{}'", val));
                } else if (key == "REFINEMENT_LEVELS") {
                    auto r = parse_int<size_t>(val);
                    if (!r) return std::unexpected("REFINEMENT_LEVELS: " + r.error());
                    ref_levels = *r;
                }
            }

            if (!type) return std::unexpected("HEADER block missing MESH_TYPE");
            out.type = *type;
            out.refinement_levels = ref_levels.value_or(0);
            return {};
        }

        std::expected<config::MeshConfig, std::string> parse_config(const std::string& content) {
            config::MeshConfig cfg;
            config::OptimizationMethods opt =
                cfg.optimization_methods.value_or(config::OptimizationMethods{});

            using Handler = std::function<std::expected<void, std::string>(std::string_view)>;

            auto as_int    = [](std::optional<int>* f)    { return [f](const std::string_view v) -> std::expected<void, std::string> { auto r = parse_int<int>(v);    if (!r) return std::unexpected(r.error()); *f = *r; return {}; }; };
            auto as_size   = [](std::optional<size_t>* f) { return [f](const std::string_view v) -> std::expected<void, std::string> { auto r = parse_int<size_t>(v); if (!r) return std::unexpected(r.error()); *f = *r; return {}; }; };
            auto as_double = [](std::optional<double>* f) { return [f](const std::string_view v) -> std::expected<void, std::string> { auto r = parse_double(v);     if (!r) return std::unexpected(r.error()); *f = *r; return {}; }; };
            auto as_bool   = [](std::optional<bool>* f)   { return [f](const std::string_view v) -> std::expected<void, std::string> { auto r = parse_bool(v);       if (!r) return std::unexpected(r.error()); *f = *r; return {}; }; };
            auto as_string = [](std::optional<std::string>* f) { return [f](const std::string_view v) -> std::expected<void, std::string> { *f = std::string(v); return {}; }; };

            const std::unordered_map<std::string_view, Handler> handlers = {
                {"refinement_levels",               as_int(&cfg.refinement_levels)},
                {"order",                           as_int(&cfg.order)},
                {"include_external_domain",         as_bool(&cfg.include_external_domain)},
                {"r_core",                          as_double(&cfg.r_core)},
                {"r_star",                          as_double(&cfg.r_star)},
                {"flattening",                      as_double(&cfg.flattening)},
                {"r_infinity",                      as_double(&cfg.r_infinity)},
                {"r_instability",                   as_double(&cfg.r_instability)},
                {"core_steepness",                  as_double(&cfg.core_steepness)},
                {"continuity_order",                as_size(&cfg.continuity_order)},
                {"surface_bdr_id",                  as_size(&cfg.surface_bdr_id)},
                {"inf_bdr_id",                      as_size(&cfg.inf_bdr_id)},
                {"core_id",                         as_size(&cfg.core_id)},
                {"envelope_id",                     as_size(&cfg.envelope_id)},
                {"vacuum_id",                       as_size(&cfg.vacuum_id)},
                {"optimization_methods-tmop",       as_bool(&opt.tmop)},
                {"optimization_methods-smoothstep", as_bool(&opt.smoothstep)},
                {"core_mapping",                    as_string(&cfg.core_mapping)},
            };

            std::istringstream iss(content);
            std::string line;
            while (std::getline(iss, line)) {
                const std::string_view t = trim(line);
                if (t.empty() || t.starts_with('#')) continue;

                const auto colon = t.find(':');
                if (colon == std::string_view::npos) continue;

                const std::string_view key = trim(t.substr(0, colon));
                const std::string_view val = trim(t.substr(colon + 1));

                const auto it = handlers.find(key);
                if (it == handlers.end()) continue;
                if (auto r = it->second(val); !r)
                    return std::unexpected(std::format("{}: {}", key, r.error()));
            }

            cfg.optimization_methods = opt;
            return cfg;
        }

        std::expected<std::unique_ptr<mfem::Mesh>, std::string> load_serial_mesh(const std::string& raw) {
            if (trim(raw).empty()) return std::unexpected("empty mesh block");
            std::istringstream iss(raw);
            try {
                return std::make_unique<mfem::Mesh>(iss);
            } catch (const std::exception& e) {
                return std::unexpected(std::string("MFEM failed to parse mesh: ") + e.what());
            }
        }

        struct ParsedMeta {
            StroidMesh  mesh;
            std::string pmesh_raw;
            std::string rmesh_raw;
            std::optional<std::string> exterior_coordinate_raw;
        };

        struct ParsedExteriorCoordinate {
            bool present{false};
            int scalar_dofs{0};
            std::vector<double> values;
        };

        std::expected<ParsedExteriorCoordinate, std::string> parse_exterior_coordinate(const std::string& content) {
            ParsedExteriorCoordinate parsed;
            std::optional<bool> present;
            std::optional<int> scalar_dofs;
            bool reading_values = false;
            std::istringstream iss(content);
            std::string line;

            while (std::getline(iss, line)) {
                const std::string_view value = trim(line);
                if (value.empty() || value.starts_with('#')) continue;

                if (reading_values) {
                    auto coordinate = parse_double(value);
                    if (!coordinate) return std::unexpected("EXTERIOR_COORDINATE value -> " + coordinate.error());
                    parsed.values.push_back(*coordinate);
                    continue;
                }

                const auto colon = value.find(':');
                if (colon == std::string_view::npos) return std::unexpected(std::format("invalid EXTERIOR_COORDINATE line '{}'.", value));

                const std::string_view key = trim(value.substr(0, colon));
                const std::string_view field_value = trim(value.substr(colon + 1));

                if (key == "PRESENT") {
                    auto result = parse_bool(field_value);
                    if (!result) return std::unexpected("EXTERIOR_COORDINATE PRESENT -> " + result.error());
                    present = *result;
                } else if (key == "NDOFS") {
                    auto result = parse_int<int>(field_value);
                    if (!result) return std::unexpected("EXTERIOR_COORDINATE NDOFS -> " + result.error());
                    scalar_dofs = *result;
                } else if (key == "VALUES") {
                    if (!field_value.empty()) return std::unexpected("EXTERIOR_COORDINATE VALUES must not contain an inline value.");
                    reading_values = true;
                } else {
                    return std::unexpected(std::format("unknown EXTERIOR_COORDINATE key '{}'.", key));
                }
            }

            if (!present.has_value()) return std::unexpected("EXTERIOR_COORDINATE block is missing PRESENT.");
            parsed.present = *present;

            if (!parsed.present) {
                if (scalar_dofs.has_value() || !parsed.values.empty()) return std::unexpected("An absent exterior coordinate cannot contain NDOFS or VALUES.");
                return parsed;
            }

            if (!scalar_dofs.has_value() || *scalar_dofs < 0) return std::unexpected("EXTERIOR_COORDINATE block has an invalid or missing NDOFS.");
            if (static_cast<int>(parsed.values.size()) != *scalar_dofs) {
                return std::unexpected(std::format("EXTERIOR_COORDINATE expected {} values but found {}.", *scalar_dofs, parsed.values.size()));
            }

            parsed.scalar_dofs = *scalar_dofs;
            return parsed;
        }

        std::expected<void, std::string> restore_exterior_coordinate(StroidMesh& mesh, const std::optional<std::string>& raw) {
            fourdst::config::Config<config::MeshConfig> config;
            config.mutate([&mesh](config::MeshConfig& value) { value = mesh.config; });

            try {
                mesh.exterior_coordinate = topology::BuildExteriorCoordinate(*mesh.mesh, *mesh.reference_mesh, config);
            } catch (const std::exception& exception) {
                return std::unexpected(std::string("failed to reconstruct exterior coordinate: ") + exception.what());
            }

            if (!raw.has_value()) return {};

            auto parsed = parse_exterior_coordinate(*raw);
            if (!parsed) return std::unexpected(parsed.error());

            const bool include_external_domain = mesh.config.include_external_domain.value_or(true);
            if (!parsed->present) {
                if (include_external_domain) return std::unexpected("EXTERIOR_COORDINATE is absent even though the mesh includes an external domain.");
                if (mesh.exterior_coordinate) return std::unexpected("An exterior-coordinate field was reconstructed for a mesh without an external domain.");
                return {};
            }

            if (!include_external_domain) return std::unexpected("EXTERIOR_COORDINATE is present for a mesh without an external domain.");
            if (!mesh.exterior_coordinate || !mesh.exterior_coordinate->space || !mesh.exterior_coordinate->values) {
                return std::unexpected("Unable to allocate the exterior-coordinate field while loading the mesh.");
            }
            if (parsed->scalar_dofs != mesh.exterior_coordinate->space->GetNDofs()) {
                return std::unexpected(std::format("EXTERIOR_COORDINATE contains {} DOFs but the reconstructed space has {}.", parsed->scalar_dofs, mesh.exterior_coordinate->space->GetNDofs()));
            }

            constexpr double consistency_tolerance = 1.0e-12;

            for (int dof = 0; dof < parsed->scalar_dofs; ++dof) {
                const double stored_coordinate = parsed->values[static_cast<size_t>(dof)];
                const double reconstructed_coordinate = (*mesh.exterior_coordinate->values)(dof);

                if (!std::isfinite(stored_coordinate) || stored_coordinate < 0.0 || stored_coordinate > 1.0) {
                    return std::unexpected(std::format("EXTERIOR_COORDINATE DOF {} has invalid stored value {}.", dof, stored_coordinate));
                }
                if (std::abs(stored_coordinate - reconstructed_coordinate) > consistency_tolerance) {
                    return std::unexpected(std::format("EXTERIOR_COORDINATE DOF {} is inconsistent with the reference mesh: stored value {}, reconstructed value {}.", dof, stored_coordinate, reconstructed_coordinate));
                }

                (*mesh.exterior_coordinate->values)(dof) = stored_coordinate;
            }

            return {};
        }

        std::expected<ParsedMeta, std::string> parse_metadata(std::istream& is) {
            auto blocks = extract_blocks(is);
            if (!blocks) return std::unexpected(blocks.error());

            auto need = [&](std::string_view name) -> std::expected<std::string, std::string> {
                const auto it = blocks->find(std::string(name));
                if (it == blocks->end())
                    return std::unexpected(std::format("missing required block '{}'", name));
                return it->second;
            };

            ParsedMeta pm{};

            const auto header = need("HEADER");
            if (!header) return std::unexpected(header.error());
            if (auto r = parse_header(*header, pm.mesh); !r) return std::unexpected(r.error());

            const auto config = need("CONFIG");
            if (!config) return std::unexpected(config.error());
            auto cfg = parse_config(*config);
            if (!cfg) return std::unexpected("CONFIG block -> " + cfg.error());
            pm.mesh.config = std::move(*cfg);

            const auto pmesh = need("PMESH");
            if (!pmesh) return std::unexpected(pmesh.error());
            pm.pmesh_raw = *pmesh;

            const auto rmesh = need("RMESH");
            if (!rmesh) return std::unexpected(rmesh.error());
            pm.rmesh_raw = *rmesh;

            if (const auto exterior_coordinate = blocks->find("EXTERIOR_COORDINATE"); exterior_coordinate != blocks->end()) {
                pm.exterior_coordinate_raw = exterior_coordinate->second;
            }

            return pm;
        }

    }

    void SaveStroidMesh(const StroidMesh &mesh, const std::string &filename, const std::string &comment) {
        std::ofstream ofs(filename);

        // First Write a header with some information
        std::string header = format_header(mesh, comment);

        std::string pmesh = format_primary_mesh(mesh);
        std::string rmesh = format_reference_mesh(mesh);

        std::string config = format_config(mesh);
        std::string exterior_coordinate = format_exterior_coordinate(mesh);

        ofs << header << "\n";
        ofs << pmesh << "\n";
        ofs << rmesh << "\n";
        ofs << config << "\n";
        ofs << exterior_coordinate << "\n";
    }

    void SaveMesh(const mfem::Mesh& mesh, const std::string& filename) {
        std::ofstream ofs(filename);
        ofs.precision(std::numeric_limits<double>::max_digits10);
        mesh.Print(ofs);
    }

    void SaveMesh(const stroid::StroidMesh &mesh, const std::string &filename) {
        SaveMesh(*mesh.mesh, filename);
    }

    void SaveVTU(mfem::Mesh &mesh, const std::string &exportName) {
        mfem::ParaViewDataCollection pd(exportName, &mesh);
        pd.SetDataFormat(mfem::VTKFormat::BINARY);
        pd.SetHighOrderOutput(true);
        pd.Save();
    }

    void SaveVTU(const stroid::StroidMesh &mesh, const std::string &exportName) {
        SaveVTU(*mesh.mesh, exportName);
    }

    void ViewMesh(mfem::Mesh &mesh, const std::string& title, const VISUALIZATION_MODE mode, const std::string &vishost, int visport) {
        mfem::socketstream sol_sock(vishost.c_str(), visport);
        if (!sol_sock.is_open()) {
            std::cerr << "Unable to connect to GLVis server at "
                      << vishost << ':' << visport << std::endl;
            return;
        }

        mfem::L2_FECollection fec(0, mesh.Dimension());
        mfem::FiniteElementSpace fes(&mesh, &fec);
        mfem::GridFunction attr_gf(&fes);
        attr_gf = 0.0;

        switch (mode) {
            case VISUALIZATION_MODE::ELEMENT_ID:
                for (int i = 0; i < mesh.GetNE(); i++) {
                    attr_gf(i) = static_cast<double>(mesh.GetAttribute(i));
                }
                break;
            case VISUALIZATION_MODE::BOUNDARY_ELEMENT_ID:
                attr_gf = 0.0;
                for (int i = 0; i < mesh.GetNBE(); i++) {
                    int elem_index, side_index;
                    mesh.GetBdrElementAdjacentElement(i, elem_index, side_index);
                    attr_gf(elem_index) = static_cast<double>(mesh.GetBdrAttribute(i));
                }
                break;

            case VISUALIZATION_MODE::NONE:
            default:
                break;
        }

        sol_sock.precision(8);
        sol_sock << "solution\n" << mesh << attr_gf;
        sol_sock << "window_title '" << title << "'\n";
        sol_sock << "keys iMj\n";
        sol_sock << std::flush;
    }

    void ViewMesh(const stroid::StroidMesh &mesh, const std::string &title, VISUALIZATION_MODE mode, const std::string &vishost, int visport) {
        ViewMesh(*mesh.mesh, title, mode, vishost, visport);
    }

    void VisualizeFaceValence(mfem::Mesh& mesh, const std::string &vishost, int visport) {
        mfem::L2_FECollection fec(0, 3);
        mfem::FiniteElementSpace fes(&mesh, &fec);
        mfem::GridFunction valence_gf(&fes);

        for (int i = 0; i < mesh.GetNBE(); i++) {
            int f, o;
            mesh.GetBdrElementFace(i, &f, &o);

            int e1, e2;
            mesh.GetFaceElements(f, &e1, &e2);

            int valence = (e2 >= 0) ? 2 : 1;
            valence_gf(i) = static_cast<double>(valence);
        }

        // View in GLVis
        mfem::socketstream sol_sock(vishost.c_str(), visport);
        if (sol_sock.is_open()) {
            sol_sock << "solution\n" << mesh << valence_gf;
            sol_sock << "window_title 'Boundary Valence: 1=Surface, 2=Internal'\n";
            sol_sock << "keys am\n" << std::flush;
        }
    }

    void VisualizeFaceValence(const stroid::StroidMesh &mesh, const std::string &vishost, int visport) {
        VisualizeFaceValence(*mesh.mesh, vishost, visport);
    }

    std::expected<StroidMesh, std::string> ParseStroidMesh(std::istream& is) {
        auto pm = parse_metadata(is);
        if (!pm) return std::unexpected(pm.error());

        if (pm->mesh.type != MFEM_MESH_TYPE::SERIAL) {
            return std::unexpected(
                "parsed a PARALLEL StroidMesh, but ParseStroidMesh(std::istream&) can only "
                "reconstruct serial meshes; use the MPI-aware overload "
                "ParseStroidMesh(std::istream&, MPI_Comm) (requires MFEM_USE_MPI)");
        }

        auto m = load_serial_mesh(pm->pmesh_raw);
        if (!m) return std::unexpected("PMESH -> " + m.error());
        auto rm = load_serial_mesh(pm->rmesh_raw);
        if (!rm) return std::unexpected("RMESH -> " + rm.error());

        pm->mesh.mesh           = std::move(*m);
        pm->mesh.reference_mesh = std::move(*rm);
        if (auto result = restore_exterior_coordinate(pm->mesh, pm->exterior_coordinate_raw); !result) return std::unexpected(result.error());
        return std::move(pm->mesh);
    }

    std::expected<StroidMesh, std::string> LoadStroidMesh(const std::string& filename) {
        std::ifstream ifs(filename);
        if (!ifs.is_open())
            return std::unexpected(std::format("could not open file '{}'", filename));
        return ParseStroidMesh(ifs);
    }

#ifdef MFEM_USE_MPI
    std::expected<StroidMesh, std::string> ParseStroidMesh(std::istream& is, MPI_Comm comm) {
        auto pm = parse_metadata(is);
        if (!pm) return std::unexpected(pm.error());

        auto build = [&](const std::string& raw)
            -> std::expected<std::unique_ptr<mfem::Mesh>, std::string> {
            if (trim(raw).empty()) return std::unexpected("empty mesh block");
            std::istringstream iss(raw);
            try {
                if (pm->mesh.type == MFEM_MESH_TYPE::PARALLEL)
                    return std::unique_ptr<mfem::Mesh>(new mfem::ParMesh(comm, iss));
                return std::make_unique<mfem::Mesh>(iss);
            } catch (const std::exception& e) {
                return std::unexpected(std::string("MFEM failed to parse mesh: ") + e.what());
            }
        };

        auto m = build(pm->pmesh_raw);
        if (!m) return std::unexpected("PMESH -> " + m.error());
        auto rm = build(pm->rmesh_raw);
        if (!rm) return std::unexpected("RMESH -> " + rm.error());

        pm->mesh.mesh           = std::move(*m);
        pm->mesh.reference_mesh = std::move(*rm);
        if (auto result = restore_exterior_coordinate(pm->mesh, pm->exterior_coordinate_raw); !result) return std::unexpected(result.error());
        return std::move(pm->mesh);
    }

    std::expected<StroidMesh, std::string> LoadStroidMesh(const std::string& filename, MPI_Comm comm) {
        std::ifstream ifs(filename);
        if (!ifs.is_open())
            return std::unexpected(std::format("could not open file '{}'", filename));
        return ParseStroidMesh(ifs, comm);
    }
#endif // MFEM_USE_MPI


}
