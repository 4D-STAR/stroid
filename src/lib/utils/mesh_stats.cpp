#include "stroid/utils/mesh_stats.h"

#include <cmath>
#include <limits>
#include <format>
#include <string>
#include <algorithm>

namespace stroid::stats {
    namespace {
        double SpheroidRadius(const double ux, const double uy, const double uz,
                              const double r_star, const double flattening) {
            const double a = r_star;
            const double c = r_star * (1.0 - flattening);
            const double inv = (ux*ux + uy*uy) / (a*a) + (uz*uz) / (c*c);
            return (inv > 0.0) ? 1.0 / std::sqrt(inv) : 0.0;
        }
    }

    MeshStats ComputeMeshStats(const StroidMesh& sm, MeshStatFeatures features, int sample_order) {
        MeshStats out;
        out.computed = features;

        auto mesh_or = sm.as_mesh();
        if (!mesh_or) {
            out.errors.push_back(mesh_or.error());
            return out;
        }
        mfem::Mesh* mesh = *mesh_or;
        if (!mesh) {
            out.errors.push_back("StroidMesh has no stored computational mesh to compute stats against");
            return out;   // BUGFIX: was falling through to a null deref
        }

        const auto& cfg = sm.config;
        const double r_star     = cfg.r_star.value_or(1.0);
        const double flattening = cfg.flattening.value_or(0.0);
        const int    surf_bdr   = static_cast<int>(cfg.surface_bdr_id.value_or(-99));
        const int    inf_bdr    = static_cast<int>(cfg.inf_bdr_id.value_or(-99));
        const int    core_id    = static_cast<int>(cfg.core_id.value_or(-99));
        const int    env_id     = static_cast<int>(cfg.envelope_id.value_or(-99));
        const int    vac_id     = static_cast<int>(cfg.vacuum_id.value_or(-99));

        int geom_order = -99;
        if (mesh->GetNodes()) {
            geom_order = mesh->GetNodes()->FESpace()->GetMaxElementOrder();
        }
        const int sorder = (sample_order > 0) ? sample_order : (2 * geom_order + 4);

        if (has_feature(features, MeshStatFeatures::CONFIG_META)) {
            ConfigMeta meta;
            meta.r_core             = cfg.r_core.value_or(-99.99);
            meta.r_star             = cfg.r_star.value_or(-99.99);
            meta.r_infinity         = cfg.r_infinity.value_or(-99.99);
            meta.flattening         = flattening;
            meta.geom_order         = geom_order;
            meta.refinement_levels  = sm.refinement_levels;
            meta.has_external_domain = cfg.include_external_domain.value_or(false);
            out.config_meta = meta;
        }

        if (has_feature(features, MeshStatFeatures::CONFORMITY)) {
            ConformityStats conformity;
            conformity.hierarchy_enabled = mesh->ncmesh != nullptr;
            for (int f = 0; f < mesh->GetNFaces(); ++f) {
                if (mesh->GetFaceInformation(f).IsNonconformingFine()) {
                    ++conformity.n_nonconforming_faces;
                }
            }
            conformity.conforming = conformity.n_nonconforming_faces == 0;
            out.conformity = conformity;
        }

        if (has_feature(features, MeshStatFeatures::REFINEMENT)) {
            RefinementStats refinement;
            auto update_depth = [](RegionRefinementStats& region, const int depth) {
                if (region.min_depth < 0) region.min_depth = depth;
                region.min_depth = std::min(region.min_depth, depth);
                region.max_depth = std::max(region.max_depth, depth);
            };
            for (int e = 0; e < mesh->GetNE(); ++e) {
                const int depth = mesh->ncmesh ? mesh->ncmesh->GetElementDepth(e) :
                    static_cast<int>(sm.refinement_levels);
                update_depth(refinement.all, depth);
                const int attr = mesh->GetAttribute(e);
                if      (attr == core_id) update_depth(refinement.core, depth);
                else if (attr == env_id)  update_depth(refinement.envelope, depth);
                else if (attr == vac_id)  update_depth(refinement.vacuum, depth);
            }
            if (const auto* fes = mesh->GetNodalFESpace()) {
                refinement.geometry_dofs = fes->GetNDofs();
                refinement.geometry_true_dofs = fes->GetNConformingDofs();
            }
            out.refinement = refinement;
        }

        // ============================ SURFACE PASS ============================
        const bool needs_surface =
            has_feature(features, MeshStatFeatures::RADIUS) ||
            has_feature(features, MeshStatFeatures::AXES) ||
            has_feature(features, MeshStatFeatures::ELLIPTICITY) ||
            has_feature(features, MeshStatFeatures::BOWING) ||
            has_feature(features, MeshStatFeatures::VOLUME_AREA);

        if (needs_surface) {
            double r_min = std::numeric_limits<double>::max();
            double r_max = std::numeric_limits<double>::lowest();
            double r_sum = 0.0, r_sum_sq = 0.0;
            double a_eq = 0.0, c_pol = 0.0;
            double bow_in = 0.0, bow_out = 0.0, bow_sum_sq = 0.0, bow_worst_r = 0.0;
            double bow_worst_mag = 0.0;
            double area = 0.0;
            long   n_samples = 0;
            mfem::Vector phys;

            for (int b = 0; b < mesh->GetNBE(); ++b) {
                if (mesh->GetBdrAttribute(b) != surf_bdr) continue;
                mfem::ElementTransformation* T = mesh->GetBdrElementTransformation(b);
                const mfem::IntegrationRule& ir = mfem::IntRules.Get(T->GetGeometryType(), sorder);

                for (int q = 0; q < ir.GetNPoints(); ++q) {
                    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
                    T->SetIntPoint(&ip);
                    T->Transform(ip, phys);

                    const double x = phys(0), y = phys(1), z = phys(2);
                    const double r = phys.Norml2();
                    const double rho = std::sqrt(x*x + y*y);
                    const double w = T->Weight() * ip.weight;

                    r_min = std::min(r_min, r);
                    r_max = std::max(r_max, r);
                    r_sum += r; r_sum_sq += r*r;
                    a_eq = std::max(a_eq, rho);
                    c_pol = std::max(c_pol, std::abs(z));
                    area += w;
                    ++n_samples;

                    if (has_feature(features, MeshStatFeatures::BOWING) && r > 1e-14) {
                        const double rt = SpheroidRadius(x/r, y/r, z/r, r_star, flattening);
                        const double dev = r - rt;
                        bow_in  = std::min(bow_in, dev);
                        bow_out = std::max(bow_out, dev);
                        bow_sum_sq += dev * dev;
                        if (std::abs(dev) > bow_worst_mag) {
                            bow_worst_mag = std::abs(dev);
                            bow_worst_r = r;
                        }
                    }
                }
            }

            if (n_samples == 0) {
                out.warnings.push_back("No samples were collected from the surface boundary. Check that the "
                                       "surface boundary ID is correct. Lacking a surface pass prevents the "
                                       "following from reporting accurate results: "
                                       "[RADIUS, AXES, ELLIPTICITY, BOWING, VOLUME_AREA]");
            } else {
                if (has_feature(features, MeshStatFeatures::RADIUS)) {
                    const double mean = r_sum / n_samples;
                    const double var  = std::max(0.0, r_sum_sq / n_samples - mean * mean);
                    out.radius = RadiusStats{
                        .min = r_min, .max = r_max, .mean = mean,
                        .stddev = std::sqrt(var), .n_samples = n_samples
                    };
                }
                if (has_feature(features, MeshStatFeatures::AXES)) {
                    out.axes = AxisStats{ .semi_major = a_eq, .semi_minor = c_pol };
                }
                if (has_feature(features, MeshStatFeatures::ELLIPTICITY)) {
                    out.ellipticity = EllipticityStats{
                        .flattening        = (a_eq > 0) ? (a_eq - c_pol) / a_eq : 0.0,
                        .polar_equatorial  = (a_eq > 0) ? c_pol / a_eq : 1.0,
                        .radius_uniformity = (r_max > 0) ? r_min / r_max : 1.0,
                    };
                }
                if (has_feature(features, MeshStatFeatures::BOWING)) {
                    out.bowing = BowingStats{
                        .max_inward = bow_in, .max_outward = bow_out,
                        .rms = std::sqrt(bow_sum_sq / n_samples), .worst_at_radius = bow_worst_r,
                    };
                }
                if (has_feature(features, MeshStatFeatures::VOLUME_AREA)) {
                    const double a = r_star, c = r_star * (1.0 - flattening);
                    out.volume = VolumeAreaStats{
                        .surface_area = area,
                        .analytic_area = (flattening == 0) ? 4.0 * M_PI * a * a : -99.99
                    };
                }
            }
        }

        // ============================ VOLUME PASS ============================
        const bool needs_volume =
            has_feature(features, MeshStatFeatures::VOLUME_AREA) ||
            has_feature(features, MeshStatFeatures::JACOBIAN) ||
            has_feature(features, MeshStatFeatures::ELEMENT_COUNT) ||
            has_feature(features, MeshStatFeatures::MESH_SIZE) ||
            has_feature(features, MeshStatFeatures::CENTROID) ||
            has_feature(features, MeshStatFeatures::BOUNDING_BOX);   // BUGFIX: was a dangling ';'

        if (needs_volume) {
            const bool need_bbox = has_feature(features, MeshStatFeatures::BOUNDING_BOX);
            const bool need_jac  = has_feature(features, MeshStatFeatures::JACOBIAN);

            // Bounding-box accumulators (seeded inverted so empty regions stay invalid).
            double cxmin=+std::numeric_limits<double>::max(), cxmax=-std::numeric_limits<double>::max();
            double cymin=cxmin, cymax=cxmax, czmin=cxmin, czmax=cxmax;                            // core
            double sxmin=cxmin, sxmax=cxmax, symin=cxmin, symax=cxmax, szmin=cxmin, szmax=cxmax;  // stellar
            double vxmin=cxmin, vxmax=cxmax, vymin=cxmin, vymax=cxmax, vzmin=cxmin, vzmax=cxmax;  // vacuum
            long n_core_box=0, n_stel_box=0, n_vac_box=0;
            mfem::Vector bphys;

            // Per-region Jacobian accumulators.
            struct JacAccum {
                double detJ_min = std::numeric_limits<double>::max();
                double detJ_max = -std::numeric_limits<double>::max();
                double min_ratio = 1.0;
                long   n_flipped = 0;
                long   n_elem = 0;
                double worst_ratio_r = -1.0;
                double min_detJ_r = -1.0;
            };
            JacAccum all_acc, stel_acc, vac_acc;

            auto jac_update = [](JacAccum& a, double dmin, double dmax, bool flip, double r) {
                ++a.n_elem;
                if (dmin < a.detJ_min) { a.detJ_min = dmin; a.min_detJ_r = r; }
                if (dmax > a.detJ_max)   a.detJ_max = dmax;
                if (dmax > 1e-30) {
                    const double ratio = dmin / dmax;
                    if (ratio < a.min_ratio) { a.min_ratio = ratio; a.worst_ratio_r = r; }
                }
                if (flip) ++a.n_flipped;
            };

            double vol = 0.0, cx = 0.0, cy = 0.0, cz = 0.0;
            double h_min = std::numeric_limits<double>::max();
            double h_max = -std::numeric_limits<double>::max();
            double h_sum = 0.0, h_sum_sq = 0.0;
            long   n_core = 0, n_env = 0, n_vac = 0, n_other = 0;
            mfem::Vector phys;

            for (int e = 0; e < mesh->GetNE(); ++e) {
                const int attr = mesh->GetAttribute(e);
                if      (attr == core_id) ++n_core;
                else if (attr == env_id)  ++n_env;
                else if (attr == vac_id)  ++n_vac;
                else                      ++n_other;

                const bool stellar = (attr == core_id || attr == env_id);

                if (has_feature(features, MeshStatFeatures::MESH_SIZE) && stellar) {
                    const double h = mesh->GetElementSize(e);
                    h_min = std::min(h_min, h);
                    h_max = std::max(h_max, h);
                    h_sum += h; h_sum_sq += h*h;
                }

                mfem::ElementTransformation* T = mesh->GetElementTransformation(e);
                const mfem::IntegrationRule& ir = mfem::IntRules.Get(T->GetGeometryType(), sorder);

                double e_detmin = std::numeric_limits<double>::max();
                double e_detmax = -std::numeric_limits<double>::max();
                bool   e_flip = false;

                for (int q = 0; q < ir.GetNPoints(); ++q) {
                    const mfem::IntegrationPoint& ip = ir.IntPoint(q);
                    T->SetIntPoint(&ip);

                    if (need_bbox) {
                        T->Transform(ip, bphys);
                        const double X = bphys(0), Y = bphys(1), Z = bphys(2);
                        if (attr == core_id) {
                            cxmin=std::min(cxmin,X); cxmax=std::max(cxmax,X);
                            cymin=std::min(cymin,Y); cymax=std::max(cymax,Y);
                            czmin=std::min(czmin,Z); czmax=std::max(czmax,Z);
                            ++n_core_box;
                        }
                        if (stellar) {   // stellar = core U envelope
                            sxmin=std::min(sxmin,X); sxmax=std::max(sxmax,X);
                            symin=std::min(symin,Y); symax=std::max(symax,Y);
                            szmin=std::min(szmin,Z); szmax=std::max(szmax,Z);
                            ++n_stel_box;
                        }
                        if (attr == vac_id) {
                            vxmin=std::min(vxmin,X); vxmax=std::max(vxmax,X);
                            vymin=std::min(vymin,Y); vymax=std::max(vymax,Y);
                            vzmin=std::min(vzmin,Z); vzmax=std::max(vzmax,Z);
                            ++n_vac_box;
                        }
                    }

                    const double dJ = T->Jacobian().Det();
                    e_detmin = std::min(e_detmin, dJ);
                    e_detmax = std::max(e_detmax, dJ);
                    if (dJ < 0.0) e_flip = true;

                    if (stellar && (has_feature(features, MeshStatFeatures::VOLUME_AREA) ||
                                    has_feature(features, MeshStatFeatures::CENTROID))) {
                        const double w = std::abs(dJ) * ip.weight;
                        vol += w;
                        if (has_feature(features, MeshStatFeatures::CENTROID)) {
                            T->Transform(ip, phys);
                            cx += w*phys(0); cy += w*phys(1); cz += w*phys(2);
                        }
                    }
                }

                if (need_jac) {
                    // Representative element radius (center) for locating the worst element.
                    const mfem::IntegrationPoint& cip = mfem::Geometries.GetCenter(T->GetGeometryType());
                    T->SetIntPoint(&cip);
                    mfem::Vector cpt;
                    T->Transform(cip, cpt);
                    const double er = cpt.Norml2();

                    jac_update(all_acc, e_detmin, e_detmax, e_flip, er);
                    if (stellar)             jac_update(stel_acc, e_detmin, e_detmax, e_flip, er);
                    else if (attr == vac_id) jac_update(vac_acc,  e_detmin, e_detmax, e_flip, er);
                }
            }

            if (has_feature(features, MeshStatFeatures::ELEMENT_COUNT)) {
                out.element_counts = ElementCounts{
                    .total = n_core + n_env + n_vac + n_other,
                    .core = n_core, .envelope = n_env, .vacuum = n_vac, .other = n_other,
                    .n_vertices = mesh->GetNV(),
                };
            }

            if (need_jac) {
                auto finalize = [](const JacAccum& a) {
                    JacobianStats j;
                    j.n_elements           = a.n_elem;
                    j.detJ_min             = (a.n_elem > 0) ? a.detJ_min : 0.0;
                    j.detJ_max             = (a.n_elem > 0) ? a.detJ_max : 0.0;
                    j.min_detJ_ratio       = a.min_ratio;
                    j.n_flipped            = a.n_flipped;
                    j.worst_ratio_at_radius = a.worst_ratio_r;
                    j.detJ_min_at_radius    = a.min_detJ_r;
                    return j;
                };
                out.jacobian = finalize(all_acc);
                if (stel_acc.n_elem > 0) out.jacobian_stellar = finalize(stel_acc);
                if (vac_acc.n_elem  > 0) out.jacobian_vacuum  = finalize(vac_acc);
            }

            if (has_feature(features, MeshStatFeatures::MESH_SIZE)) {
                const long ns = n_core + n_env;
                const double mean = (ns > 0) ? h_sum / ns : 0.0;
                const double var  = (ns > 0) ? std::max(0.0, h_sum_sq / ns - mean * mean) : 0.0;
                out.mesh_size = MeshSizeStats{
                    .h_min = h_min, .h_max = h_max, .h_mean = mean, .h_stddev = std::sqrt(var),
                };
            }

            if (has_feature(features, MeshStatFeatures::VOLUME_AREA)) {
                if (!out.volume) out.volume.emplace();   // BUGFIX: surface pass may not have created it
                out.volume->stellar_volume = vol;
                const double a = r_star, c = r_star * (1.0 - flattening);
                out.volume->analytic_volume = (4.0 / 3.0) * M_PI * a * a * c;
            }

            if (has_feature(features, MeshStatFeatures::CENTROID)) {
                if (vol <= 0) {
                    out.warnings.push_back("Stellar volume is zero or negative, cannot compute centroid.");
                } else {
                    const double ccx = cx / vol, ccy = cy / vol, ccz = cz / vol;  // BUGFIX: normalize
                    out.centroid = CentroidStats{
                        .x = ccx, .y = ccy, .z = ccz,
                        .offset = std::sqrt(ccx*ccx + ccy*ccy + ccz*ccz)
                    };
                }
            }

            if (need_bbox) {
                BoundingBoxStats bb;
                auto fill = [](BoundingBox& box, long n,
                               double xmn,double xmx,double ymn,double ymx,double zmn,double zmx) {
                    if (n > 0) {
                        box.valid = true;
                        box.xMin=xmn; box.xMax=xmx;
                        box.yMin=ymn; box.yMax=ymx;
                        box.zMin=zmn; box.zMax=zmx;
                    }
                };
                fill(bb.core,   n_core_box, cxmin,cxmax,cymin,cymax,czmin,czmax);
                fill(bb.star,   n_stel_box, sxmin,sxmax,symin,symax,szmin,szmax);
                fill(bb.vacuum, n_vac_box,  vxmin,vxmax,vymin,vymax,vzmin,vzmax);
                out.bounding_box = bb;
            }
        }

        // ============================ OUTER BOUND PASS ============================
        if (has_feature(features, MeshStatFeatures::OUTER_BOUNDS)) {
            double r_min = std::numeric_limits<double>::max();
            double r_max = std::numeric_limits<double>::lowest();
            double r_sum = 0.0;
            long   n_samples = 0;
            mfem::Vector phys;

            for (int b = 0; b < mesh->GetNBE(); ++b) {
                if (mesh->GetBdrAttribute(b) != inf_bdr) continue;
                mfem::ElementTransformation* T = mesh->GetBdrElementTransformation(b);
                const mfem::IntegrationRule& ir = mfem::IntRules.Get(T->GetGeometryType(), sorder);
                for (int q = 0; q < ir.GetNPoints(); ++q) {
                    T->SetIntPoint(&ir.IntPoint(q));
                    T->Transform(ir.IntPoint(q), phys);
                    const double r = phys.Norml2();
                    r_min = std::min(r_min, r); r_max = std::max(r_max, r);
                    r_sum += r; ++n_samples;
                }
            }

            if (n_samples == 0) {
                out.warnings.push_back("No samples found on the outer boundary, cannot compute outer bounds.");
            } else {
                out.outer_bounds = OuterBoundsStats{
                    .min = r_min, .max = r_max, .mean = r_sum / n_samples, .n_samples = n_samples
                };
            }
        }

        return out;
    }

    std::string to_string(const MeshStats& s) {
        std::string o = "MeshStats:\n";
        auto line = [&](const std::string& l){ o += " =>" + l + "\n"; };

        if (s.config_meta) {
            const auto& m = *s.config_meta;
            line(std::format(
                "config: r_core={:0.4f}, r_star={:0.4f}, r_inf={:0.4f}, flattening={:0.4f}, "
                "geometric order={}, refinement levels={}",
                m.r_core, m.r_star, m.r_infinity, m.flattening, m.geom_order, m.refinement_levels));
        }
        if (s.radius) {
            const auto& r = *s.radius;
            line(std::format("radius: min={:.6f} max={:.6f} mean={:.6f} std={:.3E} (n={})",
                 r.min, r.max, r.mean, r.stddev, r.n_samples));
        }
        if (s.axes) {
            line(std::format("axes: semi_major(eq)={:.6f} semi_minor(pol)={:.6f}",
                 s.axes->semi_major, s.axes->semi_minor));
        }
        if (s.ellipticity) {
            const auto& e = *s.ellipticity;
            line(std::format("ellipticity: flattening={:.5f} c/a={:.5f} r_min/r_max={:.5f}",
                 e.flattening, e.polar_equatorial, e.radius_uniformity));
        }
        if (s.bowing) {
            const auto& b = *s.bowing;
            line(std::format("bowing: max_inward={:.3E} max_outward={:.3E} rms={:.3E}",
                 b.max_inward, b.max_outward, b.rms));
        }
        if (s.conformity) {
            const auto& c = *s.conformity;
            line(std::format("conforming: {} (hierarchy={}, hanging face patches={})",
                 c.conforming, c.hierarchy_enabled, c.n_nonconforming_faces));
        }
        if (s.refinement) {
            const auto& r = *s.refinement;
            line(std::format(
                "refinement depth: all=[{},{}] core=[{},{}] env=[{},{}] vac=[{},{}]",
                r.all.min_depth, r.all.max_depth, r.core.min_depth, r.core.max_depth,
                r.envelope.min_depth, r.envelope.max_depth, r.vacuum.min_depth, r.vacuum.max_depth));
            line(std::format("scalar geometry DOFs: total={} true={}",
                 r.geometry_dofs, r.geometry_true_dofs));
        }

        auto jac_line = [&](const std::string& label, const JacobianStats& j) {
            line(std::format(
                "jacobian[{}]: detJ=[{:.3E},{:.3E}] min_ratio={:.3E} (@r={:.4f}) "
                "min_detJ@r={:.4f} flipped={} n={}",
                label, j.detJ_min, j.detJ_max, j.min_detJ_ratio, j.worst_ratio_at_radius,
                j.detJ_min_at_radius, j.n_flipped, j.n_elements));
        };
        if (s.jacobian)         jac_line("all",     *s.jacobian);
        if (s.jacobian_stellar) jac_line("stellar", *s.jacobian_stellar);
        if (s.jacobian_vacuum)  jac_line("vacuum",  *s.jacobian_vacuum);

        if (s.volume) {
            const auto& v = *s.volume;
            line(std::format("volume={:.6f} (analytic {:.6f})  area={:.6f}",
                 v.stellar_volume, v.analytic_volume, v.surface_area));
        }
        if (s.element_counts) {
            const auto& c = *s.element_counts;
            line(std::format("elements: total={} core={} env={} vac={} other={} NV={}",
                 c.total, c.core, c.envelope, c.vacuum, c.other, c.n_vertices));
        }
        if (s.mesh_size) {
            line(std::format("h: min={:.4E} max={:.4E} mean={:.4E} std={:.4E}",
                 s.mesh_size->h_min, s.mesh_size->h_max, s.mesh_size->h_mean, s.mesh_size->h_stddev));
        }
        if (s.bounding_box) {
            const auto& bb = *s.bounding_box;
            auto bline = [&](const char* nm, const BoundingBox& b){
                if (b.valid)
                    line(std::format("bbox[{}]: x[{:.4f},{:.4f}] y[{:.4f},{:.4f}] z[{:.4f},{:.4f}]",
                         nm, b.xMin,b.xMax, b.yMin,b.yMax, b.zMin,b.zMax));
                else
                    line(std::format("bbox[{}]: <absent>", nm));
            };
            bline("core", bb.core);
            bline("star", bb.star);
            bline("vacuum", bb.vacuum);
        }
        if (s.outer_bounds) {
            line(std::format("outer: min={:.4f} max={:.4f} mean={:.4f}",
                 s.outer_bounds->min, s.outer_bounds->max, s.outer_bounds->mean));
        }
        if (s.centroid) {
            line(std::format("centroid: x={:.6f} y={:.6f} z={:.6f} offset={:.6f}",
                 s.centroid->x, s.centroid->y, s.centroid->z, s.centroid->offset));
        }

        for (const auto& w : s.warnings) line("WARNING: " + w);
        for (const auto& e : s.errors)   line(std::format("ERROR: {}", e));
        return o;
    }
}
