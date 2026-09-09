#include "mfem.hpp"

#include <thread>
#include <atomic>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "stroid/topology/optimize.h"


#include <clocale>
#include <cstdlib>
#include <cstring>
#include <string>

namespace stroid::utils::term_support {
    inline bool locale_name_looks_utf8(const char* localeName) {
        if (localeName == nullptr) return false;

        const std::string localeString(localeName);

        return localeString.find("UTF-8") != std::string::npos ||
               localeString.find("utf-8") != std::string::npos ||
               localeString.find("utf8")  != std::string::npos ||
               localeString.find("UTF8")  != std::string::npos;
    }

    inline bool unicode_output_is_usable() {
        const char* ctypeLocale = std::setlocale(LC_CTYPE, "");
        if (locale_name_looks_utf8(ctypeLocale)) {
            return true;
        }

        const char* lcAllEnv = std::getenv("LC_ALL");
        if (locale_name_looks_utf8(lcAllEnv)) {
            return true;
        }

        const char* lcCtypeEnv = std::getenv("LC_CTYPE");
        if (locale_name_looks_utf8(lcCtypeEnv)) {
            return true;
        }

        const char* langEnv = std::getenv("LANG");
        if (locale_name_looks_utf8(langEnv)) {
            return true;
        }

        return false;
    }
}

namespace stroid::topology {
class TMOPProgressBar : public mfem::IterativeSolverMonitor {
    private:
        double r0_ = -1.0;
        double rtol_;
        int bar_width_;

        std::atomic<bool> done_{false};
        std::atomic<double> progress_{0.0};
        std::atomic<int> iter_{0};
        std::atomic<double> res_{0.0};

        std::thread spinner_thread_;

        void Spin() {
            std::vector<std::string> spin_chars;
            if (!utils::term_support::unicode_output_is_usable()) {
                spin_chars= {"|", "/", "-", "\\"};
            } else {
                spin_chars = {"▉", "▊", "▋", "▌", "▍", "▎", "▏", "▎", "▍", "▌", "▋", "▊", "▉"};
            }
            int spin_idx = 0;

            while (!done_.load()) {
                Draw(spin_chars[spin_idx]);
                spin_idx = (spin_idx + 1) % spin_chars.size();
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        }

        void Draw(const std::string& spinner) {
            const double p = progress_.load();
            const int pos = static_cast<int>(bar_width_ * p);

            std::cout << "\r[" << spinner << "] TMOP Relaxation [";
            for (int i = 0; i < bar_width_; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << std::setw(3) << static_cast<int>(p * 100.0) << "% "
                      << "(Iter: " << std::setw(2) << iter_.load()
                      << ", Res: " << std::scientific << std::setprecision(2) << res_.load() << ") " << std::flush;
        }

    public:
        TMOPProgressBar(double rel_tol, int width = 50)
            : rtol_(rel_tol), bar_width_(width) {
            spinner_thread_ = std::thread(&TMOPProgressBar::Spin, this);
        }

        ~TMOPProgressBar() override {
            if (!done_.load()) {
                done_ = true;
                if (spinner_thread_.joinable()) {
                    spinner_thread_.join();
                }
            }
        }

        void MonitorResidual(int it, double norm, const mfem::Vector &r, bool final) override {
            if (it == 0 || r0_ < 0.0) {
                r0_ = norm;
            }

            iter_ = it;
            res_ = norm;

            double p = 0.0;
            const double target_norm = r0_ * rtol_;

            if (norm <= target_norm || final) {
                p = 1.0;
            } else if (norm < r0_ && r0_ > 0.0 && target_norm > 0.0) {
                const double log_start = std::log10(r0_);
                const double log_current = std::log10(norm);
                const double log_target = std::log10(target_norm);
                p = (log_start - log_current) / (log_start - log_target);
                p = std::clamp(p, 0.0, 1.0);
            }

            progress_ = p;

            if (final) {
                done_ = true;
                if (spinner_thread_.joinable()) {
                    spinner_thread_.join();
                }

                Draw("*");
                std::cout << std::endl;
            }
        }
    };
    void ApplyTMOP(mfem::Mesh &mesh, const fourdst::config::Config<config::MeshConfig> &config) {
        const mfem::FiniteElementSpace* cfes = mesh.GetNodalFESpace();
        mfem::FiniteElementSpace* fes = const_cast<mfem::FiniteElementSpace*>(cfes);

        if (!fes) {
            std::cerr << "Error: Mesh has no nodal finite element space. Call PromoteToHighOrder first." << std::endl;
            return;
        }

        const int max_bdr_attr = mesh.bdr_attributes.Size() > 0 ? mesh.bdr_attributes.Max() : 0;
        mfem::Array<int> ess_bdr(max_bdr_attr);
        ess_bdr = 0.0;

        if (max_bdr_attr >= config->surface_bdr_id.value()) {
            ess_bdr[config->surface_bdr_id.value() - 1] = 1;
        }

        if (config->include_external_domain.value() && max_bdr_attr >= config->inf_bdr_id.value()) {
            ess_bdr[config->inf_bdr_id.value() - 1] = 1;
        }

        mfem::Array<int> ess_tdof_list;
        fes->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

        mfem::TMOP_QualityMetric* metric = new mfem::TMOP_Metric_302();
        mfem::TargetConstructor* target_c = new mfem::TargetConstructor(mfem::TargetConstructor::IDEAL_SHAPE_GIVEN_SIZE);
        target_c->SetNodes(*mesh.GetNodes());
        mfem::TMOP_Integrator* tmop_integrator = new mfem::TMOP_Integrator(metric, target_c);

        mfem::NonlinearForm a(fes);
        a.AddDomainIntegrator(tmop_integrator);
        a.SetEssentialTrueDofs(ess_tdof_list);

        mfem::GridFunction* nodes = mesh.GetNodes();
        mfem::Vector x;
        nodes->GetTrueDofs(x);
        mfem::Vector b(a.Height());
        b = 0.0;

        mfem::MINRESSolver minres;
        minres.SetMaxIter(750);
        minres.SetRelTol(1e-5);
        minres.SetAbsTol(0.0);
        minres.SetPrintLevel(-1);

        mfem::DSmoother jacobi(1, 1.0, 1);
        jacobi.SetPositiveDiagonal(true);
        minres.SetPreconditioner(jacobi);

        const int quad_order = 2 * fes->GetMaxElementOrder() + 3;
        const mfem::IntegrationRule &ir = mfem::IntRules.Get(mesh.GetTypicalElementGeometry(), quad_order);

        double min_detJ = std::numeric_limits<double>::infinity();
        for (int i = 0; i < mesh.GetNE(); i++) {
            mfem::ElementTransformation *T = mesh.GetElementTransformation(i);
            for (int j = 0; j < ir.GetNPoints(); j++) {
                T->SetIntPoint(&ir.IntPoint(j));
                min_detJ = std::min(min_detJ, T->Jacobian().Det());
            }
        }

        constexpr double newton_rtol = 1e-8;
        mfem::TMOPNewtonSolver newton(ir, 0);
        newton.SetPreconditioner(minres);
        newton.SetOperator(a);
        newton.SetMaxIter(50);
        newton.SetRelTol(newton_rtol);
        newton.SetAbsTol(0.0);
        newton.SetMinDetPtr(&min_detJ);
        newton.SetPrintLevel(0);

        TMOPProgressBar progress_bar(newton_rtol);
        newton.SetMonitor(progress_bar);

        std::cout << "Applying TMOP optimization to mesh. Note this may take a long time. Depending on your mesh resolution expect to wait up to the order of 10s of minutes..." << std::endl;
        newton.Mult(b, x);
        nodes->SetFromTrueDofs(x);

        mesh.NodesUpdated();

        delete metric;
        delete target_c;
    }

    void OptimizeMesh(mfem::Mesh& mesh, const fourdst::config::Config<config::MeshConfig> &cfg) {
        if (cfg->optimization_methods.has_value() && cfg->optimization_methods.value().tmop.has_value() && cfg->optimization_methods.value().tmop.value()) {
            ApplyTMOP(mesh, cfg);
        }
    }
}
