#include "mfem.hpp"
#include "stroid/config/config.h"
#include "stroid/IO/mesh.h"

#include <fstream>
#include <iostream>
#include <cstdint>

namespace stroid::IO {


    void SaveMesh(const mfem::Mesh& mesh, const std::string& filename) {
        std::ofstream ofs(filename);
        ofs.precision(8);
        mesh.Print(ofs);
    }

    void SaveVTU(mfem::Mesh &mesh, const std::string &exportName) {
        mfem::ParaViewDataCollection pd(exportName, &mesh);
        pd.SetDataFormat(mfem::VTKFormat::BINARY);
        pd.SetHighOrderOutput(true);
        pd.Save();
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
    void VisualizeFaceValence(mfem::Mesh& mesh) {
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
        char vishost[] = "localhost";
        int  visport   = 19916;
        mfem::socketstream sol_sock(vishost, visport);
        if (sol_sock.is_open()) {
            sol_sock << "solution\n" << mesh << valence_gf;
            sol_sock << "window_title 'Boundary Valence: 1=Surface, 2=Internal'\n";
            sol_sock << "keys am\n" << std::flush;
        }
    }
}