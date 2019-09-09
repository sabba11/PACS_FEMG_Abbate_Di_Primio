//GetFEM++ libraries
#include <gmm/gmm.h>
#include <gmm/gmm_except.h>

//Standard libraries
#include <iostream>
#include <fstream>

//Project libraries
#include "quantum_graph_problem.hpp"

namespace getfem{

// INITIALIZE THE PROBLEM

void
quantum_graph_problem::init(int argc, char *argv[])
{
     //1. Read the .param filename from standard input
     INPUT.read_command_line(argc, argv);

     //2. Import data (algorithm specifications, boundary conditions,...)
     import_data();

     //3. Build mesh for the graph
     build_mesh();

     //4. Set finite elements and integration methods
     set_im_and_fem();

     //5. Build problem parameters
     //build_param();

     //6. Build the lists of the data of the vertices
     build_vertices_lists();
}

void
quantum_graph_problem::import_data()
{
    #ifdef FEMG_VERBOSE_
    std::cout << "Importing descriptors for the problem..."<< std::endl;
    #endif

    descr.import(INPUT);

    #ifdef FEMG_VERBOSE_
    std::cout << descr;
    #endif
}

void
quantum_graph_problem::build_mesh()
{
    #ifdef FEMG_VERBOSE_
    std::cout << "Importing the mesh for the graph..."<< std::endl;
    #endif

    std::ifstream ifs(descr.MESH_FILEG);
    GMM_ASSERT1(ifs.good(),"Unable to read from file " << descr.MESH_FILEG);

    import_pts_file(ifs, meshg, BCg, n_vertices, descr.MESH_TYPEG);

    n_branches = n_vertices.size();

    ifs.close();
}

void
quantum_graph_problem::set_im_and_fem()
{
    #ifdef FEMG_VERBOSE_
    std::cout << "Setting Integration Methods for the discrete problem..."<< std::endl;
    #endif

    pintegration_method pim_g = int_method_descriptor(descr.IM_TYPEG);
    mimg.set_integration_method(meshg.convex_index(), pim_g);

    #ifdef FEMG_VERBOSE_
    std::cout << "Setting Finite Element Methods for the discrete problem..."<< std::endl;
    #endif

    bgeot::pgeometric_trans pgt_g = bgeot::geometric_trans_descriptor(descr.MESH_TYPEG);

    pfem pf_Ug = fem_descriptor(descr.FEM_TYPEG);
    pfem pf_coeffg = fem_descriptor(descr.FEM_TYPEG_DATA);

    // mf_coeffbranchg.reserve(n_branches);
    // for (size_type i = 0; i < n_branches; ++i) {
    //     mesh_fem mf_tmp(meshg);
    //     mf_tmp.set_finite_element(meshg.region(i).index(), pf_coeffg);
    //     mf_coeffbranchg.emplace_back(mf_tmp);
    //     mf_tmp.clear();
    // }

    mf_Ug.set_finite_element(meshg.convex_index(), pf_Ug);
    mf_coeffg.set_finite_element(meshg.convex_index(), pf_coeffg);

    // #ifdef FEMG_VERBOSE_
    // std::cout << "Setting FEM dimensions for the problem..." << std::endl;
    // #endif
    // dof.set(mf_Ug, mf_coeffg);
    // #ifdef FEMG_VERBOSE_
    // std::cout << std::scientific << dof;
    // #endif
}

void
quantum_graph_problem::build_param(void)
{
    #ifdef FEMG_VERBOSE_
    std::cout << "Building parameters for the problem..." << std::endl;
    #endif
    //param.build(INPUT, mf_coeffg, mf_coeffbranchg);
    #ifdef FEMG_VERBOSE_
    //std::cout << param;
    #endif
}

void
quantum_graph_problem::build_vertices_lists(void)
{
    //to be defined
}

}// end of namespace
