// Source file containing definitions for class eigen_problem.

// GetFEM++ libraries
#include <getfem/getfem_assembling.h>
#include <gmm/gmm_condition_number.h>
#include <gmm/gmm_except.h>
#include <gmm/gmm_iter_solvers.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm.h>

// Standard libraries
#include <algorithm>
#include <map>
#include <complex>
#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/filesystem.hpp>

// Project headers
#include "eigen_problem.hpp"

namespace getfem {

	void
	eigen_problem::init(int argc, char *argv[])
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
	     // build_vertices_lists();
	}

	void
	eigen_problem::import_data()
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "Importing descriptors for the problem..."<< std::endl;
	    #endif

	    descr.import_all(INPUT);

	    #ifdef FEMG_VERBOSE_
	    std::cout << descr;
	    #endif
	}

	void
	eigen_problem::build_mesh()
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "Importing the mesh for the graph..."<< std::endl;
	    #endif

	    std::ifstream ifs(descr.MESH_FILEG);
	    GMM_ASSERT1(ifs.good(),"Unable to read from file " << descr.MESH_FILEG);
      //ACHTUNG!: for the moment we don't give mesh step as an input parameter
			// if we reimplement it back you should give it
	    import_pts_file(ifs, meshg, BCg, n_origvert, n_vertices, descr.MESH_TYPEG);

	    n_branches = n_vertices.size();
  		n_totalvert = meshg.nb_points();

	    ifs.close();
	}

	void
	eigen_problem::set_im_and_fem()
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

	    mf_Ug.set_finite_element(meshg.convex_index(), pf_Ug);
	    mf_coeffg.set_finite_element(meshg.convex_index(), pf_coeffg);
	}

	void
	eigen_problem::build_param(void)
	{
	    #ifdef FEMG_VERBOSE_
	    std::cout << "Building parameters for the problem..." << std::endl;
	    #endif
	    //param.build(INPUT, mf_coeffg, mf_coeffbranchg);
	    #ifdef FEMG_VERBOSE_
	    //std::cout << param;
	    #endif
	}


  // Possibly, build this function
	// void
	// eigen_problem::build_vertices_lists(void)
	// {
	//     //to be defined
	// }

	void
	eigen_problem::assembly(void)
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling matrices..." << std::endl;
		#endif

		assembly_matL();
		assembly_matM();
		assembly_matA();

		return;
	}

	void
	eigen_problem::assembly_matL()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling stiffness matrix L..." << std::endl;
		#endif

		L.resize(n_totalvert, n_totalvert);
		vector_type coeff_val(n_totalvert, 1);
		getfem::asm_stiffness_matrix_for_laplacian(L, mimg, mf_Ug, mf_coeffg, coeff_val);
		gmm::clean(L, 1E-10);

		return;
	}

	void
	eigen_problem::assembly_matM()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling mass matrix M..." << std::endl;
		#endif

		M.resize(n_totalvert, n_totalvert);
		getfem::asm_mass_matrix(M, mimg, mf_Ug, mf_coeffg);
		gmm::clean(M, 1E-10);

		return;
	}

	void
	eigen_problem::assembly_matA()
	{
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Assembling operator matrix A..." << std::endl;
		#endif

		A.resize(n_totalvert, n_totalvert);
		if (descr.OPERATOR == "Laplacian")
			gmm::copy(L, A);

		else if (descr.OPERATOR == "Hamiltonian")
			gmm::add(L, M, A);

		else
			std::cout << "[eigen_problem] Invalid operator descriptor." << std::endl;

		return;
	}

	bool
	eigen_problem::solve(void)
	{
		if (descr.COMP_METHOD == "QR") {
			#ifdef FEMG_VERBOSE_
			std::cout << "[eigen_problem] Solving the problem via QR algorithm..." << std::endl;
			#endif
			try {
				clock_t time_inverse, time_eig, time_mult;
				dense_matrix_type eigvects(n_totalvert, n_totalvert);
				// Initial guess for all eigenvalues is 1.
				complex_vector_type eigvals(n_totalvert, 1);
				dense_matrix_type inverse_mass(n_totalvert, n_totalvert);
				dense_matrix_type aux_H(n_totalvert, n_totalvert);
				gmm::copy(M, inverse_mass);
				gmm::copy(A, aux_H);
				scalar_type cond_number = gmm::condition_number(M);
				std::cout << "[eigen_problem] The mass matrix M has condition number " << cond_number << std::endl;
				time_inverse = clock();
				gmm::lu_inverse(inverse_mass);
				time_inverse = clock() - time_inverse;
				log_data.push_back(std::make_pair("Time for matrix inversion (seconds): ", static_cast<float>(time_inverse)/CLOCKS_PER_SEC));
				gmm::mult(inverse_mass, aux_H, A);
				std::cout << static_cast<float>(time_inverse)/CLOCKS_PER_SEC << std::endl;
				scalar_type tol = descr.TOL;
				std::cout << "[eigen_problem] Starting QR routine..." << std::endl;
				dense_matrix_type dense_A(n_totalvert, n_totalvert);
				gmm::copy(A, dense_A);
				time_eig = clock();
				gmm::geev_interface_right(dense_A, eigvals, eigvects);
				time_eig = clock() - time_eig;
				log_data.push_back(std::make_pair("Time to compute eigencouples (seconds): ", static_cast<float>(time_eig)/CLOCKS_PER_SEC));
				#ifdef FEMG_VERBOSE_
				for (unsigned i = 0; i < log_data.size(); i++)
					std::cout << "[eigen_problem] " << log_data[i].first << log_data[i].second << std::endl;
				#endif
				for (unsigned i = 0; i < n_totalvert; i++) {
					vector_type aux_v;
					for (unsigned j = 0; j < n_totalvert; j++)
						aux_v.push_back(eigvects(j, i));
					if (abs(eigvals[i].imag()) < 1e-10) {
						auto eigpair = std::make_pair(eigvals[i].real(), aux_v);
						eigpairs.insert(eigpair);
					}
					else
						std::cout << "[eigen_problem] Warning: complex eigenvalue ignored. Value: " << eigvals[i].real() << " + " << eigvals[i].imag() << "i" << std::endl;
				}
			}
			GMM_STANDARD_CATCH_ERROR;
		}

		else if (descr.COMP_METHOD == "pencil_thm") {
			// to be defined
		}

		else {
			std::cout << "[eigen_problem] Invalid computational method descriptor." << std::endl;
		}

		return true;
	}

	void
	eigen_problem::sol_export(void)
	{
	// Temporary code.
		std::cout << "[eigen_problem] Printing eigenvalues..." << std::endl;
		for(auto it = eigpairs.begin(); it != eigpairs.end(); it++)
		   std::cout << it->first << std::endl;
		// Number of vertices in the graph without the mesh
		std::cout << "Number of original vertices: " << n_origvert << std::endl;
		// Number of vertices in the extended graph
		std::cout << "Number of total vertices: " << n_totalvert << std::endl;

		// Exporting routine to MATLAB interface.
		// Executing this code requires GetFEM++-MATLAB interface and the boost::filesystem library.
		// sudo nano ~/.bashrc
		// ultima riga export LD_LIBRARY_PATH=path/to/boost/lib:$LD_LIBRARY_PATH
		std::string output = descr.OUTPUT;
		unsigned n_digits = std::to_string(n_totalvert).size();
		unsigned i = 1;
		std::ostringstream dir_name_builder;
		dir_name_builder << output << "export/" << descr.OPERATOR << "/"<< std::to_string(n_totalvert) + " point-mesh/" << descr.COMP_METHOD;
		boost::filesystem::path dir(dir_name_builder.str().c_str());
		boost::system::error_code ec;
		bool status = boost::filesystem::create_directories(dir, ec);
		if( status )
			std::cout << "[eigen_problem] Directory " << dir_name_builder.str() << " created." << std::endl;
		else {
			if ( ec.message() == "Success" )
				std::cout << "[eigen_problem] Directory " << dir_name_builder.str() << " already exists." << std::endl;
			else {
				std::cerr << "[eigen_problem] Error while creating export directory: " << ec.message() << std::endl;
				return;
			}
		}
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Exporting eigenvalues and eigenvector files..." << std::endl;
		#endif
		std::ostringstream eigval_name_builder, log_name_builder;
		std::string s(n_digits, '0');
		eigval_name_builder << dir_name_builder.str() << "/" << s << "-eigenvalues.U";
		log_name_builder << dir_name_builder.str() << "/" << "export.log";
		std::fstream log(log_name_builder.str(), std::ios::out);
		if (!log) {
			std::cerr << "[eigen_problem] Unable to open log file " << log_name_builder.str() << std::endl;
		}
		std::fstream eig(eigval_name_builder.str(), std::ios::out);
		if (!eig) {
			log << "Exporting failed.\n";
			log << "Unable to open file " << eigval_name_builder.str() << " to export eigenvalues.\n";
			return;
		}
		for (auto it = eigpairs.begin(); it != eigpairs.end(); it++) {
			unsigned n_zeros = n_digits - std::to_string(i).size();
			std::string zeros(n_zeros, '0');
			scalar_type eigen = it->first;
			std::ostringstream file_name_builder;
	 		file_name_builder << dir_name_builder.str() << "/" << zeros << i << "-eigenvector-" << eigen << ".U";
	 		std::fstream f(file_name_builder.str(), std::ios::out);
			if (!f) {
				log << "Exporting failed.\n";
				log << "Unable to open file " << file_name_builder.str() << " to export eigenvector " << i << "\n";
				return;
			}
	 		for (unsigned j = 0; j < gmm::vect_size(it->second); ++j)
	   			f << it->second[j] << "\n";
			eig << eigen << "\n";
			i++;
			f.close();
		}
		eig.close();
		// Exports log file if files are created successfully.
		log << "Exporting successful.\n";
		for (auto it = log_data.begin(); it != log_data.end(); it++) {
			log << it->first << it->second << "\n";
		}
		// Exports mesh and mesh_fem objects to text files.
		meshg.write_to_file(dir_name_builder.str() + "/mesh.mh" );
	 	// when the 2nd arg is true, the mesh is saved with the |mf|
		mf_Ug.write_to_file(dir_name_builder.str() + "/solution.mf", true);
		return;
	}
}// end of namespace
