// Source file containing definitions for class eigen_problem.

// GetFEM++ libraries
#include <getfem/getfem_assembling.h>
#include <gmm/gmm_condition_number.h>
#include <gmm/gmm_except.h>
#include <gmm/gmm_iter_solvers.h>

// Standard libraries
#include <algorithm>
#include <map>

// Project headers
#include "eigen_problem.hpp"

namespace getfem {
void eigen_problem::assembly(void) {
	#ifdef FEMG_VERBOSE_
	std::cout << "[eigen_problem] Assembling matrices..." << std::endl;
	#endif

	assembly_matL();
	if (descr.OPERATOR == "Hamiltonian" || descr.RESCALING_FACTOR == "mass_matrix")
		assembly_matM();
	assembly_matA();

	return;
}

void eigen_problem::assembly_matL() {
	#ifdef FEMG_VERBOSE_
	std::cout << "[eigen_problem] Assembling stiffness matrix L..." << std::endl;
	#endif

	L.resize(n_totalvert, n_totalvert);
	vector_type coeff_val(n_totalvert, 1);
	getfem::asm_stiffness_matrix_for_laplacian(L, mimg, mf_Ug, mf_coeffg, coeff_val);
	gmm::clean(L, 1E-10);

	return;
}

void eigen_problem::assembly_matM() {
	#ifdef FEMG_VERBOSE_
	std::cout << "[eigen_problem] Assembling mass matrix M..." << std::endl;
	#endif

	M.resize(n_totalvert, n_totalvert);
	getfem::asm_mass_matrix(M, mimg, mf_Ug, mf_coeffg);
	gmm::clean(M, 1E-10);

	return;
}

void eigen_problem::assembly_matA() {
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

bool eigen_problem::solve(void) {
	if (descr.COMP_METHOD == "QR") {
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Solving the problem via QR algorithm..." << std::endl;
		#endif
		try {
			dense_matrix_type eigvects(n_totalvert, n_totalvert);
			// Initial guess for all eigenvalues is 1.
			vector_type eigvals(n_totalvert, 1);

			// Check if this can be done a priori.
			if (descr.RESCALING_FACTOR == "mesh_step") {
				gmm::scale(A, 1/mesh_step);
			}
			else if (descr.RESCALING_FACTOR == "mass_matrix") {
				dense_matrix_type inverse_mass(n_totalvert, n_totalvert);
				dense_matrix_type aux_H(n_totalvert, n_totalvert);
				gmm::copy(M, inverse_mass);
				gmm::copy(A, aux_H);
				scalar_type cond_number = gmm::condition_number(M);
				std::cout << "[eigen_problem] The mass matrix M has condition number " << cond_number << std::endl;
				gmm::lu_inverse(inverse_mass);
				gmm::mult(inverse_mass, aux_H, A);
			}

			else
				std::cout << "[eigen_problem] Invalid rescaling factor descriptor." << std::endl;

			scalar_type tol = descr.TOL;
			gmm::implicit_qr_algorithm(A, eigvals, eigvects, tol);
			for (int i = 0; i < n_totalvert; i++) {
				vector_type aux_v;
				for (int j = 0; j < n_totalvert; j++)
					aux_v.push_back(eigvects(j, i));
				auto eigpair = std::make_pair(eigvals[i], aux_v);
				eigpairs.insert(eigpair);
			}
			// Temporary printing routine.
			auto it = eigpairs.begin();
			for(int k = 0; k < 6; k++) {
				std::cout << "Eigenvalue: " << it->first << std::endl;
				std::cout << "Eigenvector: " << std::endl;
				for (auto it2 : it->second)
					std::cout << it2 << std::endl;
				it++;
			}
		}
		GMM_STANDARD_CATCH_ERROR;
	}

	else if (descr.COMP_METHOD == "pencil_thm") {
		// to be defined
	}

	else {
		std::cout << "[eigen_problem] Invalid computational method." << std::endl;
	}

	return true;
}

void eigen_problem::sol_export(const std::string & suff) {
// Temporary code.
	std::cout << "[eigen_problem] Printing eigenvalues..." << std::endl;
	for(auto it = eigpairs.begin(); it != eigpairs.end(); it++)
	   std::cout << it->first << std::endl;
	// Number of vertices in the graph without the mesh
	std::cout << "Number of original vertices: " << n_origvert << std::endl;
	// Number of vertices in the extended graph
	std::cout << "Number of total vertices: " << n_totalvert << std::endl;

	// Exporting routine to MATLAB interface.
	// Executing this code requires GetFEM++-MATLAB interface.
 	unsigned i = 0;
 	for (auto it = eigpairs.begin(); it != eigpairs.end(); it++) {
		std::ostringstream file_name_builder;
		std::string filename = descr.MESH_FILEG;
		size_type idx = filename.rfind('/');
		std::string graph_name = "";
		while (filename[idx] != '.') {
			graph_name += filename[idx];
			idx++;
		}
		//boost::filesystem::create_directory("export");
 		file_name_builder << "export/eigenvector-" << it->first << "-" << i + 1 << "_" << descr.COMP_METHOD << ".U";
		#ifdef FEMG_VERBOSE_
		std::cout << "[eigen_problem] Exporting eigenvector file " << i + 1 << "..." << std::endl;
		#endif
 		std::fstream f(file_name_builder.str(), std::ios::out);
		if (!f)
			std::cerr << "Error opening file " << i + 1 << std::endl;
 		for (unsigned j = 0; j < gmm::vect_size(it->second); ++j) {
   			f << it->second[j] << "\n";
		}
		i++;
		f.close();
	}

 	// when the 2nd arg is true, the mesh is saved with the |mf|
	mf_Ug.write_to_file("export/solution.mf", true);
	return;
}
}// end of namespace
