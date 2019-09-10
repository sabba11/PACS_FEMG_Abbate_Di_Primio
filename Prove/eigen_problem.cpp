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
	std::cout << "[eigen_problem] Assembling matrices... 0/2 completed" << std::endl;
	#endif

	assembly_matA();
	assembly_matM();

	return;
}

void eigen_problem::assembly_matA() {
	A.resize(n_totalvert, n_totalvert);
	vector_type coeff_val(n_totalvert, 1);
	getfem::asm_stiffness_matrix_for_laplacian(A, mimg, mf_Ug, mf_coeffg, coeff_val);

	#ifdef FEMG_VERBOSE_
	std::cout << "[eigen_problem] Assembling matrices... 1/2 completed" << std::endl;
	#endif

	return;
}

void eigen_problem::assembly_matM() {
	M.resize(n_totalvert, n_totalvert);
	getfem::asm_mass_matrix(M, mimg, mf_Ug, mf_coeffg);

	#ifdef FEMG_VERBOSE_
	std::cout << "[eigen_problem] Assembling matrices... 2/2 completed" << std::endl;
	#endif

	return;
}

bool eigen_problem::solve(void) {
	// Solving via inversion of the mass matrix.
	#ifdef FEMG_VERBOSE_
	std::cout << "[eigen_problem] Solving the problem via mass matrix inversion..." << std::endl;
	#endif
	try {
		//Same caveat as above.
		gmm::clean(A, 1E-10);
		gmm::clean(M, 1E-10);
		dense_matrix_type inverse_mass(n_totalvert, n_totalvert);
		gmm::copy(M, inverse_mass);
		scalar_type cond_number = gmm::condition_number(M);
		std::cout << "[eigen_problem] The mass matrix M has condition number " << cond_number << std::endl;
		gmm::lu_inverse(inverse_mass);
		dense_matrix_type H (n_totalvert, n_totalvert);
		gmm::add(A, M, H);
		eig_M.resize(n_totalvert, n_totalvert);
		gmm::mult(inverse_mass, H, eig_M);
		double tol = 1E-6;
		// Initial guess for all eigenvalues is 1.
		vector_type aux(n_totalvert, 1);
		eigvals.swap(aux);
		eigvects.resize(n_totalvert, n_totalvert);
		gmm::scale(A, 5);
		gmm::implicit_qr_algorithm(A, eigvals, eigvects, tol);
		for (int i = 0; i < n_totalvert; i++) {
			vector_type aux_v;
			for (int j = 0; j < n_totalvert; j++)
				aux_v.push_back(eigvects(j, i));
			auto eigpair = std::make_pair(eigvals[i], aux_v);
			eigpairs.insert(eigpair);
		}
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
	return true;
}

void eigen_problem::sol_export(const std::string & suff) {
// Temporary code.
	std::cout << "[eigen_problem] Printing eigenvalues..." << std::endl;
	for(int i = 0; i < eigvals.size(); ++i)
	   std::cout << eigvals[i] << std::endl;
	std::cout << "origvert " << n_origvert << std::endl;
	std::cout << "totalvert " << n_totalvert << std::endl;

	// Exporting routine to MATLAB interface.
 	unsigned i = 0;
 	for (auto it = eigpairs.begin(); it != eigpairs.end(); it++) {
		std::ostringstream file_name_builder;
 		file_name_builder << "export/eigenvector-" << it->first << "-" << i << ".U";
 		std::fstream f(file_name_builder.str(), std::ios::out);
		if (!f)
			std::cerr << "Error opening file " << i << std::endl;
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
