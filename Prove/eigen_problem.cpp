// Source file containing definitions for class eigen_problem.

// GetFEM++ libraries
#include <getfem/getfem_assembling.h>
#include <gmm/gmm_condition_number.h>
#include <gmm/gmm_except.h>

// Standard libraries
#include <algorithm>

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
	// Caveat! It still lacks computation of n_origvert and n_extdvert
	size_type n_totalpoints = 21;
	A.resize(n_totalpoints, n_totalpoints);
	vector_type coeff_val(n_totalpoints, 1);
	getfem::asm_stiffness_matrix_for_laplacian(A, mimg, mf_Ug, mf_coeffg, coeff_val);

	#ifdef FEMG_VERBOSE_
	std::cout << "[eigen_problem] Assembling matrices... 1/2 completed" << std::endl;
	#endif

	return;
}

void eigen_problem::assembly_matM() {
	// Caveat! It still lacks computation of n_origvert and n_extdvert
	size_type n_totalpoints = 21;
	M.resize(n_totalpoints, n_totalpoints);
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
		size_type n_totalpoints = 21;
		gmm::clean(A, 1E-10);
		gmm::clean(M, 1E-10);
		dense_matrix_type inverse_mass(n_totalpoints, n_totalpoints);
		gmm::copy(M, inverse_mass);
		scalar_type cond_number = gmm::condition_number(M);
		std::cout << "[eigen_problem] The mass matrix M has condition number " << cond_number << std::endl;
		gmm::lu_inverse(inverse_mass);
		dense_matrix_type H (n_totalpoints, n_totalpoints);
		gmm::add(A, M, H);
		dense_matrix_type eig_M(n_totalpoints, n_totalpoints);
		gmm::mult(inverse_mass, H, eig_M);
		double tol = 1E-6;
		vector_type aux(n_totalpoints, 1);
		eigvals.swap(aux);
		gmm::implicit_qr_algorithm(eig_M, eigvals, tol);
		std::sort(eigvals.begin(), eigvals.end());
	}
	GMM_STANDARD_CATCH_ERROR;
	return true;
}

void eigen_problem::sol_export(const std::string & suff) {
	// Temporary code.
	std::cout << "[eigen_problem] Printing eigenvalues..." << std::endl;
	for(int i = 0; i < eigvals.size(); ++i)
	   std::cout << eigvals[i] << std::endl;
	return;
}

}// end of namespace
