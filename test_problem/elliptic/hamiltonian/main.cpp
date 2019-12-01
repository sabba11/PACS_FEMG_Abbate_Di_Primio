#include <elliptic_problem.hpp>

using namespace getfem;

scalar_type f_source(const base_node & p);
scalar_type f_potential(const base_node & p);

int main (int argc, char* argv[]) {

	try {
		elliptic_problem EP;

		//Initialize the problem
		EP.init(argc, argv);

		// Setting source
		EP.set_source(f_source);

		// Setting potential for building Hamiltonian
		std::vector<std::string> s{"potential"};
		std::vector<function_type> fv{f_potential};
		EP.set_coefficients(fv,s);

		//Assembling the problem
		EP.assembly();

		//Solving the problem
		bool success = EP.solve();

		//Export the solution (test vrs.), just prints eigvals
		EP.sol_export();
	}

	GMM_STANDARD_CATCH_ERROR;
	return 0;
}

// Source simply at 1
scalar_type f_source(const base_node & p){
	return 1;
}
// Potential at 1 for building Hamiltonian
scalar_type f_potential(const base_node & p){
	return 1;
}
