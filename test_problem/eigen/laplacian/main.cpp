#include <eigen_problem.hpp>

using namespace getfem;

int main (int argc, char* argv[]) {

	try {
		eigen_problem EP;

		//Initialize the problem
		EP.init(argc, argv);

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
