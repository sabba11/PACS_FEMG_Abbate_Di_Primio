// Header for the derived class to describe a generalized eigenvalue problem.

#ifndef FEMG_EIGEN_PROBLEM_HPP
#define FEMG_EIGEN_PROBLEM_HPP

#include "quantum_graph_problem.hpp"
#include <map>

namespace getfem {
class eigen_problem : public getfem::quantum_graph_problem {
public:
	//Constructor
	eigen_problem(void) : quantum_graph_problem() {}

	//Assembles matrices and vectors for the problem.
	virtual void assembly(void) override;

	//Solves the problem. Returns true if successful.
	virtual bool solve(void) override;

	// Possibly other solve methods here, or solve(void) should eventually
	// contain all solve methods and select one according to descriptors

	//Exports the solution for external postprocessing.
	virtual void sol_export(const std::string & suff = "") override;

private:
	// Mass matrix for the discrete problem
	sparse_matrix_type M;

	// Laplacian matrix of the extended graph
	sparse_matrix_type L;

	// Auxiliary matrix (final LHS)
	dense_matrix_type eig_M;

	// Map to save eigenpairs
	std::multimap< scalar_type, vector_type > eigpairs;

	//Auxiliary methods for assembling procedure
	void assembly_matA(void);
	void assembly_matM(void);
	void assembly_matL(void);
};
} //end of namespace
#endif
