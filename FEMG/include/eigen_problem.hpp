// Header for the derived class to describe a generalized eigenvalue problem.

#ifndef FEMG_EIGEN_PROBLEM_HPP
#define FEMG_EIGEN_PROBLEM_HPP

#include "quantum_graph_problem.hpp"
#include "eigen_descr_qg.hpp"
#include <map>

namespace getfem {
class eigen_problem : public getfem::quantum_graph_problem {
public:
	//Constructor.
	eigen_problem(void) : quantum_graph_problem() {}

	//Sets up problem parameters and data from input file.
	virtual void init(int argc, char *argv[]) override;

	//Assembles matrices and vectors for the problem.
	virtual void assembly(void) override;

	//Solves the problem. Returns true if successful.
	virtual bool solve(void) override;

	// Possibly other solve methods here, or solve(void) should eventually
	// contain all solve methods and select one according to descriptors

	//Exports the solution for external postprocessing.
	virtual void sol_export(void) override;

private:
	// Algorithm descriptors
	eigen_descr_qg descr;

	// Mass matrix for the discrete problem
	sparse_matrix_type M;

	// Laplacian matrix of the extended graph
	sparse_matrix_type L;

	// Map to save eigenpairs
	std::multimap< scalar_type, vector_type > eigpairs;

	//Auxiliary methods for init procedure
	void import_data(void);
	void build_mesh(void);
	void set_im_and_fem(void);
	void build_param(void);
	//void build_vertices_lists(void);

	//Auxiliary methods for assembling procedure
	void assembly_matA(void);
	void assembly_matM(void);
	void assembly_matL(void);
};
} //end of namespace
#endif
