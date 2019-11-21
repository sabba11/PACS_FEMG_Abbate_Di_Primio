// Header for the derived class to describe a generalized eigenvalue problem.

#ifndef FEMG_ELLIPTIC_PROBLEM_HPP
#define FEMG_ELLIPTIC_PROBLEM_HPP

#include "quantum_graph_problem.hpp"
#include "elliptic_descr_qg.hpp"

#include <functional>

namespace getfem {
class elliptic_problem : public quantum_graph_problem {
public:
	//Constructor.
	elliptic_problem(void) : quantum_graph_problem() {}

	//Sets up problem parameters and data from input file.
	virtual void init(int argc, char *argv[]) override;

	//Builds coefficient vectors
	void set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec, const unsigned & n_mean_points);
	void set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec);

  //Builds source term
  void set_source(const function_type & f, const unsigned & n_mean_points);
  void set_source(const function_type & f);

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
	elliptic_descr_qg descr;

	// Laplacian matrix of the extended graph
	sparse_matrix_type L;

	// Potential mass matrix
	sparse_matrix_type V;

	// Weight matrix for mass
	vector_type F_source;

	// Weight matrix for Laplacian
	vector_type weights;

	// Potential matrix
  vector_type potential;

	// Solution
  vector_type U;

	// Iterative solvers convergence infos
	int n_iteration;
	bool converged_by_tol;

	// Conditioning number for SuperLU;
	double condest;

	// Vector structure to save log data
	std::vector<std::pair<std::string, scalar_type>>  log_data;

	//Auxiliary methods for init procedure
	void import_data(void);
	void build_mesh(void);
	void set_im_and_fem(void);
	void build_param(void);

	//Auxiliary methods to build coefficients
	scalar_type compute_circular_mean(const unsigned & n_mean_points, const scalar_type & radius, const base_node & point, const vector_type & tg_vector, const std::function<scalar_type(base_node)> & f);

	//Auxiliary methods for assembling procedure
	void assembly_matA(void);
	void assembly_matL(void);
};
} //end of namespace
#endif
