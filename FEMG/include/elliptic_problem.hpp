#ifndef FEMG_ELLIPTIC_PROBLEM_HPP
#define FEMG_ELLIPTIC_PROBLEM_HPP

#include "quantum_graph_problem.hpp"
#include "elliptic_descr_qg.hpp"

#include <functional>

namespace getfem {
	//! Derived class to define an elliptic differential problem.
	/*!
		The elliptic_problem class, child of quantum_graph_problem, defines a differential
	 	problem of the kind
		\f[
			 -\mathrm{div}(a(x)\nabla u) + v(x)u(x) = f(x)
		\f]
		where a, v, f are given scalar functions defined on the whole graph.
		All coefficients should respect well-posedness costraints, namely:
	*/
	class elliptic_problem : public quantum_graph_problem {
	public:
		//! Constructor. Using quantum_graph_problem constructor to link GetFEM++ objects.
		elliptic_problem(void) : quantum_graph_problem() {}

		//! Method to evaluate known coefficients on the extended graph nodes.
		/*!
			The function evaluates the elements of f_vec on the points of the mesh and fills data structures
			according to the keywords in s_vec.
			\param f_vec a vector of functions.
			\param s_vec a vector of keywords to identify elements of f_vec, possible keywords: left, right, potential for a, p and v respectively.
		*/
		void set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec);

		//! Overloaded method to evaluate known coefficients on the extended graph nodes.
		/*!
			The function evaluates data on the points of the mesh and fills data structures
			according to the keywords in s_vec. The computation is made thinking of every
			branch as an axis of a cylinder. The value of the coefficients in a point is given by the
			arithmetic mean of the values of the elements of f_vec in n_mean_points points equally
			distributed along the cylinder's circular cross section taken at that point.
			\warning IMPORT_RADIUS should be 1 to use this method.
			\param f_vec a vector of functions.
			\param s_vec a vector of keywords to identify elements of f_vec, possible keywords: left, right, potential for a, p and v respectively.
			\param n_mean_points number of points in each circle.
		*/
		void set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec, const unsigned & n_mean_points);

	  	//! Additional method to evaluate the source term on the extended graph nodes.
		/*!
			The function evaluates f on the points of the mesh and fills a data structure
			accordingly.
			\param f an std::function.
		*/
		void set_source(const function_type & f);

		//! Overloaded method to evaluate the source term on the extended graph nodes.
		/*!
			The function operates in the same fashion as the overloaded set_coefficients
			method.
			\param f an std::function.
			\param n_mean_points number of points in each circle.
		*/
	  	void set_source(const function_type & f, const unsigned & n_mean_points);

		//! Overridden method to assemble matrices and vectors for the problem.
		virtual void assembly(void) override;

		//! Overridden method to solve the discrete problem.
		/*!
			\return True if successful, False if failed.
		*/
		virtual bool solve(void) override;

		//! Overridden method to export the solution for external postprocessing.
		virtual void sol_export(void) override;

	private:
		/*
		+------------------------------------------------------+
		| 1. Descriptors struct 							   |
		+------------------------------------------------------+
		*/
		//! Struct to contain all descriptors.
		elliptic_descr_qg descr;

		/*
		+------------------------------------------------------+
		| 2. Additional structures for the discrete problem    |
		+------------------------------------------------------+
		*/
		//! Laplacian matrix of the extended graph.
		sparse_matrix_type L;

		//! Potential matrix.
		sparse_matrix_type V;
		
		//! Source term at rhs of the differential problem.
		vector_type F_source;

		//! Weight vector for Laplacian (function a(x)).
		vector_type weights;

		//! Potential vector (function v(x)).
	  	vector_type potential;

		// Solution
	  	vector_type U;

		/*
		+------------------------------------------------------+
		| 3. Numerical solvers variables					   |
		+------------------------------------------------------+
		*/
		//! Number of iterations to achieve convergence (at most ITER).
		int n_iteration;

		//! Boolean to assert if a numerical method achieved convergence.
		bool converged_by_tol;

		//! Condition number (for SuperLU);
		scalar_type condest;

		/*
		+------------------------------------------------------+
		| 4. Auxiliary methods								   |
		+------------------------------------------------------+
		*/
		// 4.1 Auxiliary methods for init procedure
		//! Overridden method to import data from input file.
		virtual void import_data(void) override;

		//! Overridden method to read .pts file to create the mesh.
		virtual void build_mesh(void) override;

		//! Overridden method to set Finite Element Method and Integration method.
		virtual void set_im_and_fem(void) override;

		//! Overridden method to build other problem parameters (for example, radii).
		virtual void build_param(void) override;

		//! Overridden method to set default values for known coefficients.
		virtual void set_default_coefficients(void) override;

		// 4.2 Auxiliary methods to compute coefficients
		//! Auxiliary method to compute means on cross sections.
		/*!
			\return value of coefficient at given point through mean on cross section.
			\param n_mean_points number of points in each circle.
			\param radius radius of the branch.
			\param point coordinates of the point.
			\param tg_vector direction of the branch.
			\param f coefficient function.
		*/
		scalar_type compute_circular_mean(const unsigned & n_mean_points, const scalar_type & radius, const base_node & point, const vector_type & tg_vector, const std::function<scalar_type(base_node)> & f);

		// 4.3 Auxiliary methods for assembly routines
		//! Auxiliary method to assembly the final matrix.
		void assembly_matA(void);

		//! Auxiliary method to assembly the Laplace matrix.
		void assembly_matL(void);
	};
} //end of namespace
#endif
