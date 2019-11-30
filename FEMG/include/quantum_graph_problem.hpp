#ifndef FEMG_QUANTUMGRAPHPROBLEM_HPP
#define FEMG_QUANTUMGRAPHPROBLEM_HPP

//GetFEM++ libraries
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_mesh_im.h>

//Standard libraries
#include <vector>

//Project headers
#include "type_aliases.hpp"
#include "mesh_1d.hpp"
#include "node.hpp"

namespace getfem {
	//! Abstract class to define a differential problem on a quantum graph.
	/*!
		The quantum_graph_problem class defines a generic differential problem
		on a class. It is an interface with four pure virtual methods, that should
		be called in the presented order (init, assembly, solve, sol_export) in the
		source code generating the final executable.
		It provides protected data members to save the imported graph's properties,
		it initializes and links GetFEM++ finite elements and integration methods
		to the mesh, it contains data structures to handle the resulting
		discrete problem and finally provides auxiliary methods for init procedures.
		To define a differential problem, inherit from this class.
	*/
	class quantum_graph_problem {
	public:
		//! Virtual method to set up the problem and various parameters.
		/*!
			init, by default, calls a set of pure virtual methods to set up the
			problem.
			\param argc number of arguments passed via terminal execution, should be 2.
			\param argv strings passed via terminal execution, should be the name of the input file.
		*/
		virtual void init(int argc, char *argv[]) {
			//1. Read the .param filename from standard input
   	     	INPUT.read_command_line(argc, argv);

   	     	//2. Import data (algorithm specifications, boundary conditions,...)
   	     	import_data();

   	     	//3. Build mesh for the graph
   	     	build_mesh();

   	     	//4. Set finite elements and integration methods
   	     	set_im_and_fem();

   	     	//5. Build problem parameters
   	     	build_param();

   		 	//6. Set default values for coefficients
   	 	 	set_default_coefficients();

   	     	//7. Build the lists of the data of the vertices
   	     	// build_vertices_lists();

			return;
		}

		//! Method to evaluate known coefficients on the extended graph nodes.
		/*!
			The function evaluates the elements of f_vec on the points of the mesh and fills data structures
			according to the keywords in s_vec.
			\param f_vec a vector of functions.
			\param s_vec a vector of keywords to identify elements of f_vec, possible keywords: left, right, potential for a, p and v respectively.
		*/
		virtual void set_coefficients(const vector_function_type & f_vec, const vector_string_type & s_vec) = 0;
		
		//! Virtual method to assembly problem matrices and vectors.
		virtual void assembly(void) = 0;

		//! Solves the discrete problem.
		/*!
			\return True if successful, False if failed.
		*/
		virtual bool solve(void) = 0;

		//! Exports the solution for external postprocessing.
		virtual void sol_export(void) = 0;

		//! Virtual destructor.
		virtual ~quantum_graph_problem() {}

	protected:
		//! Protected constructor.
		/*!
			The constructor links integration methods and finite elements to the
			mesh meshg. See GetFEM++ documentation for details.
		*/
		quantum_graph_problem(void) :
			mimg(meshg), mf_Ug(meshg), mf_coeffg(meshg) {}
		/*
		+------------------------------------------------------+
		| 1. GetFEM++ objects (mesh, integration method, fems) |
		+------------------------------------------------------+
		*/
		//! Graph mesh.
		mesh meshg;

		//! Integration method used on the edges of the graph.
		mesh_im mimg;

		//! Finite Element Method for the unknown function.
		mesh_fem mf_Ug;

		//! Finite Element Method for known coefficients defined on the graph.
		mesh_fem mf_coeffg;

		//! Finite Element Method for known coefficients from branch to branch.
		std::vector<mesh_fem> mf_coeffbranchg;

		/*
		+------------------------------------------------------+
		| 2. Input file										   |
		+------------------------------------------------------+
		*/
		//! Input file.
		ftool::md_param INPUT;

		/*
		+------------------------------------------------------+
		| 3. Graph properties								   |
		+------------------------------------------------------+
		*/
		//! Dimensions of the problem.
	  	unsigned dim_prob;

		//! Vector of tangent vectors (one for branch).
	  	std::vector<vector_type> tg_vectors;

		//! Number of vertices in each branch.
		vector_size_type n_vertices;

		//! Total number of branches.
		size_type n_branches;

		//! Number of vertices in the original graph.
		size_type n_origvert;

		//! Number of vertices in the extended graph (i.e. including discretization nodes).
		size_type n_totalvert;

		// Discretization parameter vector (one for each branch).
	 	// vector_type mesh_step;

		//! Vector of BC nodes of the graph.
		std::vector<node> BCg;

		//! Vector of branch radii.
		/*!
			Used only if IMPORT_RADIUS = 1. Stores data contained in RFILE.
		*/
		vector_type radii;

		/*
		+------------------------------------------------------+
		| 4. Data structures for the discrete problem          |
		+------------------------------------------------------+
		*/
		//! Matrix for the LHS of the discrete problem.
		sparse_matrix_type A;

		//! Array of unknowns for the discrete problem.
		vector_type U;

		/*
		+------------------------------------------------------+
		| 5. Pure virtual auxiliary methods for init procedure |
		+------------------------------------------------------+
		*/
		//! Method to import data from input file.
		virtual void import_data(void) = 0;

		//! Read .pts file to create the mesh.
		virtual void build_mesh(void) = 0;

		//! Set Finite Element Method and Integration method.
		virtual void set_im_and_fem(void) = 0;

		//! Build other problem parameters (for example, radii).
		virtual void build_param(void) = 0;

		//! Set default values for known coefficients.
		virtual void set_default_coefficients(void) = 0;

		// Build the lists of the data of the vertices
		// void build_vertices_lists(void);

		/*
		+------------------------------------------------------+
		| 6. Data structure for log data					   |
		+------------------------------------------------------+
		*/
		//! Vector structure to save log data
		std::vector<std::pair<std::string, scalar_type>> log_data;

	}; /* end of class quantum_graph_problem */

} /* end of namespace */

#endif
