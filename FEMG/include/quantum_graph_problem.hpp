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
// Abstract class of a generic differential problem defined on a graph
class quantum_graph_problem {
public:
	//Sets up problem parameters and data from input file.
	virtual void init(int argc, char *argv[]) = 0;

	//Assembles problem matrices and vectors.
	virtual void assembly(void) = 0;

	//Solves the differential problem. Returns true if successful.
	virtual bool solve(void) = 0;

	//Exports the solution for external postprocessing.
	virtual void sol_export(void) = 0;

	//Virtual destructor.
	virtual ~quantum_graph_problem() {}

protected:
	//Protected constructor: not public since abstract class objects cannot
	//be instantiated, links the mesh to the finite element and integration methods
	quantum_graph_problem(void) :
		mimg(meshg), mf_Ug(meshg), mf_coeffg(meshg) {}

	// VARIABLES FOR ALL MESHES, INTEGRATION METHODS AND FEMS

	// Mesh on the graph:
	mesh meshg;

	// Integration Method used on the edges of the graph
	mesh_im mimg;

	// Finite Element Method on the function u
	mesh_fem mf_Ug;

	// Finite Element Method for PDE coefficients defined on the graph
	mesh_fem mf_coeffg;

	// Finite Element Method for coefficient from branch to branch
	std::vector<mesh_fem> mf_coeffbranchg;


	// INPUT DATA

	// Input File
	ftool::md_param INPUT;


	// GRAPH PARAMETERS

	// Dimensions of the points
  	unsigned dim_prob;

	// Tangent vectors
  	std::vector<vector_type> tg_vectors;

	// Number of vertices in each branch
	vector_size_type n_vertices;

	// Number of branches
	size_type n_branches;

	// Number of original vertices
	size_type n_origvert;

	// Number of total (original and extended) vertices
	size_type n_totalvert;

	// Discretization parameter vector (one for each branch)
 	// vector_type mesh_step;

	// Possibly useful to save branch lengths?


	// OBJECTS OF THE GRAPH

	// List of BC nodes of the graph
	std::vector<node> BCg;


	// MATRIXES AND VECTORS

	// Matrix for the LHS of the discrete problem
	sparse_matrix_type A;

	// Array of unknowns for the discrete problem
	vector_type U;

  // Vector of radii
  vector_type radii;

	// AUX METHOD FOR INIT

	// Import algorithm specification
	void import_data(void);

	// Build the mesh on the imported DATA
	void build_mesh(void);

	// Set Finite Element Method and Integration Method
	void set_im_and_fem(void);

	// Build problem Parameters
	void build_param(void);

	// Build the lists of the data of the vertices
	// void build_vertices_lists(void);

}; /* end of class quantum_graph_problem */

} /* end of namespace */

#endif
