// Prima bozza dell'header del problema, la sto scrivendo a partire da problema 3d1d

//L'idea è di scrivere una classe padre che non sappiamo ancora se sarà
//astratta che ci permetta di avere come figli il problema agli autovalori
// e i problemi ellittici, in tal modo se mai dovessimo programmare anche
// per un accoppiato potremmo averlo come figlio
#ifndef FEMG_QUANTUMGRAPHPROBLEM_HPP
#define FEMG_QUANTUMGRAPHPROBLEM_HPP

//GetFEM++ libraries
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>
//#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_mesh_im.h> //was not present and the code worked...

//Standard libraries
#include <vector>

//Project headers
#include "type_aliases.hpp"
#include "mesh1d_prova_new_order.hpp"
#include "node.hpp"
#include "descr_qg.hpp"
//#include "param_qg.hpp" to be defined

namespace getfem {
// Abstract class of a generic differential problem defined on a graph
class quantum_graph_problem {
public:
	//Sets up problem parameters and data from input file.
	void init(int argc, char *argv[]);

	//Assembles problem matrices and vectors.
	virtual void assembly(void) = 0;

	//Solves the differential problem. Returns true if successful.
	virtual bool solve() = 0;

	//Exports the solution for external postprocessing.
	virtual void sol_export(const std::string & suff = "") = 0;

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

	// Algorithm description strings (mesh files, FEM types, solver info...)
	descr_qg descr;


	// GRAPH PARAMETERS

	// Number of vertices in each branch
	vector_size_type n_vertices;

	// Number of branches
	size_type n_branches;

	// Number of original vertices
	size_type n_origvert;

	// Number of total (original and extended) vertices
	size_type n_totalvert;

	// Discretization parameter
	scalar_type mesh_step;


	// OBJECTS OF THE GRAPH

	// List of BC nodes of the graph
	std::vector<node> BCg;

	// List of original vertices of the graph
	// std::vector<node> OGg;

	// List of extended vertices of the graph
	// std::vector<node> EXg;


	// MATRIXES AND VECTORS

	// Matrix for the LHS of the discrete problem
	sparse_matrix_type A;

	// Array of unknowns for the discrete problem
	vector_type U;


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
	void build_vertices_lists(void);

}; /* end of class quantum_graph_problem */

} /* end of namespace */

#endif
