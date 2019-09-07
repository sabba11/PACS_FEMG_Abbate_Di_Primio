// Prima bozza dell'header del problema, la sto scrivendo a partire da problema 3d1d

//L'idea è di scrivere una classe padre che non sappiamo ancora se sarà
//astratta che ci permetta di avere come figli il problema agli autovalori
// e i problemi ellittici, in tal modo se mai dovessimo programmare anche
// per un accoppiato potremmo averlo come figlio
#ifndef FEMG_QUANTUMGRAPHPROBLEM_HPP
#define FEMG_QUANTUMGRAPHPROBLEM_HPP

// per ora non metto librerie le aggiungerò via via mentre testo
// o programmo

namespace getfem {
// Abstract class of a generic differential problem defined on a graph
class quantum_graph_problem {
public:
	//Sets up problem parameters and data from input file.
	void init(int argc, char *argv[]);

	//Assembles problem matrices and vectors.
	void assembly(void) = 0;

	//Solves the differential problem. Returns true if successful.
	bool solve() = 0;

	//Exports the solution for external postprocessing.
	void export(const string & suff = "") = 0;

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

	// Physical parameters (dimensionless)
	param_qg param;


	// GRAPH PARAMETERS

	// Number of vertices per branch
	vector_size_type n_vertices;

	// Number of branches
	size_type n_branches;

	//// Number of boundary vertices
	// size_type n_boundvert;
	//
	//// Number of internal vertices
	// size_type n_intvert;

	// Number of original vertices
	size_type n_origvert;

	// Number of extended vertices
	size_type n_extdvert;

	//// Number of degrees of freedom
	// dofg dof;


	// OBJECTS OF THE GRAPH

	// // List of BC nodes of the graph
	// std::vector<node> BCg;
	//
	// // List of junction nodes of the graph
	// std::vector<node> Jg;

	// List of original vertices of the graph
	std::vector<node> OGg;

	// List of extended vertices of the graph
	std::vector<node> EXg;


	// MATRIXES AND VECTORS

	// Matrix for the LHS of the discrete problem
	sparse_matrix_type A;

	// Mass Matrix of the discrete problem
	sparse_matrix_type M;

	// Array of unknowns for the discrete problem
	vector_type U;

	// Source term of the problem (RHS)
	vector_type F;


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

	// Build the matrices and vertices
	void assembly_matA(void);
	void assembly_matM(void);
	void assembly_rhs(void);

}; /* end of class quantum_graph_problem */

} /* end of namespace */

#endif
