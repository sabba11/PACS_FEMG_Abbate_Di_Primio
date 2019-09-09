// Struct to contain GetFEM++ descriptors.

#ifndef FEMG_DESCR_QG_HPP
#define FEMG_DESCR_QG_HPP

#include <string>

namespace getfem {

	struct descr_qg {
		//Path to meshg file
		std::string MESH_FILEG;

		//Mesh descriptor
		std::string MESH_TYPEG;

		//FEM descriptors for the solution
		std::string FEM_TYPEG;

		//FEM descriptor for the coefficients
		std::string FEM_TYPEG_DATA;

		//IM descriptor
		std::string IM_TYPEG;

		//Output directory
		std::string OUTPUT;

		//Solver descriptor
		std::string SOLVE_METHOD;

		//Max number of iterations (for iterative solvers)
		size_type MAX_ITER;

		//Residue tolerance (for iterative solvers)
		scalar_type TOL;

		//Input file
		ftool::md_param FILE_;

		friend std::ostream & operator << (
		std::ostream & out, const descr_qg & descr
		)
		{
			std::cout << "--------------------------------------------------" << std::endl;
			std::cout << " GetFEM++ descriptors:                            " << std::endl;
			std::cout << "--------------------------------------------------" << std::endl;
			std::cout << " MESH FILE (GRAPH)         : " << descr.MESH_FILEG  << std::endl;
			std::cout << " IM TYPE (GRAPH)           : " << descr.IM_TYPEG	  << std::endl;
			std::cout << " FEM TYPE (SOLUTION)       : " << descr.FEM_TYPEG   << std::endl;
			std::cout << " FEM TYPE (COEFFICIENTS)   : " << descr.FEM_TYPEG_DATA << std::endl;
			std::cout << "--------------------------------------------------" << std::endl;
			std::cout << " SOLVER                    : " << descr.SOLVE_METHOD << std::endl;
			std::cout << " MAX ITERATIONS            : " << descr.MAX_ITER    << std::endl;
			std::cout << " TOLERANCE                 : " << descr.TOL         << std::endl;
			std::cout << "--------------------------------------------------" << std::endl;
		return out;
		} //end of operator<<


		//! Import algorithm specifications from file .param
		void import(ftool::md_param & fname)
		{
			FILE_ = fname;
			MESH_FILEG  = FILE_.string_value("MESH_FILEG", "1D points file");
			MESH_TYPEG  = FILE_.string_value("MESH_TYPEG", "1D mesh type");
			FEM_TYPEG   = FILE_.string_value("FEM_TYPEG", "FEM type - solution");
			FEM_TYPEG_DATA = FILE_.string_value("FEM_TYPEG_DATA", "FEM type - coefficients");
			IM_TYPEG 	= FILE_.string_value("IM_TYPEG","Name of integration method");

			// Set default values if empty strings are passed.
			if(FEM_TYPEG_DATA == "") FEM_TYPEG_DATA = "FEM_PK(1,1)";

			SOLVE_METHOD = FILE_.string_value("SOLVE_METHOD", "Name of solver");
			if (SOLVE_METHOD != "SuperLU") { // iterative solver
				MAX_ITER  = FILE_.int_value("MAX_ITER", "Max number of sub-iterations");
			TOL = FILE_.real_value("TOL"); if (TOL == 0.) TOL = 2.0e-10;
			}
			OUTPUT = FILE_.string_value("OUTPUT","Output Directory");
		} // end of import
	}; // end of descr_qg struct

} // end of namespace
#endif
