#ifndef FEMG_ELLIPTIC_DESCR_QG_HPP
#define FEMG_ELLIPTIC_DESCR_QG_HPP

#include "descr_qg.hpp"
#include "type_aliases.hpp"

namespace getfem {

	struct elliptic_descr_qg : descr_qg {
		//File from which we have to import the radii
    std::string RFILE;

		// Flag for importing the radius
		bool IMPORT_RADIUS;

		//Computational method
		std::string COMP_METHOD;

		//Output directory for export routines
		std::string OUTPUT;

		//Residue tolerance (for iterative solvers)
		scalar_type TOL = 0;

		// Bool that say if a maximum iteratio is given
		bool BY_ITERATION;

		//Maximum iteration (for iterative solvers)
		scalar_type ITER = 0;

		//Imports eigen_problem specific descriptors
		virtual void import(ftool::md_param & fname) override;

		//Checks validity of non-GetFEM++ descriptors
		bool elliptic_check_validity(void) const;

		//Prints eigen_problem specific descriptors
		virtual void print(void) const override;
	};

}
#endif
