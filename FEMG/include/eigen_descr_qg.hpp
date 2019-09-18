#ifndef FEMG_EIGEN_DESCR_QG_HPP
#define FEMG_EIGEN_DESCR_QG_HPP

#include "descr_qg.hpp"
#include "type_aliases.hpp"

namespace getfem {

	struct eigen_descr_qg : descr_qg {
		//Operator whose eigenvalues are computed
		std::string OPERATOR;

		//Computational method
		std::string COMP_METHOD;

		//Output directory for export routines
		std::string OUTPUT;

		//Residue tolerance (for iterative solvers)
		scalar_type TOL;

		//Imports eigen_problem specific descriptors
		virtual void import(ftool::md_param & fname) override;

		//Checks validity of non-GetFEM++ descriptors
		bool eigen_check_validity(void) const;

		//Prints eigen_problem specific descriptors
		virtual void print(void) const override;
	};

}
#endif
