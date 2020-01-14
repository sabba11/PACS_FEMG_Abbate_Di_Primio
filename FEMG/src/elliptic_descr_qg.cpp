/*!
	\file elliptic_descr_qg.cpp
	\author Stefano Abbate
	\author Andrea Di Primio
	\brief Source file containing definitions of methods from the elliptic_descr_qg derived class.
*/

#include "elliptic_descr_qg.hpp"

namespace getfem {

	void
	elliptic_descr_qg::import(ftool::md_param & fname)
	{
		COMP_METHOD  = FILE_.string_value("COMP_METHOD", "Computational method");
		if (COMP_METHOD != "LU"){
			TOL = FILE_.real_value("TOL", "Residue tolerance");
			BY_ITERATION = FILE_.int_value("BY_ITERATION", "Import Maximum iteration Flag");
			if (BY_ITERATION)
				ITER = FILE_.int_value("ITER", "Maximum iteration");
		}
		GMM_ASSERT1(elliptic_check_validity(), "Invalid COMP_METHOD, TOL or ITER in .param file.");
		return;
	}

	bool
	elliptic_descr_qg::elliptic_check_validity(void) const
	{
		bool check1 = (COMP_METHOD == "LU" || COMP_METHOD == "CG" || COMP_METHOD == "GMRES" || COMP_METHOD == "QMR");
		bool check2 = true;
		if (COMP_METHOD != "LU")
			check2 = (TOL > 0);
		bool check3 = true;
		if (COMP_METHOD != "LU" && BY_ITERATION)
				check3 = (ITER > 0);
		return (check1 & check2 & check3);
	}

	void
	elliptic_descr_qg::print(void) const
	{
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " Elliptic problem specifications:                    " << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " COMPUTATIONAL METHOD      : " << COMP_METHOD << std::endl;
		if (COMP_METHOD != "LU")
			std::cout << " TOLERANCE                 : " << TOL << std::endl;
		if (COMP_METHOD != "LU" && BY_ITERATION)
			std::cout << " ITERATIONS                 : " << ITER << std::endl;
		return;
	}

}
