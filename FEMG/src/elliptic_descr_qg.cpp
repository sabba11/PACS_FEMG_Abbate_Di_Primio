#include "elliptic_descr_qg.hpp"

namespace getfem {

	void
	elliptic_descr_qg::import(ftool::md_param & fname)
	{
    IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS", "Import Radius Flag");
    RFILE  = FILE_.string_value("RFILE", "Input file for radii");
		COMP_METHOD  = FILE_.string_value("COMP_METHOD", "Computational method");
		OUTPUT = FILE_.string_value("OUTPUT","Output directory");
		if (COMP_METHOD != "LU"){
			TOL = FILE_.real_value("TOL", "Residue tolerance");
			BY_ITERATION = FILE_.int_value("BY_ITERATION", "Import Maximum iteration Flag");
			if (BY_ITERATION)
				ITER = FILE_.real_value("ITER", "Maximum iteration");
		}

		if (!elliptic_check_validity())
			std::cerr << "Invalid COMP_METHOD or TOL in .param file." << std::endl;
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
		std::cout << " Eigen problem specifications:                    " << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " COMPUTATIONAL METHOD      : " << COMP_METHOD << std::endl;
		if (COMP_METHOD != "LU")
			std::cout << " TOLERANCE                 : " << TOL << std::endl;
		if (COMP_METHOD != "LU" && BY_ITERATION)
			std::cout << " ITERATIONS                 : " << ITER << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::string s = (IMPORT_RADIUS == 1)? "Yes" : "No";
		std::cout << " IMPORT RADIUS DATA?       : " << s << std::endl;
		if (s == "Yes")
			std::cout << " RADIUS DATA FILE          : " << RFILE << std::endl;
		return;
	}

}
