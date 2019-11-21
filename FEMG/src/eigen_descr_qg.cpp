#include "eigen_descr_qg.hpp"

namespace getfem {

	void
	eigen_descr_qg::import(ftool::md_param & fname)
	{
    IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS", "Import Radius Flag");
    RFILE  = FILE_.string_value("RFILE", "Input file for radii");
		COMP_METHOD  = FILE_.string_value("COMP_METHOD", "Computational method");
		OUTPUT = FILE_.string_value("OUTPUT","Output directory");
		TOL = FILE_.real_value("TOL", "Residue tolerance");
		if (!eigen_check_validity())
			std::cerr << "Invalid COMP_METHOD or TOL in .param file." << std::endl;
		return;
	}

	bool
	eigen_descr_qg::eigen_check_validity(void) const
	{
		bool check1 = (COMP_METHOD == "QR" || COMP_METHOD == "QZ");
		bool check2 = (TOL > 0);
		return (check1 & check2);
	}

	void
	eigen_descr_qg::print(void) const
	{
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " Eigen problem specifications:                    " << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " COMPUTATIONAL METHOD      : " << COMP_METHOD << std::endl;
		std::cout << " TOLERANCE                 : " << TOL << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::string s = (IMPORT_RADIUS == 1)? "Yes" : "No";
		std::cout << " IMPORT RADIUS DATA?       : " << s << std::endl;
		if (s == "Yes")
			std::cout << " RADIUS DATA FILE          : " << RFILE << std::endl;
		return;
	}

}
