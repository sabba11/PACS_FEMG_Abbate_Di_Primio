#include "eigen_descr_qg.hpp"

namespace getfem {

	void
	eigen_descr_qg::import(ftool::md_param & fname)
	{
		OPERATOR  = FILE_.string_value("OPERATOR", "Operator");
		COMP_METHOD  = FILE_.string_value("COMP_METHOD", "Computational method");
		OUTPUT = FILE_.string_value("OUTPUT","Output directory");
		TOL = FILE_.real_value("TOL", "Residue tolerance");
		return;
	}

	bool
	eigen_descr_qg::eigen_check_validity(void) const
	{
		bool check1 = (OPERATOR == "Laplacian" || OPERATOR == "Hamiltonian");
		bool check2 = (COMP_METHOD == "QR");
		bool check3 = (TOL > 0);
		return (check1 & check2 & check3);
	}

	void
	eigen_descr_qg::print(void) const
	{
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " Eigen problem specifications:                    " << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " OPERATOR                  : " << OPERATOR    << std::endl;
		std::cout << " COMPUTATIONAL METHOD      : " << COMP_METHOD << std::endl;
		return;
	}

}
