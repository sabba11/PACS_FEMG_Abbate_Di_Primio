#include "descr_qg.hpp"

#include <iostream>

namespace getfem {

	void
	descr_qg::import_all(ftool::md_param & fname)
	{
		FILE_ = fname;
		MESH_FILEG  = FILE_.string_value("MESH_FILEG", "1D points file");
		MESH_TYPEG  = FILE_.string_value("MESH_TYPEG", "1D mesh type");
		FEM_TYPEG   = FILE_.string_value("FEM_TYPEG", "FEM type - solution");
		FEM_TYPEG_DATA = FILE_.string_value("FEM_TYPEG_DATA", "FEM type - coefficients");
		IM_TYPEG 	= FILE_.string_value("IM_TYPEG","Name of integration method");
		import(fname);
		return;
	}

	void
	descr_qg::print_all(void) const
	{
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " GetFEM++ descriptors:                            " << std::endl;
		std::cout << "--------------------------------------------------" << std::endl;
		std::cout << " MESH FILE (GRAPH)         : " << MESH_FILEG  << std::endl;
		std::cout << " FEM TYPE (SOLUTION)       : " << FEM_TYPEG   << std::endl;
		std::cout << " FEM TYPE (COEFFICIENTS)   : " << FEM_TYPEG_DATA << std::endl;
		std::cout << " IM TYPE (GRAPH)           : " << IM_TYPEG	  << std::endl;
		print();
		std::cout << "--------------------------------------------------" << std::endl;
		return;
	}

}
