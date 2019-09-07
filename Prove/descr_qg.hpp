// Struct to contain GetFEM++ descriptors. To be completed.

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

		//TBA: solve descriptors, import routine

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

		return out;
		} //end of operator<<
	} // end of descr_qg struct

} // end of namespace
#endif
