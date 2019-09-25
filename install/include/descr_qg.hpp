// Struct to contain GetFEM++ descriptors common to all problems.

#ifndef FEMG_DESCR_QG_HPP
#define FEMG_DESCR_QG_HPP

#include <string>
#include <getfem/bgeot_ftool.h>

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

		//Input file
		ftool::md_param FILE_;

		//Virtual routines for derived classes.
		virtual void print(void) const = 0;
		virtual void import(ftool::md_param & fname) = 0;

		// Imports all descriptors from
		void import_all(ftool::md_param & fname);
		void print_all(void) const;
	};

	inline std::ostream &
	operator<< (std::ostream & out, const descr_qg & descr)
	{
		descr.print_all();
		return out;
	}
}
#endif
