#ifndef FEMG_DESCR_QG_HPP
#define FEMG_DESCR_QG_HPP

#include <string>
#include <getfem/bgeot_ftool.h>

namespace getfem {
	//! Struct to contain descriptors common to all kind of problems.
	/*!
		The descr_qg struct is an abstract struct containing all the GetFEM++
		descriptors. It provides methods to import and print them, to be inherited.
		Inherit from this class to define problem-specific descriptors.
	*/
	struct descr_qg {
		//! Path to .pts mesh file.
		std::string MESH_FILEG;

		//! GetFEM++ mesh descriptor.
		std::string MESH_TYPEG;

		//! GetFEM++ FEM descriptor.
		std::string FEM_TYPEG;

		//! GetFEM++ FEM descriptor for coefficients.
		std::string FEM_TYPEG_DATA;

		//! GetFEM++ integration method descriptor.
		std::string IM_TYPEG;

		//! Input file from which descriptors are read.
		ftool::md_param FILE_;

		//! Virtual method to print descriptors for derived classes.
		virtual void print(void) const = 0;

		//! Virtual method to import descriptors for derived classes.
		/*!
			\param fname the name of the file.
		*/
		virtual void import(ftool::md_param & fname) = 0;

		//! Method to print all descriptors to screen.
		void print_all(void) const;

		//! Method to import all descriptors from file.
		/*!
			\param fname the name of the file.
		*/
		void import_all(ftool::md_param & fname);

	};

	//! Overload of operator<< to print descriptors.
	/*!
		\param out reference to ostream.
		\param descr the descriptor to be printed.
	*/
	inline std::ostream &
	operator<< (std::ostream & out, const descr_qg & descr)
	{
		descr.print_all();
		return out;
	}
}
#endif
