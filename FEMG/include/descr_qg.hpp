/*!
	\file descr_qg.hpp
	\author Stefano Abbate
	\author Andrea Di Primio
	\brief Header containing the descr_qg abstract class.
*/

#ifndef FEMG_DESCR_QG_HPP
#define FEMG_DESCR_QG_HPP

#include <string>
#include <getfem/bgeot_ftool.h>

namespace getfem {
	//! Struct to contain descriptors common to all kinds of problems.
	/*!
		The descr_qg struct is an abstract struct containing all the GetFEM++
		descriptors as std::strings and a ftool::md_param object. The latter is
		provided by GetFEM++ to handle data importing and it is not a
		descriptor.
		To define problem-specific descriptors, create a child struct.
		The class provides protected methods to import and print all descriptors,
		as well as virtual methods to define problem-specific descriptors.
		The descriptors should all be contained in the same input file.
	*/
	struct descr_qg {
	public:
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

		//! Boolean flag to import radii data (not a GetFEM++ descriptor).
		/*!
			If IMPORT_RADIUS = 1, each branch of the graph is intended as axis of a
			cylinder with the specified radius. This may only affect the computation
			of known coefficients. See eigen_problem.hpp for details.
		*/
		bool IMPORT_RADIUS;

		//! Path to file containing radii data.
		/*!
			Only relevant if IMPORT_RADIUS = 1. It should be a list of positive, real numbers.
		*/
		std::string RFILE;

		//! Path to output directory for export routines.
		std::string OUTPUT;

		//! Input file from which descriptors are read.
		ftool::md_param FILE_;

		//! Method to import all descriptors from file.
		/*!
			Call this method to import all descriptors. It calls import().
			\param fname file containing input data.
			\warning fname should contain a list of assignments of the type TAG = VALUE,
					 where TAG is a descriptor name (use the std::string labels as TAGs).
					 ftool::md_param FILE_ is not a descriptor.
		*/
		void import_all(ftool::md_param & fname);

		//! Overload of operator<< to print descriptors.
		/*!
			Same as print_all(). Use this operator to print desciptors.
			\return reference to ostream.
			\param out reference to ostream.
			\param descr the descriptor to be printed.
		*/
		friend inline std::ostream &
		operator<< (std::ostream & out, const descr_qg & descr)
		{
			descr.print_all();
			return out;
		}

	protected:
		//! Method to print all descriptors to screen, called by operator<<.
		/*!
			Use operator<< to print all descriptors.
		*/
		void print_all(void) const;

		//! Virtual method to print descriptors for derived classes.
		virtual void print(void) const = 0;

		//! Virtual method to import descriptors for derived classes.
		/*!
			The derived class should have as public objects all the custom
			descriptors. To import them, use the ftool::md_param methods
			real_value, string_value, int_value, defined in the bgeot_ftool.h header.
			They have the structure return_type type_value(TAG, description), where
			description is a std::string briefly explaining what the variable represents.
			\param fname the name of the ftool::md_param object.
			\warning fname should contain a list of assignments of the type TAG = VALUE,
					 where TAG is a descriptor name (use the std::string labels as TAGs).
					 ftool::md_param FILE_ is not a descriptor.

		*/
		virtual void import(ftool::md_param & fname) = 0;
	};
}
#endif
