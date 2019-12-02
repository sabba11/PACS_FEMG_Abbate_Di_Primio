/*!
	\file elliptic_descr_qg.hpp
	\author Stefano Abbate
	\author Andrea Di Primio
	\brief Header containing the elliptic_descr_qg derived class.
*/

#ifndef FEMG_ELLIPTIC_DESCR_QG_HPP
#define FEMG_ELLIPTIC_DESCR_QG_HPP

#include "descr_qg.hpp"
#include "type_aliases.hpp"

namespace getfem {
	//! Struct to contain elliptic problem-specific descriptors.
	/*!
		The elliptic_descr_qg struct contains all the custom
		descriptors for an elliptic differential problem of the kind
		\f[
			 -\mathrm{div}(a(x)\nabla u) + v(x)u(x) = f(x)
		\f]
		where a, f are given functions defined on the whole graph.
		It provides methods to import, check and print them.
		Execute import_all() and print_all() (inherited methods) in init subroutines.
	*/
	struct elliptic_descr_qg : descr_qg {
		//! Numerical method employed in the solve routine.
		/*!
			For the elliptic_problem class, there are four valid computational methods: SuperLU (direct),
			CG, GMRES, QMR (iterative).
		*/
		std::string COMP_METHOD;

		//! Residue tolerance (for iterative solvers).
		scalar_type TOL = 0;

		//! Boolean to know if a maximum iteration is given.
		bool BY_ITERATION;

		//! Maximum number of iterations (for iterative solvers).
		scalar_type ITER = 0;

	protected:
		//! Overridden import() method to define custom descriptors.
		/*!
			The derived class should have as public objects all the custom
			descriptors. To import them, use the ftool::md_param methods
			real_value, string_value, int_value, defined in the bgeot_ftool.h header.
			They have the structure return_type type_value(TAG, description), where
			description is a std::string briefly explaining what the variable represents.
			\param fname fname the name of the ftool::md_param object.
			\warning fname should contain a list of assignments of the type TAG = VALUE,
					 where TAG is a descriptor name (use the std::string labels as TAGs).
		*/
		virtual void import(ftool::md_param & fname) override;

		//! Checks validity of custom descriptors.
		bool elliptic_check_validity(void) const;

		//! Overridden print() method to print elliptic_problem specific descriptors.
		virtual void print(void) const override;
	};

}
#endif
