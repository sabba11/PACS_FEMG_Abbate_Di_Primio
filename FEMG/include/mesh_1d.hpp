/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015

       Copyright (C) 2015 Domenico Notaro, 2014 Laura Cattaneo
======================================================================*/
/*!
  @file   mesh1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   April 2015.
  @brief  Miscelleanous handlers for 1D network mesh.
  @details
  Create the edge sequence and build the related 1D mesh. \n
  The input stream ist is used to read a file with the following format (gambit-like): \n

		BEGIN_LIST
		BEGIN_ARC
		BC KEYWA [VALUES]
		BC KEYWB [VALUES]
		 106       0.4421      -1.6311       2.5101		start
		 107       0.4421      -1.6311       7.5101		end
		 108       0.3546      -1.6524       2.5539		point
		 109       0.2621      -1.6695       2.5880		point
		...
		END_ARC
		...
		BEGIN_ARC
		...
		END_ARC
		...
		END_LIST

  1. The list of points of the IS ordered as follows:
     - first we have the coordinates of TWO ENDS (A,B) (i.e. A=start and B=end)
     - then we have all the remaining nodes of the arc, from A to B
  2. If a node is shared by more arcs, the arcs will be CONNECTED in the resulting 1D mesh.
  3. BC KEYWA [VALUES] and BC KEYWB [VALUES] are keywords/values related to boundary conditions
     to be imposed at nodes A, B. Each KEYW [VALUES] entry can be one of the following: \n
     - BC DIR 1.1
     - BC MIX
     - BC INT \n
     Correspondingly, each node will be marked with the associated boundary condition, that are:
     - Dirichlet node (pressure = 1.1)
     - Mixed     node (flux = coef*(pressure - p0))
     - Internal  node
     At this stage, this is only meant to assign such BC labels to each node.
     If one end is INTERNAL, the corresponding BC will be ignored
     (for clarity, please write the INT keyword).
*/
/*!
	\defgroup geom Problem geometry
 */

#ifndef FEMG_MESH_1D_HPP_
#define FEMG_MESH_1D_HPP_

#include <unordered_set>
#include <node.hpp>
#include <cmath>

namespace getfem {

// bool operator==(const base_node & lhs, const base_node & rhs);
// bool operator<(const base_node & lhs, const base_node & rhs);


/*!
	Import the network points from the file of points (pts) and build the mesh.

	\ingroup geom
 */
//! \note It also build vessel mesh regions (#=0 for branch 0, #=1 for branch 1, ...).
template<typename VEC>
void
import_pts_file(
		std::istream & ist,
		getfem::mesh & mh1D,
		std::vector<getfem::node> &  BCList,
		size_type & n_original,
    //vector_type & mesh_step,
		VEC & Nn,
		const std::string & MESH_TYPE
		)
{

	size_type Nb = 0; // nb of branches
	size_type last_added = 0; // auxiliary counter
	Nn.resize(0); Nn.clear();
	mh1D.clear();

	ist.precision(16);
	ist.seekg(0); ist.clear();
	GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"),
		"This seems not to be a data file");

	size_type globalBoundaries = 0;

  std::vector<base_node> rpoints;
  std::vector<base_node> spoints;
  std::vector<size_type> BC_check;

	while (bgeot::read_until(ist, "BEGIN_ARC")) {

		Nb++;
		Nn.emplace_back(0);

		std::vector<base_node> lpoints;

		dal::dynamic_array<scalar_type> tmpv;
		std::string tmp, BCtype, value;
		bool thend = false;
		size_type bcflag = 0;
		size_type bcintI = 0, bcintF = 0;
		node BCA, BCB;


		// Read an arc from data file and write to lpoints
		while (!thend) {
			bgeot::get_token(ist, tmp, 1023);
			if (tmp.compare("END_ARC") == 0) {
				thend = true;
			}
			else if (ist.eof()) {
				GMM_ASSERT1(0, "Unexpected end of stream");
			}
			else if (tmp.compare("BC") == 0) {
				bcflag++;
				bgeot::get_token(ist, BCtype, 4);
				if (BCtype.compare("DIR") == 0) {
					bgeot::get_token(ist, value, 1023);
					if (bcflag == 1) {
						BCA.label = BCtype;
						BCA.value = stof(value);
						//BCA.ind = globalBoundaries;
						globalBoundaries++;
					}
					else if (bcflag == 2) {
						BCB.label = BCtype;
						BCB.value = stof(value);
						//BCB.ind = globalBoundaries;
						globalBoundaries++;
					}
					else
						GMM_ASSERT1(0, "More than 2 BC found on this arc!");
				}
				else if (BCtype.compare("MIX") == 0) {
					bgeot::get_token(ist, value, 1023);
					if (bcflag == 1) {
						BCA.label = BCtype;
						BCA.value = stof(value);
						//BCA.ind = globalBoundaries;
						globalBoundaries++;
					}
					else if (bcflag == 2) {
						BCB.label = BCtype;
						BCB.value = stof(value);
						//BCB.ind = globalBoundaries;
						globalBoundaries++;
					}
				}
				else if (BCtype.compare("INT") == 0) {
					if (bcflag == 1) {
						bcintI++;
						BCA.label = BCtype;
						//BCA.value = stof(value); //Error: no number to read
					}
					else if (bcflag == 2) {
						bcintF++;
						BCB.label = BCtype;
						//BCB.value = stof(value); //Error: no number to read
					}
					else
						GMM_ASSERT1(0, "More than 2 BC found on this arc!");
				}
				else
					GMM_ASSERT1(0, "Unknown Boundary condition");

			} /* end of "BC" case */
			else if (tmp.size() == 0) {
				GMM_ASSERT1(0, "Syntax error in file, at token '"
							 << tmp << "', pos=" << std::streamoff(ist.tellg()));
			}
			else { /* "point" case */
				Nn[Nb-1]++;
				int d = 0;
				while ( (isdigit(tmp[0]) != 0) || tmp[0] == '-' || tmp[0] == '+' || tmp[0] == '.'){
					tmpv[d++] = stof(tmp);
					bgeot::get_token(ist, tmp, 1023);
				}
                if (d != 4) GMM_ASSERT1(0, "Points must have 3 coordinates - number of coordinates:" << d);

				base_node tmpn(tmpv[1], tmpv[2]);
				//base_node tmpn(tmpv[1], tmpv[2], tmpv[3]);
				lpoints.push_back(tmpn);
				if (tmp.compare("END_ARC") == 0) { thend = true; Nn[Nb-1]--; }
			}

		} /* end of inner while */

		// Insert the arc into the 1D mesh and build a new region for the corresponding branch
		// Check validity of branch region
		GMM_ASSERT1(mh1D.has_region(Nb-1)==0, "Overload in meshv region assembling!");


    //storing real points, their successor and the vector storing info for then assigning the right boundary condition
    rpoints.push_back(lpoints[0]);
    rpoints.push_back(lpoints[1]);
    spoints.push_back(lpoints[2]);
    size_type index_last = lpoints.size()-1;
    spoints.push_back(lpoints[index_last]);

    if ((bcflag>0) && (bcintI==0)){
			BCList.push_back(BCA);
      BC_check.push_back(BCList.size());
		}
    else
      BC_check.push_back(0);


	  if ((bcflag>1) && (bcintF==0)){
			BCList.push_back(BCB);
      BC_check.push_back(BCList.size());
		}
    else
      BC_check.push_back(0);

		//adding to the mesh internal points and sub-arcs only between internal points
		for (size_type i=2; i<lpoints.size()-1; ++i ){

			std::vector<size_type> ind(2);
			ind[0] = mh1D.add_point(lpoints[i]);
			ind[1] = mh1D.add_point(lpoints[i+1]);
			size_type cv;
			cv = mh1D.add_convex(bgeot::simplex_geotrans(1,1), ind.begin());

			// Build branch regions
			mh1D.region(Nb-1).add(cv);

		} /* end of inner for */
		//mesh_step.push_back(0); // it will be the (Nb-1)-th element
		//for (unsigned k = 0; k < rpoints[last_added].size(); k++) {
		//	mesh_step[Nb-1] += (rpoints[last_added][k]-spoints[last_added][k])*(rpoints[last_added][k]-spoints[last_added][k]);
		//}
		//mesh_step[Nb-1] = std::sqrt(mesh_step[Nb-1]);
		//last_added += 2; // skip sub-arc of ending point
	} /* end of outer while (end of an branch) */
	std::set<base_node> realpoints(rpoints.begin(), rpoints.end());
	n_original = realpoints.size();

  //adding to the mesh points and sub-arcs that have a real point
  for (size_type i=0; i<rpoints.size(); ++i){

      std::vector<size_type> ind(2);
			ind[0] = mh1D.add_point(rpoints[i]);
			ind[1] = mh1D.add_point(spoints[i]);
			size_type cv;
			cv = mh1D.add_convex(bgeot::simplex_geotrans(1,1), ind.begin());
            mh1D.region(Nb).add(cv);

      //assigning boundary condition to the corresponding point
      if (BC_check[i]>0)
          BCList[BC_check[i]-1].idx = ind[0];

	} /*end of rpoints for */

} /* end of import_pts_file */

} // end of namespace
#endif
