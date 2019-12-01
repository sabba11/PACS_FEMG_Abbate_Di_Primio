#include "quantum_graph_problem.hpp"

namespace getfem {
	void quantum_graph_problem::import_pts_file(std::istream & ist, std::istream & rad, const bool & IMPORT_RADIUS) {
		n_branches = 0; // nb of branches
		vector_type Rdata; // data structure to read radii
		size_type nb_branches = 0; // check number of radii is correct
		n_vertices.resize(0); n_vertices.clear();
		meshg.clear();
		bool dim_check = true;

		//opening mesh data file
		ist.precision(16);
		ist.seekg(0); ist.clear();
		GMM_ASSERT1(bgeot::read_until(ist, "BEGIN_LIST"),
			"This seems not to be a data file");

		//opening and reading radii data file
		if (IMPORT_RADIUS) {
			rad.precision(16);
			rad.seekg(0); rad.clear();
			GMM_ASSERT1(bgeot::read_until(rad, "BEGIN_LIST"), "This seems not to be a data file");
			std::string line;
			bool thend = false;
			while (!thend){
				bgeot::get_token(rad, line, 1023);
				thend = (line=="END_LIST");
				if (!thend){
					Rdata.emplace_back(stof(line));
					nb_branches++;
				}
			}
		}

		size_type globalBoundaries = 0;

		std::vector<base_node> rpoints;
		std::vector<base_node> spoints;
		std::vector<size_type> BC_check;
		std::vector<vector_type> branch_tg;
		vector_type radii_aux;
		//Reading mesh data file
		while (bgeot::read_until(ist, "BEGIN_ARC")) {

			n_branches++;
			n_vertices.emplace_back(0);

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
							globalBoundaries++;
						}
						else if (bcflag == 2) {
							BCB.label = BCtype;
							BCB.value = stof(value);
							globalBoundaries++;
						}
						else
							GMM_ASSERT1(0, "More than 2 BC found on this arc!");
				}
				else if (BCtype.compare("NEU") == 0) {
						bgeot::get_token(ist, value, 1023);
						if (bcflag == 1) {
							BCA.label = BCtype;
							BCA.value = stof(value);
							globalBoundaries++;
						}
						else if (bcflag == 2) {
							BCB.label = BCtype;
							BCB.value = stof(value);
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
					n_vertices[n_branches-1]++;
					int d = 0;
					while ( (isdigit(tmp[0]) != 0) || tmp[0] == '-' || tmp[0] == '+' || tmp[0] == '.'){
						tmpv[d++] = stof(tmp);
						bgeot::get_token(ist, tmp, 1023);
					}
					if (dim_check){
						if (d == 4 || d == 3)
							dim_prob = d;
						else
							GMM_ASSERT1(0, "Points must have 3 or 2 coordinates - number of coordinates:" << d);
						dim_check = false;
					}
					base_node tmpn;
					if (d == dim_prob) {
						if (dim_prob == 3) {
							base_node dim_tmpn(tmpv[1], tmpv[2]);
							tmpn = dim_tmpn;
						}
						else if (dim_prob == 4) {
							base_node dim_tmpn(tmpv[1], tmpv[2], tmpv[3]);
							tmpn = dim_tmpn;
						}
					}
					else
						GMM_ASSERT1(0, "Points do not have costant number of coordinates - first point n. of coordinates:" << dim_prob << " - current point n. of coordinates" << d);

					//base_node tmpn(tmpv[1], tmpv[2], tmpv[3]);
					lpoints.push_back(tmpn);
					if (tmp.compare("END_ARC") == 0){
						thend = true;
						n_vertices[n_branches-1]--;
					}
				}
			} /* end of inner while */

			// Insert the arc into the 1D mesh and build a new region for the corresponding branch
			// Check validity of branch region
			GMM_ASSERT1(meshg.has_region(n_branches-1)==0, "Overload in meshv region assembling!");

			// Compute normalized tangent vector
			branch_tg.resize(n_branches);
			for (unsigned i = 0; i<dim_prob; i++){
				branch_tg[n_branches-1].push_back(lpoints[1][i]-lpoints[0][i]);
			}
			scalar_type norm = sqrt(branch_tg[n_branches-1][0]*branch_tg[n_branches-1][0]+branch_tg[n_branches-1][1]*branch_tg[n_branches-1][1]+branch_tg[n_branches-1][2]*branch_tg[n_branches-1][2]);
			for (unsigned i = 0; i<dim_prob; i++){
				branch_tg[n_branches-1][i] = branch_tg[n_branches-1][i]/norm;
			}

			//Storing real points, their successor and the vector storing info for then assigning the right boundary condition
			rpoints.push_back(lpoints[0]); // start
			rpoints.push_back(lpoints[1]); // end
			spoints.push_back(lpoints[2]); // first discretization node
			size_type index_last = lpoints.size()-1;
			spoints.push_back(lpoints[index_last]); // last discretization node

			//Storing radii relative to real points
			radii_aux.push_back(Rdata[n_branches-1]);
			radii_aux.push_back(Rdata[n_branches-1]);

			//Updating boundary condition list
			if ((bcflag>0) && (bcintI==0)){
				BCg.push_back(BCA);
				BC_check.push_back(BCg.size());
			}
			else
				BC_check.push_back(0);

			if ((bcflag>1) && (bcintF==0)){
				BCg.push_back(BCB);
				BC_check.push_back(BCg.size());
			}
			else
				BC_check.push_back(0);

			//adding to the mesh internal points and sub-arcs only between internal points
			std::vector<size_type> ind(2);
			size_type cv;

			for (size_type i=2; i<lpoints.size()-1; ++i ){
				ind[0] = meshg.add_point(lpoints[i]);
				ind[1] = meshg.add_point(lpoints[i+1]);

				cv = meshg.add_convex(bgeot::simplex_geotrans(1,1), ind.begin());

				// Build branch regions
				meshg.region(n_branches-1).add(cv);

			} /* end of inner for */

		} /* end of outer while (end of a branch) */

		GMM_ASSERT1(nb_branches == n_branches, "Number of given radii is not equal to the number of branches!");
		//Indexes of additional regions to store BCs
		size_type neumann_bc_region = n_branches + 1;
		size_type dirichlet_bc_region = n_branches + 2;

		//Updating tangent vectors and radii structures to count discretization nodes
		std::set<base_node> realpoints(rpoints.begin(), rpoints.end());
		n_origvert = realpoints.size();
		for (size_type b=0; b<n_branches; ++b)
			for (unsigned i = 0; i < meshg.region(b).nb_convex()+1; ++i)
				tg_vectors.push_back(branch_tg[b]);
		if (IMPORT_RADIUS)
			for (size_type b=0; b<n_branches; ++b)
				for (unsigned i = 0; i < meshg.region(b).nb_convex()+1; ++i)
					radii.push_back(Rdata[b]);

		//adding to the mesh points and sub-arcs that have a real point
		GMM_ASSERT1(meshg.has_region(n_branches)==0 & meshg.has_region(n_branches+1)==0 & meshg.has_region(n_branches+2)==0 , "Overload in meshg region assembling!");
		for (size_type i=0; i<rpoints.size(); ++i){
			std::vector<size_type> ind(2);
			ind[0] = meshg.add_point(rpoints[i]);
			ind[1] = meshg.add_point(spoints[i]);
			size_type cv;
			cv = meshg.add_convex(bgeot::simplex_geotrans(1,1), ind.begin());
			meshg.region(n_branches).add(cv);

			vector_type v_tg;
			for (unsigned j = 0; j<dim_prob; j++){
				v_tg.push_back(rpoints[i][j]-spoints[i][j]);
			}
			scalar_type norm = sqrt(v_tg[0]*v_tg[0]+v_tg[1]*v_tg[1]+v_tg[2]*v_tg[2]);
			for (unsigned j = 0; j<dim_prob; j++){
				v_tg[j] = v_tg[j]/norm;
			}
			tg_vectors.push_back(v_tg);
			radii.push_back(radii_aux[i]);
			//assigning boundary condition to the corresponding point
			if (BC_check[i]>0) {
				BCg[BC_check[i]-1].idx = ind[0];
				if (BCg[BC_check[i]-1].label == "DIR")
					meshg.region(dirichlet_bc_region).add(cv, 1);
				else if (BCg[BC_check[i]-1].label == "NEU")
					meshg.region(neumann_bc_region).add(cv, 1);
				else if (BCg[BC_check[i]-1].label == "INT")
					meshg.region(neumann_bc_region).add(cv, 1);
				else
					GMM_ASSERT1(false, "Invalid boundary condition. Valid boundary conditions are: DIR (Dirichlet), NEU (Neumann-Kirchhoff).");
			}
		} /*end of rpoints for */
		//checking that boundary conditions related to the same node are coherent
		check_boundary_conditions();
		return;
	}

	bool quantum_graph_problem::check_boundary_conditions(void) {
		for (unsigned i = 0; i < BCg.size(); i++)
			std::cout << BCg[i] << std::endl;
		for (unsigned i = 0; i < BCg.size()-1; i++)
			for (unsigned j = i + 1; j < BCg.size(); j++)
				if (BCg[i].idx == BCg[j].idx)
					GMM_ASSERT1(BCg[i].label == BCg[j].label & BCg[i].value == BCg[j].value, "Incoherent boundary conditions. All BC data of the same node should be equal.");
		return true;
	}
}
