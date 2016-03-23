#include "GST.h"
#include <fstream>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <iterator>
#include <queue>


#define DEBUG

dreal INF = std::numeric_limits<dreal>::infinity();

// Boost libraries
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/iterator/counting_iterator.hpp>

#if defined(__GNUC__)
#include <boost/filesystem.hpp>
#endif

#include <boost/heap/priority_queue.hpp>

using namespace std;


#ifdef LOG
#if defined( _MSC_VER ) 
	#include <direct.h>
	void mkDir(const char* path)
	{
		_mkdir(path);
	}

#else 
	#include <sys/stat.h>
	void mkDir(const char* path)
	{
		mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
#endif
#endif

// #undef VERBOSE


typedef std::pair<int, int> edge;

// implementation of header file

GST_instance::GST_instance() { }


GST_solver::GST_solver(GST_instance& instance) :I(instance), frac_OPT(INF), OPT(INF) { mkDir("./logs"); }




GST_instance::GST_instance(filename stp_file_path)
{
	/*
	reading input from file.
	*/
	read_status = read_stp(stp_file_path);

	if (read_status == FILE_MISSING)
	{
		 std::cout << "Input file is missing, current input path is set to " << INSTANCE_PATH << std::endl;
		 
		 exit(1);
	}
	
}


// define copy constructor




// IO

status_t GST_instance::read_stp(filename stp_file)
{
	
	// catch an exception if file doesn't exist
	ifstream stp_handler;
	stp_handler.open(stp_file, ios::in);
	if (!stp_handler.is_open()) return FILE_MISSING;

#if defined(__GNUC__)
	boost::filesystem::path p(stp_file.c_str());
	this->stp_file = p.stem().string();
#else
	std::string _filename = (stp_file.substr(stp_file.find_last_of("\\") + 1));
	this->stp_file = _filename.substr(0, _filename.find_last_of("."));
#endif
	

	string line;
	integer V, E;
	integer u, v;
	sreal w;

	// edges as pairs
	vector<edge> edges;
	vector<float> wt;

	int i = 0;
	int g = 0; // _index of NodeList


	while (getline(stp_handler, line))
	{
		boost::algorithm::trim(line);
		std::vector<std::string> words;
		// split by delimiter " " 
		boost::split(words, line, boost::is_any_of(" "), boost::token_compress_on);
		// process data from input

		if (words[0] == std::string("Nodes")) V = atoi(words[1].c_str());
		if (words[0] == std::string("Edges")) { E = atoi(words[1].c_str()); edges.resize(E); wt.resize(E); }

		if (words[0] == std::string("E"))
		{
			u = atoi(words[1].c_str()); v = atoi(words[2].c_str()); w = (sreal)atof(words[3].c_str());
			edges[i].first = u, edges[i].second = v;
			wt[i] = w;
			i++;
		}

		if (words[0] == std::string("T"))
		{

			for (size_t t = 1; t < words.size(); t++)
			{
				int node = atoi(words[t].c_str());
				nodes_of_group[g].insert(node);
				groups_of_node[node].insert(g);
				//terminal_nodes.insert(node);
			}
			g++;
		}
		


		// prize collecting addition
		// ...

		if (words[0] == std::string("EOF")) break;
	}

	// define graph
	// make edges, nodes are _indexed (0,1,...,V-1)
	for (size_t i = 0; i < edges.size(); i++)
	{
		WeightedEdge *e = new WeightedEdge(edges[i].first, edges[i].second, wt[i]);
		G.insert_edge(e);
	}

#ifdef GRAPH_VERBOSE
	GraphX::IO<WeightedGraph, WeightedEdge>::print(G);

	NodeList adj_vertices = G.neighbourhood(1);

	cout << "Neighbourhood of " << 1 << ": ";
	for (auto c : adj_vertices) std::cout << c << ' ';

	cout << endl;
	cout << "Terminals/Groups: " << endl;
	for (int i = 0; i < nodes_of_group.size(); i++)
	{
		for (auto x : nodes_of_group[i]) cout << x << " ";
		std::cout << std::endl;
	}
#endif
	

	return FILE_OK;
}

bool GST_instance::isInstanceConnected()
{
	return G.isConnected();
}

NodeList& GST_instance::get_nodes(_index g)
{
	return nodes_of_group[g];
}

NodeList& GST_instance::get_groups(_index v)
{
	return groups_of_node[v];
}

WeightedEdgeList GST_instance::get_edges()
{
	return G.edges();
}


NodeList GST_instance::get_neighbours(_index v)
{
	return G.neighbourhood(v);
}

GroupList& GST_instance::get_all_groups()
{
	return nodes_of_group;
}




integer GST_instance::size() const
{
	return G.V();
}

integer GST_instance::number_of_groups() const { return (integer)nodes_of_group.size(); }

integer GST_instance::number_of_edges() const { return G.E(); }

integer GST_instance::bound_on_group_size() const 
{
	
	auto m = std::max_element(
		nodes_of_group.begin(),
		nodes_of_group.end(),
		[](const std::pair<_index, std::set<_index> > x, const std::pair<_index, std::set<_index> > y)
		{
			 return (x.second.size() < y.second.size());
		}
		);

	return (integer)(m->second).size();
}

bool GST_instance::in_group(_index g, _index v)
{
	for (auto gr: get_groups(v)) if (gr == g) return true;
	return false;
}




/* *************  GST solver implementation ************************** */


bool GST_solver::check_feasibility(gst_tree& T) const
{
	/*
		Test the feasibility of solution for given tree T		
	*/

	std::map< _index, bool > _checked;
	int checked = 0;

	for (auto &v : T.nodes())
	{
		for (auto& g : I.get_groups(v))  if (!_checked[g]) { _checked[g] = true; checked++; }
	}

	return (checked == I.number_of_groups());
}

status_t GST_solver::solve_Steiner_approx()
{


	/*
		Steiner Tree problem solver

	*/
	SteinerTreeGraph* alggraph = new SteinerTreeGraph(I);
	dreal radiusOfMoat;
	while (alggraph->get_number_of_cc()>1)
	{
		//sample all tight edges
		radiusOfMoat = alggraph->sampleAllTightEdges();
		alggraph->blowUpEdges(radiusOfMoat);// updates a dual constraint value for all affected edges
		alggraph->mergeWithSampledEdges();
		alggraph->refreshSampler();
	}
	approximate_solution = alggraph;
	return FEAS;
}



status_t GST_solver::GRB_solve_netflow_GR(_index root_group, sol_t solve_type)
{
	/*
	MILP solver for unrooted version of problem.

	*/

	// status checks
	status_t solve_status;
	if (I.read_status == FILE_MISSING) return FILE_MISSING;

	
	// define immutables
	integer E = I.number_of_edges();
	integer V = I.size();
	integer k = I.number_of_groups();

	// edge list
	WeightedEdgeList edges = I.get_edges(); // edges of graph
	WeightedEdgeList ind_edges;
	


	//  typedefs
	typedef std::map<_index, std::map<_index, std::map<_index, GRBVar> > > GRB_flow_var_t;
	typedef std::map< WeightedEdge*, GRBVar > GRB_edge_var_t;


	// define lambda 
	auto group_node = [V](_index i) { return (V + i); };

	// LP variables
	GRB_edge_var_t var_edge;  // edge variable 
	GRB_edge_var_t var_ind_edge; // induced edges for root group
	GRB_flow_var_t var_flow; // flow variable 

	// solution for induced edges
	edge_var_t sol_ind_edge;


	// induced edges for root group
	for (auto& v : I.get_nodes(root_group))
	{
		_index g = group_node(root_group);
		WeightedEdge *e = new WeightedEdge(v, g, 0.0);
		ind_edges.insert(e);
	}


	// define complement of rooted group

	NodeList non_root_groups;
	for (_index i = 0; i < k; i++) { if (i != root_group) non_root_groups.insert(i); }
	

#ifndef DEBUG
	try { 	// test for Gurobi exceptions
#endif

		// define GRB model (default: minimization)
		GRBEnv   env = GRBEnv();
		GRBModel model = GRBModel(env);
		
		
#ifndef VERBOSE
		model.getEnv().set(GRB_IntParam_OutputFlag, 0); // set verbose off
#endif  

		// ILP or MILP solver?
		switch (solve_type)
		{
		case ILP:
			// initialize edge vars, induced edges	
			for (auto& e : edges)      var_edge[e] = model.addVar(0.0, 1.0, e->wt(), GRB_BINARY);
			for (auto& e : ind_edges)  var_ind_edge[e] = model.addVar(0.0, 1.0, e->wt(), GRB_BINARY);


			// initialize flow vars
			// *** define flow on edges from root group to its root from all other groups
			for (auto& e : ind_edges)
			{
				Node v = e->v(); Node g = e->w();
				for (auto& i : non_root_groups) var_flow[i][v][g] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
			}


			// *** define flow from non root groups to root group
			for (auto& i : non_root_groups)
			{
				// add flow on (undirected) edge
				for (auto& e : edges)
				{
					_index s = e->v(); _index t = e->w();
					var_flow[i][s][t] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					var_flow[i][t][s] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				}

				// define flow from group node g to all nodes in group	(if nodes belong to more groups)
				_index g = group_node(i);
				for (auto& v : I.get_nodes(i))
				{
					var_flow[i][g][v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
				}
			}
			break;

		
		case LP:
			// initialize edge vars, induced edges	
			for (auto& e : edges)  var_edge[e] = model.addVar(0.0, 1.0, e->wt(), GRB_CONTINUOUS);
					
			for (auto& e : ind_edges)	var_ind_edge[e] = model.addVar(0.0, 1.0, e->wt(), GRB_CONTINUOUS);
			


			// initialize flow vars
			// *** define flow on edges from root group to its root from all groups
			for (auto& e : ind_edges)
			{
				Node v = e->v(); Node g = e->w();
				for (auto& i : non_root_groups) var_flow[i][v][g] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
			}


			// define flow from non root groups to root group
			for (auto& i : non_root_groups)
			{
				// add flow per (undirected) edge
				for (auto& e : edges)
				{
					_index s = e->v(); _index t = e->w();
					var_flow[i][s][t] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
					var_flow[i][t][s] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);

				}

				// define flow from group node g to all nodes in group	(if nodes belong to more groups)
				_index g = group_node(i);
				for (auto v : I.get_nodes(i))
				{
					var_flow[i][g][v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);

				}
			}
			break;

		}

		// UPDATE MODEL
		model.update();

		// use expression object to represent lhs of constraints
		GRBLinExpr expr;

		// CONSTRAINTS
		// *** one edge chosen from root group to the root
		for (auto& e : ind_edges) 	expr += var_ind_edge[e];
		model.addConstr(expr == 1.0);

		// iterate per group
		for (auto& i : non_root_groups)
		{
			_index g = group_node(i); // _index of group node
			_index g_r = group_node(root_group); // _index of root group

			// *** capacity constraint
			for (auto& e : edges)  // for edges of graph
			{
				_index s = e->v(); _index t = e->w();
				model.addConstr(var_flow[i][s][t] + var_flow[i][t][s] <= var_edge[e]);
			}
			for (auto& e : ind_edges)  // for induced edges
			{
				Node s = e->v(); Node t = e->w();
				model.addConstr(var_flow[i][s][t] <= var_ind_edge[e]);
			}

			// condition: flow conservation		
			for (_index v = 0; v < V; v++)
			{
				expr = GRBLinExpr();
				// ***add flow from graph nodes		
				for (auto& t : I.get_neighbours(v))
				{
					expr += var_flow[i][v][t]; // if non-root group node send flow from v 
					expr -= var_flow[i][t][v]; // send flow back to v
				}
				// receive from group node
				if (I.in_group(i, v)) expr -= var_flow[i][g][v];
				// add flow to root
				if (I.in_group(root_group, v)) expr += var_flow[i][v][g_r]; 
							
				model.addConstr(expr == 0.0);
				
			}
			
			// add flow from group nodes
			expr = GRBLinExpr();
			for (auto& t : I.get_nodes(i)) expr += var_flow[i][g][t]; // add flow to v from graph vertices
			model.addConstr(expr == 1.0);


			// add flow from root group to the root			
			expr = GRBLinExpr();
			for (auto& t : I.get_nodes(root_group)) expr -= var_flow[i][t][g_r]; // add flow to v from graph vertices
			model.addConstr(expr == -1.0);

			

		} // end for i ...
				  

		// OPTIMIZATION
		model.optimize();
			


		// check Gurobi status
		int optimstatus = model.get(GRB_IntAttr_Status);


		// choose if OPT or fractional
		dreal      *obj_val = NULL;
		edge_var_t *sol_x   = NULL;
		flow_var_t *sol_f   = NULL;

		switch (solve_type)
		{
		
		case ILP:
			obj_val =   &(this->OPT);
			sol_x   =   &(this->x_OPT);
			sol_f   =   &(this->f_OPT);
			break;

		case LP:
			obj_val = &(this->frac_OPT);
			sol_x   = &(this->x);
			sol_f   = &(this->f);
			break;			

		}

		
		switch (optimstatus)
		{
		case GRB_OPTIMAL:
			solve_status = FEAS;
			
			sol_f->clear(); sol_x->clear();
			*obj_val = model.get(GRB_DoubleAttr_ObjVal);

			// store solution
			for (auto& i : non_root_groups)
			{
				// flow_vars
				for (auto& v : I.get_nodes(i))
				{
					_index g = group_node(i);
					(*sol_f)[i][v] = var_flow[i][g][v].get(GRB_DoubleAttr_X);
				}
			}
			// edge vars
			for (auto& e : I.get_edges()) (*sol_x)[e] = var_edge[e].get(GRB_DoubleAttr_X);
			for (auto& e : ind_edges) sol_ind_edge[e] = var_ind_edge[e].get(GRB_DoubleAttr_X);
			
			
			// induced edge vars

			if (solve_type == ILP)
			{
				for(auto& e : ind_edges) 
					opt_log << "  " << e->v() << " " << e->w() << ": " << sol_ind_edge[e] << std::endl;
			}
			
		
			if (solve_type == LP)
			{
				dreal treshold = 0.0;
				int cnt = 0; // count number of fractional values larger than alpha
				
				for (auto& e : ind_edges)
				{
						approx_log << "" << e->v() << " " << e->w() << ": " << sol_ind_edge[e] << std::endl;
						if (sol_ind_edge[e] > treshold) ++cnt;
				}

					
				std::cout << std::endl;
				std::cout << " -> checking only fractional values of nodes in root group ..." << std::endl;
				std::cout << " *** No. of nodes with large enough fractional values: " << cnt << std::endl;
				std::cout << " *** Optimal fractional value: " << *obj_val << std::endl;


				return NO_ROUND;
			}

			
			break;

		case GRB_INF_OR_UNBD:
			solve_status = INFEAS_OR_UNBND;
			break;
		case GRB_INFEASIBLE:
			solve_status = INFEAS;
			std::cout << "Model is infeasible. " << std::endl;
			break;
		}

#ifndef DEBUG
	}

	catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
		return UNDEF_ERROR;


	}
	catch (...)
	{
		cout << "Undefiend exception during optimization" << endl;
		return UNDEF_ERROR;
	}


	// if done, flag the status_t, use Gurobi flags
	

#endif
		

	return solve_status;
}

status_t GST_solver::GRB_solve_netflow_NR(_index r,sol_t solve_type)
{
	/* 
		ILP/LP compact model of GST problem.
	paramters:
		_index r          ... root of the Steiner Tree
		sol_t solve_type ... type of linear program (ILP,LP)

	return:
		returns status report of solver
	*/

	// status checks
	status_t solve_status;
	if (I.read_status == FILE_MISSING) return FILE_MISSING;

		
	// define immutables
	integer E = I.number_of_edges();
	integer V = I.size();
	integer k = I.number_of_groups();
	
	// edge list
	WeightedEdgeList edges = I.get_edges(); // edges of graph
	
	//  typedefs
	typedef std::map<_index,std::map<_index,std::map<_index,GRBVar> > > GRB_flow_var_t;
	typedef std::map< WeightedEdge*, GRBVar > GRB_edge_var_t;
	typedef std::map< WeightedEdge*, GRBVar > GRB_edge_var_t;
	
	// define lambda 
	auto group_node = [V](_index i) { return (V + i); };

	// LP variables
	GRB_edge_var_t var_edge;  // edge variable 
	GRB_flow_var_t var_flow; // flow variable 
	

	// discard groups with root
	NodeList root_groups = I.get_groups(r);
	non_root_groups.clear();
	NodeList all_groups(boost::counting_iterator<_index>(0), boost::counting_iterator<_index>(k));
	boost::set_difference(all_groups, root_groups, std::inserter(non_root_groups, non_root_groups.begin()));
	
	try { // test for Gurobi exceptions

		// define GRB model (default: minimization)
		GRBEnv   env = GRBEnv();
		GRBModel model = GRBModel(env);
#ifndef VERBOSE
		model.getEnv().set(GRB_IntParam_OutputFlag, 0); // set verbose off
#endif

		// initialize edge vars	
		switch (solve_type)
		{
		case ILP:
			for (auto& e : edges)  var_edge[e] = model.addVar(0.0, 1.0, e->wt(), GRB_BINARY);

			// initialize flow vars
			for (auto& i : non_root_groups)
			{
				// add flow per (undirected) edge
				for (auto& e : edges)
				{
					_index s = e->v(); _index t = e->w();
					var_flow[i][s][t] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
					var_flow[i][t][s] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				}


				// define flow from group node g to all nodes in group	(if nodes belong to more groups)	

				_index g = group_node(i);
				for (auto& v : I.get_nodes(i))
				{
					var_flow[i][g][v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

				}
			}

			break;
		case LP:
			for (auto& e : edges)  var_edge[e] = model.addVar(0.0, 1.0, e->wt(), GRB_CONTINUOUS);

			// initialize flow vars
			for (auto i : non_root_groups)
			{
				// add flow per (undirected) edge
				for (auto& e : edges)
				{
					_index s = e->v(); _index t = e->w();
					var_flow[i][s][t] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
					var_flow[i][t][s] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);

				}


				// define flow from group node g to all nodes in group	(if nodes belong to more groups)	

				_index g = group_node(i);
				for (auto& v : I.get_nodes(i))
				{
					var_flow[i][g][v] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);

				}
			}
			break;

		default:
			return UNDEF_ERROR;
			break;
		}



		// update the model
		model.update();

		// use expression object to represent lhs of constraints
		GRBLinExpr expr;

		// define constraints
		// iterate per group
		for (auto& i : non_root_groups)
		{
			_index g = group_node(i);
			// condition: capacity constraint
			for (auto& e : edges)
			{
				_index s = e->v(); _index t = e->w();
				model.addConstr(var_flow[i][s][t] + var_flow[i][t][s] <= var_edge[e]);
			}

			// condition: flow conservation		
			for (_index v = 0; v < V; v++)
			{
				expr = GRBLinExpr();
				// -> add flow from graph nodes		
				for (auto t : I.get_neighbours(v))
				{
					if (v != r) expr += var_flow[i][v][t]; // if non-root node send flow from v 
					if (t != r) expr -= var_flow[i][t][v]; // send flow back to v
				}

				if (v != r)
				{
					if (I.in_group(i, v)) expr -= var_flow[i][g][v];
				}

				if (v == r) model.addConstr(expr == -1.0);
				else model.addConstr(expr == 0.0);

			}

			// add flow from group nodes
			expr = GRBLinExpr();
			for (auto& t : I.get_nodes(i)) expr += var_flow[i][g][t]; // add flow to v from graph vertices
			model.addConstr(expr == 1.0);

		} // end for i ...

		// optimize the model
		model.optimize();

		// choose if OPT or fractional
		dreal      *obj_val  = NULL;
		edge_var_t *sol_x    = NULL;
		flow_var_t *sol_f    = NULL;

		switch (solve_type)
		{
		case LP:
			obj_val = &(this->frac_OPT);
			sol_x = &(this->x);
			sol_f = &(this->f);

			break;

		case ILP:

			obj_val = &(this->OPT);
			sol_x = &(this->x_OPT);
			sol_f = &(this->f_OPT);
		}
		// check Gurobi status
		int optimstatus = model.get(GRB_IntAttr_Status);
		
		switch (optimstatus)
		{
		case GRB_OPTIMAL:
			solve_status = FEAS;
			temp_obj_val = model.get(GRB_DoubleAttr_ObjVal);

			if (temp_obj_val < *obj_val)
			{
				*obj_val = temp_obj_val;

				sol_f->clear(); sol_x->clear();

				// store solution
				for (auto& i : non_root_groups)
				{
					// flow_vars
					for (auto& v : I.get_nodes(i))
					{
						_index g = group_node(i);
						(*sol_f)[i][v] = var_flow[i][g][v].get(GRB_DoubleAttr_X);
					}
				}
				// edge vars

				for (auto& e : I.get_edges()) (*sol_x)[e] = var_edge[e].get(GRB_DoubleAttr_X);


				// set root
				root = r;
			}


			break;

		case GRB_INF_OR_UNBD:
			solve_status = INFEAS_OR_UNBND;
			break;
		case GRB_INFEASIBLE:
			solve_status = INFEAS;
			break;

		}
}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << endl;
		std::cout << e.getMessage() << endl;

	}
	catch (...) {
		std::cout << "Exception during optimization" << endl;
		return UNDEF_ERROR;
	}

		// if done, flag the status_t, use Gurobi flags
	return solve_status;
}


status_t GST_solver::solve_opt(model_t model)
{
	/*
	find an optimal solution to the GST problem instance.
	*/



	
#ifdef USE_TIMER
	Timer opt_time;
	opt_time.start();
#endif
	status_t GRB_status, ROUND_status;

# ifdef LOG
	std::ostringstream stringStream, logPath;
	stringStream << "opt_log_" << I.stp_file << ((model == type1) ? "_NR" : "_GR") << "_.log";
	std::cout << "Opening log file: " << stringStream.str() << std::endl;
	logPath << LOG_PATH << stringStream.str();
	
	opt_log.open(logPath.str().c_str());

	opt_log << "*** OPTIMAL SOLUTION ***" << std::endl;	
	opt_log << "** INSTANCE:" << std::endl;
	opt_log << "V " << I.G.V() << std::endl;
	opt_log << "E " << I.G.E() << std::endl;
	opt_log << "K " << I.number_of_groups() << std::endl;
	opt_log << std::endl;
	opt_log << "** LOG:" << std::endl;
		
#endif
	
	switch (model){
	case  type1:
		// finding root candidates
		std::cout << "\n------------------------------------------------------------------------" << std::endl;
		std::cout << "      Optimal algorithm based on node root (NR) ILP flow based formulation    " << std::endl;
		std::cout << "---------------------------------------------------------------------------" << std::endl;
		root_group = find_root_group();
		root_candidates = I.nodes_of_group[root_group];
		std::cout << " -> Solving ILP's with respect to root candidates: ";
		// logging
		opt_log << "ROOTS/OBJ:" << std::endl;

		// solving ILP model

		
		for (auto& r : root_candidates)
		{
			
			std::cout << r << " ";
			opt_log << "" << r << ": ";
			GRB_status = GRB_solve_netflow_NR(r, ILP);
						
			opt_log << temp_obj_val << std::endl;
		}
		
		
		std::cout << "\n\n -> Optimal solution achieved with root: " << root << std::endl;
		// logging
		opt_log << "\nOPT ROOT:      " << root << std::endl; 
		opt_log << "SOLUTION TYPE: " << "INTEGRAL" << std::endl;
		

		break;

	case type2:
		std::cout << "\n------------------------------------------------------------------------" << std::endl;
		std::cout << "      Optimal algorithm based on group root (GR) LP flow-based formulation       " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
		root_group = find_root_group();
		std::cout << " -> Solving ILP's with respect to root group: " << root_group << std::endl;

		// logging 
		opt_log << "ROOT GROUP: " << root_group << std::endl;
		opt_log << "INDUCED EDGES:" << std::endl;
		GRB_status = GRB_solve_netflow_GR(root_group, ILP);
		break;
	}

	dreal opt_time_secs = 0;

#ifdef USE_TIMER
	opt_time.stop();
	opt_time_secs = opt_time.secs();

#endif
	gst_tree T;

	switch (GRB_status)
	{
	case FEAS:
		// stdout
		std::cout << " \n " << std::endl;
		std::cout << "*** Optimal value   :				" << OPT           << "		   " << std::endl;
		std::cout << "*** Total time      :				" << opt_time_secs << "s      		   " << std::endl;
		std::cout << "\n\n";
		// logging
		opt_log << "OPTIMAL VALUE: " << OPT << std::endl;
		opt_log << "TOTAL TIME:    " << opt_time.secs() << std::endl;

	
		for (auto& e : x_OPT)
		{
			if (e.second != 0)
			{
				EdgeWithConnComp *ec = new EdgeWithConnComp(e.first);
				T.insert_edge(ec);
			}
		}

		opt_log << "** SOLUTION:" << (check_feasibility(T) ? " feasible" : " non feasible") << std::endl;
		for (auto& e : x_OPT) 
		{
			if (e.second != 0) opt_log << e.first->v() << "--" << e.first->w() << ":" << e.first->wt() << std::endl;
		}
		break;
	case INFEAS:
		std::cout << "Model is infeasible. " << std::endl;
		break;
	default:
		std::cout << "Error in model. Exiting." << std::endl;
		break;
	}

	opt_log << flush;
	opt_log.close();
	return GRB_status;

}


status_t GST_solver::solve_approx(model_t model)
{
	/*

	solving approximation of GST problem based on the specified model: type1, type2

	*/
		

	// finding root candidates
#ifdef USE_TIMER
	Timer approx_time;
	approx_time.start();

	Timer Gurobi_time, define_terminals_time, steiner_solver_time;
	
	status_t GRB_status, ROUND_status;

#endif

# ifdef LOG
	std::ostringstream stringStream, logPath;
	stringStream << "approx_log_" << I.stp_file << ((model == type1) ? "_NR" : "_GR") << "_.log";
	std::cout << "Opening log file: " << stringStream.str() << std::endl;
	logPath << LOG_PATH << stringStream.str();

	approx_log.open(logPath.str().c_str());

	approx_log << "*** APPROXIMATION ***" << std::endl;
	approx_log << "INSTANCE:" << std::endl;
	approx_log << "V " << I.G.V() << std::endl;
	approx_log << "E " << I.G.E() << std::endl;
	approx_log << "K " << I.number_of_groups() << std::endl;
	approx_log << std::endl;
	approx_log << "** LOG:" << std::endl;
#endif


	status_t POST_ROUND_status = FEAS;

	switch (model)
	{
	case  type1:
		// finding root candidates
		std::cout << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "          Approximation based on node root LP net flow formulation:         "	<< std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		root_group = find_root_group();
		root_candidates = I.nodes_of_group[root_group];

		
		// logging
		approx_log << "ROOTS/OBJ:" << std::endl;
		
		std::cout << " -> Solving LP with respect to root candidates: ";
		
		// solving LP model
		for (auto &r : root_candidates)
		{
			std::cout << r << " ";
			Gurobi_time.start();
			GRB_status = GRB_solve_netflow_NR(r, LP);
			Gurobi_time.stop();

			// logging
			approx_log << "" << r << ": " << temp_obj_val << std::endl;
		}


		
		
		std::cout << "\n\n -> Best solution achieved with root: " << root << std::endl;

		// logging
		approx_log << "\nOPT ROOT:" << root << std::endl; 
		approx_log << "OBJ:     " << frac_OPT << std::endl;
		
		break;

	case type2:
		std::cout << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "   Approximation algorithm based on group root LP flow-based formulation?   " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		root_group = find_root_group();
		std::cout << " -> Solving LP with respect to root group: " << root_group << std::endl;
		Gurobi_time.start();
		approx_log << "INDUCED EDGES:" << std::endl;
		GRB_status = GRB_solve_netflow_GR(root_group, LP);
		Gurobi_time.stop();

		// logging 
		approx_log << "ROOT GROUP: " << root_group << std::endl;
		approx_log << "OBJ:     " << frac_OPT << std::endl;
		
		break;
				
		
	}

	if (model == type1)
	{
		// updating root groups
		NodeList root_groups = I.get_groups(root);
		non_root_groups.clear();
		NodeList all_groups(boost::counting_iterator<_index>(0), boost::counting_iterator<_index>(I.number_of_groups()));
		boost::set_difference(all_groups, root_groups, std::inserter(non_root_groups, non_root_groups.begin()));
	}
	

	status_t integrality_status = check_integrality(&x, &f);
	approx_log << "SOLUTION TYPE: " << ((integrality_status == INTEGRAL) ? "INTEGRAL" : "FRACTIONAL") << std::endl;

	dreal approx_time_secs = 0;
	// cost of approximate solution
	dreal wt_sol = 0.0;

	if (integrality_status == FRACTIONAL)
	{
		// if solution is fractional, proceed to solving steiner tree problem.
		if (GRB_status != NO_ROUND)
		{
			
			define_terminals_time.start();
			define_terminals();
			define_terminals_time.stop();

			std::cout << " -> Rounding solution of LP using the best root: " << root << std::endl;

			steiner_solver_time.start();
			solve_Steiner_approx();
			deleteSteinerLeaves();
			wt_sol = weightApprox();
			steiner_solver_time.stop();

			if (check_feasibility(*approximate_solution)) POST_ROUND_status = FEAS; else POST_ROUND_status = INFEAS;
			// update integral solution = approx
			x.clear();
			for (auto& e : approximate_solution->edges()) x[e] = 1;
		}

		
	}
	else
	{	
		wt_sol = frac_OPT;
			
	}

#ifdef USE_TIMER
	approx_time.stop();
	approx_time_secs = approx_time.secs();
#endif



	std::cout << "\n                                                                         " << std::endl;
	std::cout << "*** Fractional value:				" << frac_OPT << "		                   " << std::endl;
	std::cout << "*** Approx value    :				" << wt_sol << "		                   " << std::endl;
	std::cout << "*** Total time      :				" << approx_time_secs << "s     		   " << std::endl;
	std::cout << "-----------------------------------------------------------------------------" << std::endl;
	std::cout << "\n\n";

	// logging
	approx_log << "APPROX:           " << wt_sol << std::endl;

	approx_log << "\n** TIME:" << std::endl;
	if(model==type1) 
	approx_log << "AVG. TIME for LPs:       " << Gurobi_time.secs() / root_candidates.size() << std::endl;
	else 
	approx_log << "TIME for LPs:            " << Gurobi_time.secs() << std::endl;
	approx_log << "DEFINE TERMINALS:        " << define_terminals_time.secs() << std::endl;
	approx_log << "STEINER SOLVER  :        " << steiner_solver_time.secs() << std::endl;
	approx_log << "TOTAL TIME:              " << approx_time_secs << std::endl;
	approx_log << "** SOLUTION:" << (POST_ROUND_status != INFEAS ? " feasible " : " non feasible") << std::endl;
	for (auto& e : x) if (e.second != 0) approx_log << e.first->v() << "--" << e.first->w() << ":" << e.first->wt() << std::endl;
	
	approx_log << flush;
	approx_log.close();
	return FEAS;
}



_index GST_solver::find_root_group()
{
	/*
		find smallest group which is a root group candidate or contains best node root
	*/
	typedef std::pair<_index, std::set<_index> > dict_iter;
	// find the smalllest group of nodes
	
	auto min = std::min_element
		(
			I.nodes_of_group.begin(), 
			I.nodes_of_group.end(),
			[](const dict_iter& x, const dict_iter& y){ return (x.second.size() < y.second.size()); }
		);

	_index root_group = min->first;

	return root_group;
}




status_t GST_solver::define_terminals()
{
	/*
	method represents group as terminal for large enough fractional flow value. In case, of tie ...
	*/
	
	// set of groups incident to node j
	auto voting_groups = [this](_index j)
	{
		// set of groups incident to node j
		std::set<_index> G_j;
		for (auto& i : I.get_groups(j))
		{
			if (f[i][j] >= 1.0 / I.get_nodes(i).size()) G_j.insert(i);
		}
		
		return G_j;
	};

		
	// voting function
	auto vote = [this,voting_groups](_index j)
	{
		// vote value for some node j

		auto _vote_group = voting_groups(j);
		dreal max_f = -1;
		_index i_max = *_vote_group.begin();

		for (auto& g : _vote_group)
			if (max_f < f[g][j]) { max_f = f[g][j]; i_max = g; }

		
		//find max
		
		auto val = _vote_group.size()>0 ? f[i_max][j] / _vote_group.size() : 0.0;
		
		return pair<_index,dreal>(i_max,val);
				
	};

	
	// define priority queueue for finding nodes to round
	
	auto vote_Comp = [vote](_index i, _index j)
	{
		return (vote(i).second < vote(j).second);
	};
	
	// find terminals in groups
	_index terminal;
	
	GroupList nodes_of_group;
	GroupList groups_of_node;

	// use priority queue for selecting terminals
	boost::heap::priority_queue<_index, boost::heap::compare< decltype(vote_Comp) > > _terminal_candidates(vote_Comp);

	// find terminals in non root groups
	for (auto& i : non_root_groups)
	{
		_terminal_candidates.clear();

		auto& _nodes = I.get_nodes(i);

		for (auto& j : _nodes)
		{			
			if (f[i][j] >= 1.0 / _nodes.size()) _terminal_candidates.push(j);
		}

		// extract max from priority queue
		terminal = _terminal_candidates.top();
		_terminal_candidates.pop();
		
		// update instance with group which best votes
		_index vote_group = vote(terminal).first;
				
		if (nodes_of_group[vote_group].size() == 0)
		{
			nodes_of_group[vote_group].insert(terminal);
			groups_of_node[terminal].insert(vote_group);
		}
		
	}

	nodes_of_group[root_group].insert(root);
	groups_of_node[root].insert(root_group);
	

	// copy values from tentative instance to original instance
	I.nodes_of_group = nodes_of_group;
	I.groups_of_node = groups_of_node;
	
	// add root to instance
	return UNDEF_ERROR;
}


status_t GST_solver::check_integrality(edge_var_t *x, flow_var_t *f) const
{
	dreal val, tol = 1e-5;
	for (auto e : *x)
	{
		val = e.second;
		if (fabs(val-round(val)) > tol) return FRACTIONAL;
		
	}
	return INTEGRAL;
}



EdgeWithConnComp::EdgeWithConnComp(int v, int w, double wt) : WeightedEdge(v, w, wt), dualConstrValue(0)
{
	firstCut = NULL;
	secondCut = NULL;
	indOfFirstCut = -1;
	indOfSecondCut = -1;

}

EdgeWithConnComp::EdgeWithConnComp(WeightedEdge* e) : WeightedEdge(*e){}

void EdgeWithConnComp::increaseDualConstrValue(dreal delta)
{
	dualConstrValue += delta;
}

dreal EdgeWithConnComp::calculatePrice()
{
	return (wght - dualConstrValue) / numOfDeltaSets; // may be +inf, but it is OK!!!
}

dreal EdgeWithConnComp::get_price() const
{
	return price;
}

void EdgeWithConnComp::set_price(dreal value)
{
	price = value;
}

void EdgeWithConnComp::recalculateNumOfDSets()
{
	numOfDeltaSets = (_index)(firstCut != NULL) + (_index)(secondCut != NULL);
	return;
}


_index EdgeWithConnComp::get_numOfDSets() const
{
	return numOfDeltaSets;
}

void EdgeWithConnComp::addConnectedComponent(ConnectedComponent<EdgeWithConnComp>* cc, _index ind)
{
	if (firstCut != NULL && secondCut == NULL) { secondCut = cc; indOfSecondCut = ind; return; }
	if (firstCut == NULL && secondCut != NULL) { firstCut = cc; indOfFirstCut = ind; return; }
	if (firstCut == NULL && secondCut == NULL){ firstCut = cc; indOfFirstCut = ind; return; }
	if (firstCut != NULL && secondCut != NULL){ std::cout << "Warning: Edge can be incident to at most 2 connected components." << endl; return; }

}


ConnectedComponent<EdgeWithConnComp>* EdgeWithConnComp::getConnectedComponent(_index i) const
{
	if (i == 1){ return firstCut; }
	if (i == 2){ return secondCut; }
	std::cout << "\nError: There is no connected component with idex i " << i << " only 1 or 2" << endl;
	return NULL;
}

_index EdgeWithConnComp::getIndOfConnComp(_index i) const
{
	if (i == 1){ return indOfFirstCut; }
	if (i == 2){ return indOfSecondCut; }
	std::cout << "\nError: There is no connected component with idex i " << i << " only 1 or 2" << endl;
	return -1;
}


bool EdgeWithConnComp::operator<(const EdgeWithConnComp& e) const
{
	return this->price < e.get_price();
}

bool EdgeWithConnComp::isOnCut(ConnectedComponent<EdgeWithConnComp>* cc)
{
	if (cc == firstCut || cc == secondCut) return true;
	return false;
}


EdgeWithConnComp::~EdgeWithConnComp() {}


SteinerTreeGraph::SteinerTreeGraph(GST_instance& SteinerTreeInstance) : UndirectedGraph<EdgeWithConnComp>()
{
	// get all terminals that define connected components at the start of the algorithm
	GroupList groups = SteinerTreeInstance.get_all_groups();
	numberOfConnComp = SteinerTreeInstance.number_of_groups();
	for (auto itg : groups) connComp[*itg.second.begin()] = new ConnectedComponent<EdgeWithConnComp>(itg.second);

	// create a priority queue for all edges
	edgeSampler = new PriorityQueue<EdgeWithConnComp>(SteinerTreeInstance.number_of_edges());
	std::set<WeightedEdge*> edge_list = SteinerTreeInstance.get_edges();
	for (auto ite : edge_list)
	{
		integer start = ite->v(), end = ite->w();
		EdgeWithConnComp* algedge = new EdgeWithConnComp(start, end, ite->wt());
		if (connComp.count(start))
		{
			algedge->addConnectedComponent(connComp[start], start);
			connCompToEdges[connComp[start]].insert(algedge);
		}
		else
		{
			nodeToEdges[start].insert(algedge);
		}
		if (connComp.count(end))
		{
			algedge->addConnectedComponent(connComp[end], end);
			connCompToEdges[connComp[end]].insert(algedge);
		}
		else
		{
			nodeToEdges[end].insert(algedge);
		}
		algedge->recalculateNumOfDSets();
		algedge->set_price(algedge->calculatePrice());
		edgeSampler->push(algedge);
	}
}


EdgeWithConnComp* SteinerTreeGraph::sampleEdge()
{
	EdgeWithConnComp* sampled_edge = edgeSampler->top();
	edgeSampler->pop();
	this->insert_edge(sampled_edge);
	return sampled_edge;
}

bool SteinerTreeGraph::isConnected()
{
	return connComp.size() == 1;
}

dreal SteinerTreeGraph::sampleAllTightEdges()
{
	dreal price, priceNext;
	countEdges = 0;
	do
	{
		EdgeWithConnComp* firstTight = this->edgeSampler->top();
		edgeSampler->pop();
		orderOfEdges.insert(firstTight);
		this->insert_edge(firstTight);
		price = firstTight->get_price();
		countEdges++;
		if (edgeSampler->isEmpty()) break;
		priceNext = edgeSampler->top()->get_price();
	} while (price == priceNext);
	return price;
}

void SteinerTreeGraph::mergeComponents(EdgeWithConnComp* addedEdge)
{
	ConnectedComponent<EdgeWithConnComp>* firstConn = addedEdge->firstCut;
	ConnectedComponent<EdgeWithConnComp>* secondConn = addedEdge->secondCut;
	_index indOfFirstCut = addedEdge->indOfFirstCut;
	_index indOfSecondCut = addedEdge->indOfSecondCut;
	if ((firstConn != NULL && secondConn == NULL) || (firstConn == NULL && secondConn != NULL))
	{
		ConnectedComponent<EdgeWithConnComp>* comp;
		Node newnode;
		if (firstConn != NULL) comp = firstConn;
		else comp = secondConn;

		// first component is augmented
		if (comp->isElement(addedEdge->v())) newnode = addedEdge->w();
		else newnode = addedEdge->v();

		for (auto e : nodeToEdges[newnode])
		{
			if (e->isOnCut(comp))// if edge is not on the cut of the new component it is not tight or if it is, it becomes unnecesary
				//remove it from priority queue
			{
				edgeSampler->deleteElement(e);
				connCompToEdges[comp].erase(e);
				e->firstCut = NULL;
				e->secondCut = NULL;
				//this->remove_edge(e);
			}
			else //if edge is on the cut of the old component
			{
				e->addConnectedComponent(comp, indOfFirstCut);
				connCompToEdges[comp].insert(e);
			}

		}
		comp->add_node(newnode);
		nodeToEdges.erase(newnode);
		// addedEdge is not anymore on the cut
		//connCompToEdges[comp].erase(addedEdge);
		return;
	}
	if (firstConn != NULL && secondConn != NULL) // it always adds second cut to the first
	{
		for (auto e : connCompToEdges[secondConn]) // all edges, that are on the cut defined by secondConn component, transfer to the
			// cut defined by the firstConn
			// delete all edges that are in both components. These edges will not be on the cut defined by the new merged component
		{
			if (e->isOnCut(firstConn))
			{
				connCompToEdges[firstConn].erase(e);
				e->firstCut = e->secondCut = NULL;
				e->indOfFirstCut = e->indOfSecondCut = -1;
				edgeSampler->deleteElement(e);
				//this->remove_edge(e);
				continue;
			}
			if (e->firstCut == secondConn)
			{
				e->firstCut = firstConn;
				e->indOfFirstCut = indOfFirstCut;
				connCompToEdges[firstConn].insert(e);

			}
			if (e->secondCut == secondConn)
			{
				e->secondCut = firstConn;
				e->indOfSecondCut = indOfFirstCut;
				connCompToEdges[firstConn].insert(e);
			}
		}
		//connCompToEdges[firstConn].insert(connCompToEdges[secondConn].begin(), connCompToEdges[secondConn].end()); // add edges that are on the cut defined by the secondConn
		firstConn->add_component(secondConn); // add nodes from second component to th node set of the first
		connComp.erase(indOfSecondCut);
		connCompToEdges.erase(secondConn);
		numberOfConnComp--;
		return;
	}

	if (firstConn == NULL && secondConn == NULL)
	{
		this->remove_edge(addedEdge);
		return;
	}


}

void SteinerTreeGraph::mergeWithSampledEdges()
{
	for (auto e : orderOfEdges)
	{
		mergeComponents(e);
	}
	orderOfEdges.erase(orderOfEdges.begin(), orderOfEdges.end());
}

void SteinerTreeGraph::blowUpEdges(dreal delta)
{
	for (auto cc : connComp)
	{
		for (auto e : connCompToEdges[cc.second])
		{
			e->increaseDualConstrValue(delta);
			//e->recalculateNumOfDSets();
			//edgeSampler->changeKey(e, e->calculatePrice());
		}
	} // THE COMPONENT THAT IS MERGED SHOUD BE DELETED FROM connComp
}

void SteinerTreeGraph::refreshSampler()
{
	for (auto cc : connComp)
	{
		for (auto e : connCompToEdges[cc.second])
		{
			//e->increaseDualConstrValue(delta);
			e->recalculateNumOfDSets();
			edgeSampler->changeKey(e, e->calculatePrice());
		}
	} // THE COMPONENT THAT IS MERGED SHOUD BE DELETED FROM connComp
}

integer SteinerTreeGraph::get_number_of_cc() const
{
	return numberOfConnComp;
}


void GST_solver::deleteSteinerLeaves()
{
	//std::map<_index, _index> degrees;
	/*for (auto e : orderOfEdges)
	{
	Node start = e->v();
	Node end = e->w();
	if (degrees.count(start)) degrees[start]++;
	else degrees[start] = 1;
	if (degrees.count(end)) degrees[end]++;
	else degrees[end] = 1;
	}*/
	bool goFurther = true;
	NodeList listOfNodes = approximate_solution->nodes();
	while (goFurther)
	{
		goFurther = false;
		for (auto n : listOfNodes)
		{
			if (I.get_groups(n).size() == 0 && approximate_solution->degree_of_node(n) == 1)
			{
				approximate_solution->remove_node(n);
				listOfNodes = approximate_solution->nodes();
				goFurther = true;
				break;
			}
		}
	}
}

dreal GST_solver::weightApprox()
{
	dreal sum = 0;
	std::set<EdgeWithConnComp*> edges = approximate_solution->edges();
	for (auto e : edges)
	{
		sum += e->wt();
	}
	return sum;
}

_index ant_random_sample(vector<dreal>& pheromones, map<Node,dreal>& heuristicinf, dreal alpha, dreal beta, NodeList& nghSol)
{
	
	map<_index,dreal> probab;
	dreal sumpow = 0;
	for (auto& v : nghSol) { probab[v] = pow(pheromones[v], alpha)*pow(heuristicinf[v], beta); sumpow += probab[v]; }
	for (auto& v : nghSol) { probab[v] /= sumpow; /*assert(probab[v] <= 1 && probab[v] >= 0);*/ }
	dreal unif = ((dreal)rand()) / RAND_MAX;
	dreal cumul = 0;
	_index sel = *(nghSol.rbegin());
	for (auto& v : nghSol)
	{
		cumul += probab[v];
		if (unif <= cumul) { sel = v; break; }
	}
	/*if (sel > 60)
	{
		for (auto& v : nghSol) cout << v <<" "<<probab[v] << endl;
		cout << cumul << endl;
		cout << unif << endl;
		cout << sel<<endl;
	}*/
	return sel;
}

NodeList remove_vertex(GroupList& nodes_of_group, GroupList& groups_of_node, Node v)
{
	NodeList SteinerPoints;
	if (groups_of_node.count(v) == 0) return SteinerPoints;
	for (auto& x : groups_of_node[v])
	{
		for (auto& y : nodes_of_group[x])
		{
			if (y != v) groups_of_node[y].erase(x);
			if (groups_of_node[y].size() == 0) { groups_of_node.erase(y); SteinerPoints.insert(y); }
		}
	}
	for (auto& x : groups_of_node[v]) nodes_of_group.erase(x);
	groups_of_node.erase(v);
	return SteinerPoints;
}

sreal cost_expanding(DistanceMatrix& dist, PredAccesMatrix& prlist, GroupList& nodes_of_group, GroupList& groups_of_node, map<Node,sreal>& distToTree, map<Node,Node>& brotherInTree, WeightedGraph& G, Node x)
{
	sreal hinf=0;
	
	GroupList nOfG = nodes_of_group;
	GroupList GOfN = groups_of_node;
	map<Node, sreal> dTT = distToTree;
	map<Node, Node> brT = brotherInTree;

	NodeList SteinerP = remove_vertex(nOfG, GOfN, x);
	//dTT.erase(x);
	//brT.erase(x);
	//for (auto& y : SteinerP) { dTT.erase(y); brT.erase(y); }

	size_t n = GOfN.size(); // O(1)
	if (n == 0) return hinf;
	PriorityQueue<Pair<_index>> Q(n); // veli�ina prioritetnog reda je broj terminala nakon izbacivanja x-a
	map<_index,Pair<_index>*> pairs;

	for (auto& y : GOfN) { if (dist[x][y.first] < dTT[y.first]) { dTT[y.first] = dist[x][y.first]; brT[y.first] = x; } }
	for (auto& y: GOfN) { pairs[y.first] = (new Pair<_index>(y.first, dTT[y.first])); Q.push(pairs[y.first]); }
	while (nOfG.size() > 0)
	{
		Pair<_index>* next = Q.top();
		Q.pop();
		_index target = brT[next->v];
		Node pred = prlist[target][next->v];
		
		NodeList SteinerP = remove_vertex(nOfG, GOfN, next->v);
		//dTT.erase(next->v);
		//brT.erase(next->v);
		for (auto& y : SteinerP) { /*dTT.erase(y); brT.erase(y);*/ Q.deleteElement(pairs[y]); }
		for (auto& x : GOfN) { if (dTT[x.first] > dist[next->v][x.first]) { dTT[x.first] = dist[next->v][x.first]; Q.changeKey(pairs[x.first], dTT[x.first]); brT[x.first] = next->v; } }
		hinf += G.edge(pred, next->v)->wt();
		while (pred != target)
		{
			NodeList SteinerP = remove_vertex(nOfG, GOfN, pred);
			Q.deleteElement(pairs[pred]);
			//dTT.erase(pred);
			//brT.erase(pred);
			for (auto& y : SteinerP) {/*dTT.erase(y); brT.erase(y);*/ Q.deleteElement(pairs[y]); }
			for (auto& x : GOfN) { if (dTT[x.first] > dist[pred][x.first]) { dTT[x.first] = dist[pred][x.first]; Q.changeKey(pairs[x.first], dTT[x.first]); brT[x.first] = pred; } }
			hinf += G.edge(prlist[target][pred], pred)->wt();
			pred = prlist[target][pred];
		}
	}
	return hinf;
}

void lsa(WeightedGraph& AntSol, NodeList& Sol, WeightedGraph& G, sreal& weight)
{
	NodeList nghSol = G.neighbourhood(Sol); // optimizirati
	sreal improvedWeight = weight;
	WeightedGraph minMST;
	WeightedGraph inducedSub;
	
	for (auto& v : Sol)
	{
		WeightedGraph::adjIterator iterator(G, v);
		for (WeightedEdge* e = iterator.begin(); e != NULL; e = iterator.next())
		{
			if (Sol.count(e->w()) && Sol.count(e->v())) inducedSub.insert_edge(e);
		}
	}
	
	for (auto& x : nghSol)
	{
		WeightedGraph impMST;
		inducedSub.insert_node(x);
		NodeList xneig = G.neighbourhood(x);
		for (auto& y : xneig) if (Sol.count(y)) inducedSub.insert_edge(G.edge(x, y));
		improvedWeight = inducedSub.mst_prim(impMST);
		if (improvedWeight < weight)
		{
			weight = improvedWeight;
			AntSol = impMST;
		}
		inducedSub.remove_node(x);
	}
	return;
}

bool isCover(GroupList& groups_of_node, NodeList& Sol, _index k)
{
	NodeList coveredGroups;
	for (auto& x : Sol)
	{
		if (groups_of_node.count(x) == 0) continue;
		coveredGroups.insert(groups_of_node[x].begin(), groups_of_node[x].end());
		if (coveredGroups.size() == k) return true;
	}
	return false;
}

void lsr(WeightedGraph& AntSol, NodeList& Sol, GroupList& groups_of_node, WeightedGraph& G, sreal& weight, _index k)
{
	//NodeList nghSol = G.neighbourhood(Sol); // optimizirati
	sreal improvedWeight = weight;
	WeightedGraph inducedSub;
	NodeList SolSearch = Sol;

	for (auto& v : Sol) // graf induciran skupom vrhova Sol
	{
		WeightedGraph::adjIterator iterator(G, v);
		for (WeightedEdge* e = iterator.begin(); e != NULL; e = iterator.next())
		{
			if (Sol.count(e->w()) && Sol.count(e->v())) inducedSub.insert_edge(e);
		}
	}


	for (auto& v : Sol)
	{
		WeightedGraph impMST;
		set<WeightedEdge*> remEdges = inducedSub.remove_node(v);
		SolSearch.erase(v);

		if (!inducedSub.isConnected() || (  groups_of_node.count(v)>0   && !isCover(groups_of_node, SolSearch, k))) { SolSearch.insert(v); for (auto& e : remEdges) { inducedSub.insert_edge(e); }  continue; }
		improvedWeight = inducedSub.mst_prim(impMST);
		if (improvedWeight <= weight)
		{
			weight = improvedWeight;
			AntSol = impMST;
			
		}
		else
		{
			for (auto& e : remEdges) inducedSub.insert_edge(e);
			SolSearch.insert(v);
		}
		
	}
	return;
}

status_t GST_solver::solve_heuristic(string method)
{
	if (method == "ANTCOLONY")
	{
		std::cout << "\n------------------------------------------------------------------------" << std::endl;
		std::cout << "            ANT-COLONY SEARCH heuristic GST problem solver                    " << std::endl;
		std::cout << "------------------------------------------------------------------------" << std::endl;
	}

#ifdef USE_TIMER
	Timer heuristic_time;
	heuristic_time.start();
#endif


	// logging
# ifdef LOG
	std::ostringstream stringStream, logPath;
	stringStream << "heuristic_log_" << I.stp_file << "_ANT" << "_.log";
	std::cout << "Opening log file: " << stringStream.str() << std::endl;
	logPath << LOG_PATH << stringStream.str();

	heuristic_log.open(logPath.str().c_str());

	heuristic_log << "*** APPROXIMATION ***" << std::endl;
	heuristic_log << "** INSTANCE:" << std::endl;
	heuristic_log << "V " << I.G.V() << std::endl;
	heuristic_log << "E " << I.G.E() << std::endl;
	heuristic_log << "K " << I.number_of_groups() << std::endl;
	heuristic_log << std::endl;
	heuristic_log << "** LOG:" << std::endl;

#endif


	/* parameters setting */
	_index n = this->I.size();
	_index k = this->I.number_of_groups();
	dreal alpha = 0.5, beta = 10;
	DistanceMatrix dist(n, DistanceList(n, INFINITY));
	PredAccesMatrix prlist(n, PredAccesList(n, 0));
	this->I.G.all_pairs_shorthest_paths(dist, prlist);
	std::cout << "Computation of all shorthest paths finished." << endl;
	// odre�ivanje najmanje grupe groot
	_index groot = 0;
	size_t grootsz = I.nodes_of_group[groot].size();
	for (auto& g : this->I.nodes_of_group) if (g.second.size() < grootsz) { groot = g.first; grootsz = g.second.size(); }
	//cout << groot << " " << grootsz << endl;
	sreal estOpt = INFINITY;
	NodeList grootNodes = I.nodes_of_group[groot];

	for (auto& v : grootNodes)
	{
		GroupList nodes_of_group = I.nodes_of_group; // za svaku grupu sadr�i listu vrhova koji pripadaju toj grupi
		GroupList groups_of_node = I.groups_of_node; // za svaki vrh sadr�i listu grupa kojima taj vrh pripada
		//for (auto& x : groups_of_node) if (x.second.size() == 0) { groups_of_node.erase(x.first); } // obri�i Steinerove vrhove odmah na po�etku
		remove_vertex(nodes_of_group, groups_of_node, v); // izbaci pokrivene grupe i one vrhove koji vi�e nisu terminali


		map<Node, sreal> distToTree; // mapa udaljenosti terminalnih vrhova do trenutnog stabla
		map<Node, Node> brotherInTree;
		for (auto& x : groups_of_node) { distToTree[x.first] = dist[v][x.first]; brotherInTree[x.first] = v; }
		sreal vSol = cost_expanding(dist, prlist, nodes_of_group, groups_of_node, distToTree, brotherInTree, I.G, v);
		if (vSol < estOpt) estOpt = vSol;
	}

	std::cout << "Estimation of optimum: " << estOpt << endl<<endl;

	//if (I.groups_of_node.size() > 27) { cout << "STOP-6!" << I.groups_of_node.size() << endl; }


	sreal optSol = INFINITY;
	sreal avgSol = 0;
	_index nC=0, numOfColonies = 10;
	//parameter estimation
	
	
	sreal taumax = numOfColonies * (1/estOpt);
	sreal taumin = taumax / (2*n);


	vector<dreal> pheromones(n, 1);
	WeightedGraph theBestSolution;
	vector<sreal> colSolutions(numOfColonies, INFINITY);
	map<string, dreal> heusristicInfCont;
	_index countH = 0;
	_index antId = 0;

	while (nC < numOfColonies)
	{
		std::cout << "Ant colony "<<nC + 1<<" started." << endl;
		WeightedGraph colSolution;
		antId = 0;
		while (antId < 500)
		{


			//NodeList groups;
			//for (_index i = 0; i < this->I.number_of_groups(); i++) groups.insert(i);
			WeightedGraph AntSolMin;
			sreal AntSolMinWeight = INFINITY;

			for (NodeList::iterator v = grootNodes.begin(); v != grootNodes.end(); ++v)
			{

				WeightedGraph AntSol; // inicijalizacija praznog te�inskog grafa koji predstavlja rje�enje
				sreal AntSolWeight = 0;
				AntSol.insert_node(*v); // dodavanje korjenskog terminala u rje�enje O(log n)
				//if (I.groups_of_node.size() > 27) { cout << "STOP-5!" << I.groups_of_node.size() << endl; break; }
				NodeList Sol; //samo vrhovi da ne moram svaki put vra�ati listu vrhova grafa
				//if (I.groups_of_node.size() > 27) { cout << "STOP-4" << endl; break; }
				Sol.insert(*v); // O(log n)
				assert(!Sol.empty());
				//if (I.groups_of_node.size() > 27) { cout << "STOP-3!" << endl; break; }
				GroupList nodes_of_group = I.nodes_of_group; // za svaku grupu sadr�i listu vrhova koji pripadaju toj grupi
				GroupList groups_of_node = I.groups_of_node; // za svaki vrh sadr�i listu grupa kojima taj vrh pripada
				/*for (auto& x : groups_of_node)
				{
					if (x.second.size() == 0)
					{
						groups_of_node.erase(x.first);
					} // obri�i Steinerove vrhove odmah na po�etku
				}*/
				remove_vertex(nodes_of_group, groups_of_node, *v); // izbaci pokrivene grupe i one vrhove koji vi�e nisu terminali
				//if (I.groups_of_node.size() > 27) { cout << "STOP-2!" << endl; break; }

				map<Node, sreal> distToTree; // mapa udaljenosti terminalnih vrhova do trenutnog stabla
				map<Node, Node> brotherInTree;
				for (auto& x : groups_of_node) { distToTree[x.first] = dist[*v][x.first]; brotherInTree[x.first] = *v; }// na pocetku su to udaljenosti do vrha v jer je jedino taj vrh u stablu O(n)
				//if (I.groups_of_node.size() > 27) { cout << "STOP-1!" << endl; break; }

				while (nodes_of_group.size() > 0) //O(1)
				{

					NodeList nghSol = this->I.G.neighbourhood(Sol); // neefikasno je svaki put ra�unati susjedstvo O(n)
					// calculate key
					stringstream keyH;
					for (auto& x : Sol) keyH << x << "*";
					string key = keyH.str();
				

					map<Node, dreal> heuristicinf;
					for (auto& x : nghSol)
					{
						stringstream keyHx;
						keyHx << key;
						keyHx << x;
						string keyx = keyHx.str();
						map<string, dreal>::iterator calc = heusristicInfCont.find(keyx);
						if (calc != heusristicInfCont.end())
						{
							heuristicinf[x] = calc->second;
							countH++;
							//cout << "HELP" <<antId << endl;
						}
						else
						{
							sreal mind = INFINITY; // ra�unam udaljenost x do trenutnog rje�enja - ako je x terminal onda to ve� znam - iskoristiti kasnije
							for (auto& y : Sol) if (dist[x][y] < mind)
							{
								mind = dist[x][y];
							}

							sreal cost_ex = cost_expanding(dist, prlist, nodes_of_group, groups_of_node, distToTree, brotherInTree, I.G, x);
							/*if (cost_ex == 0)  heuristicinf[x] = 1e-5;
							else heuristicinf[x] = 1.0 / (mind * cost_ex); // pozivan Cost expanding algoritam O(n * (log n + k)));*/
							heuristicinf[x] = 1.0 / (mind + cost_ex);
							if(heusristicInfCont.size()<100000000000) heusristicInfCont[keyx] = heuristicinf[x];

						}
						
					}
					
					Node next = ant_random_sample(pheromones, heuristicinf, alpha, beta, nghSol);
					
					sreal mind = INFINITY;
					_index minInd = 0;
					WeightedEdge* addEdge = NULL;
					WeightedGraph::adjIterator it(I.G, next);
					for (WeightedEdge* e = it.begin(); e != NULL; e = it.next())
					{
						if (Sol.count(e->other(next)) && e->wt() < mind) {
							mind = e->wt(); minInd = e->other(next);
							addEdge = e;
						}
					}
					

					  // u narednom bloku dodajem sve vrhove na najkra�em putu od next do stabla, azuriram udaljenosti terminalnih vrhova do novog rjesenja i ostale strukture
					//Node pred = prlist[minInd][next];
					//if (I.groups_of_node.size() > 27) { cout << "STOP1!" << endl; break; }

					//WeightedEdge* addEdge = I.G.edge(pred, next);
					//if (addEdge == NULL) { cout<< "Outside " << minInd << " " << next<< " " << endl; }
					AntSol.insert_edge(addEdge);
					AntSolWeight += addEdge->wt();
					Sol.insert(next);
					//if (I.groups_of_node.size() > 27) { cout << "STOP2!" << endl; break; }
					NodeList SteinerP = remove_vertex(nodes_of_group, groups_of_node, next);
					distToTree.erase(next);
					brotherInTree.erase(next);
					for (auto& x : SteinerP) { distToTree.erase(x); brotherInTree.erase(x); }
					for (auto& x : groups_of_node) { if (distToTree[x.first] > dist[next][x.first]) { distToTree[x.first] = dist[next][x.first]; brotherInTree[x.first] = next; } }



					//for (auto& x : this->I.groups_of_node[next]) groupsUncovered.erase(x);
					//if (I.groups_of_node.size() > 27) { cout << "STOP3!" << endl; break; }
					/*while (pred != minInd)
					{
						WeightedEdge* addEdge = I.G.edge(prlist[minInd][pred], pred);
						AntSol.insert_edge(addEdge);
						if (addEdge == NULL) { cout<< "OnPath " << minInd << " " << pred << " " << endl; }
						AntSolWeight += addEdge->wt();
						Sol.insert(pred);
						NodeList SteinerP = remove_vertex(nodes_of_group, groups_of_node, pred);
						distToTree.erase(pred);
						brotherInTree.erase(pred);
						for (auto& x : SteinerP) { distToTree.erase(x); brotherInTree.erase(x); }
						for (auto& x : groups_of_node) { if (distToTree[x.first] > dist[pred][x.first]) { distToTree[x.first] = dist[pred][x.first]; brotherInTree[x.first] = pred; } }
						pred = prlist[minInd][pred];
					}*/
					
				}
				//if (I.groups_of_node.size() > 27) { cout << "STOP4!" << endl; break; }
				//GraphX::IO<WeightedGraph, WeightedEdge>::print(AntSol);
				
				//cout << AntSol.total_weight() << endl;
				//cout << "Ant " << antId << " finished search proccess: root = " << v << "; solution cost before lsa = " << AntSolWeight << endl;
				
				//std::cout << "Ant " << antId << " finished search proccess: root = " << v << "; solution cost after ls improvement = " << AntSolWeight << endl;
				//cout << "Ant search path on root " << v << " finished." << endl;
				//lsa(AntSol, Sol, I.G, AntSolWeight);
				//if (I.groups_of_node.size() > 27) { cout << "STOP5!" << endl; break; }
				lsr(AntSol, Sol, I.groups_of_node, I.G, AntSolWeight, k);
				//if (I.groups_of_node.size() > 27) { cout << "STOP6!" << endl; break; }
				if (AntSolWeight < AntSolMinWeight) { AntSolMin = AntSol; AntSolMinWeight = AntSolWeight; }
				//if (I.groups_of_node.size() > 27) { cout << "STOP7!" << endl; break; }
				//AntSol.~UndirectedGraph();
			}
			//std::cout << endl;
			//std::cout << "\t\tAnt " << antId << " finished search proccess with solution cost: " << AntSolMinWeight << endl << endl;
			if (AntSolMinWeight < colSolutions[nC]) { colSolutions[nC] = AntSolMinWeight; colSolution = AntSolMin; }
			antId++;
		}
		if (colSolutions[nC] < optSol) { optSol = colSolutions[nC]; theBestSolution = colSolution; }


		// pheromone update
		WeightedGraph::nodeIterator it(theBestSolution);
		for (Node v = it.begin(); v != NULL; v = it.next()) {
			if (pheromones[v] + (1 / optSol) <= taumax) { pheromones[v] += 1 / optSol; }
			else pheromones[v] = taumax;
		}


		//pheromone evaporation
		//dreal rho = ((dreal)rand()) / RAND_MAX;
		//rho = (9 * rho + 1) / 100;
		dreal rho = 0.01;
		for (int i = 0; i < n; i++)
		{
			if (pheromones[i] * (1 - rho) < taumin) { pheromones[i] = taumin; }
			else { pheromones[i] *= (1 - rho); }
		}
		avgSol += colSolutions[nC];
		
		std::cout << "Ant colony " << nC + 1 << " finished. Cost of the solution: " << colSolutions[nC] << endl << endl;
		// std::cout << "Container size: " << heusristicInfCont.size() << endl;
		nC++;
	}
	// stop  timing
#ifdef USE_TIMER
		heuristic_time.stop();
#endif

	std::cout << endl;
	std::cout << "Average  ACOGS: " << avgSol/numOfColonies <<endl;
	std::cout << "The best ACOGS: " << optSol <<" "<<endl;
	std::cout << "Number of accesss to heuristic container: " << countH << endl;
	//GraphX::IO<WeightedGraph, WeightedEdge>::print(theBestSolution);
	//WeightedGraph::nodeIterator it(theBestSolution);
	//for (Node v = it.begin(); v != NULL; v = it.next()) cout << v << endl;
	//cout << theBestSolution.isConnected() << endl;
	std::cout << endl;


	// logging
	heuristic_log << "APPROX:           " << optSol << std::endl;

	heuristic_log << "\n** TIME:" << std::endl;
	heuristic_log << "ANTS PER COLONY:           " << antId << std::endl;
	heuristic_log << "ANT COLONY RUNS:           " << nC << std::endl;
	heuristic_log << "AVG. TIME for ACOGS:       " << heuristic_time.secs()/nC << std::endl;
	heuristic_log << "TOTAL TIME:                " << heuristic_time.secs() << std::endl;
	heuristic_log << "** SOLUTION:" << " feasible " << std::endl;
	for (auto& e : x) if (e.second != 0) heuristic_log << e.first->v() << "--" << e.first->w() << ":" << e.first->wt() << std::endl;

	heuristic_log << flush;
	heuristic_log.close();

	return FEAS;
}
