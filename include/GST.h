#ifndef GST_INST_H
#define GST_INST_H

#include "DataTypes.h"
//string processing
#include <boost/algorithm/string.hpp>
// Gurobi optimizer
#include "gurobi_c++.h"
#include<unordered_set>


// typedefs 
typedef UndirectedGraph< WeightedEdge > WeightedGraph;
typedef std::set< WeightedEdge* > WeightedEdgeList;
typedef enum {type1, type2, type3} model_t;

// aux data structures

// logging & verbose
#define VERBOSE 
#define LOG
void mkDir(const char* path);

// path to instances
#define INSTANCE_PATH "./instances/"
#define LOG_PATH      "./logs/"

class GST_instance
{
	/*
	abstract interface to GST instance
	*/
	
public:
	GST_instance();
	GST_instance(filename stp_file);
	


	// IO
	status_t read_stp(filename stp_file);

	// interface to instance
	integer size()				  const; //returns the number of nodes
	integer number_of_groups()	  const; //returns the number of groups, or in the Steiner tree case, number of terminals
	integer number_of_edges()	  const; //returns the number of edges	
	integer bound_on_group_size() const; // return the upper bound on group size  
	

	// attributes
	bool	 isInstanceConnected();		// check if underlying graph is connected

	NodeList& get_groups(_index v);		// returns vector of groups in which  node v belongs
	NodeList& get_nodes(_index g);		// returns vector of nodes which belong to group g
	GroupList& get_all_groups();
	
	
	WeightedEdgeList	get_edges();
	NodeList			get_neighbours(_index v);
	bool				in_group(_index g, _index v);
	status_t			read_status;
	filename			stp_file;
	
	
	// data
protected:

	WeightedGraph G; 
	GroupList		  nodes_of_group;
	GroupList		  groups_of_node;
	//NodeList		  terminal_nodes; // nodes that belong to some group
	

	
	friend class GST_solver;
		
};


class EdgeWithConnComp : public WeightedEdge
{
protected:

	dreal dualConstrValue;
	dreal price;
	_index numOfDeltaSets;

public:
	ConnectedComponent<EdgeWithConnComp>* firstCut;
	ConnectedComponent<EdgeWithConnComp>* secondCut;// a set cuts that are defined by connected conponent
	integer indOfFirstCut;
	integer indOfSecondCut;

	EdgeWithConnComp() {}
	EdgeWithConnComp(int v, int w, double wt);
	EdgeWithConnComp(WeightedEdge* e);
	dreal calculatePrice();
	dreal get_price() const;
	void set_price(dreal value);
	void recalculateNumOfDSets();
	_index get_numOfDSets() const;
	void increaseDualConstrValue(dreal delta);
	void addConnectedComponent(ConnectedComponent<EdgeWithConnComp>* cc, _index ind);
	ConnectedComponent<EdgeWithConnComp>* getConnectedComponent(_index i) const;
	_index getIndOfConnComp(_index i) const;
	bool operator<(const EdgeWithConnComp& e) const;
	bool isOnCut(ConnectedComponent<EdgeWithConnComp>* cc);
	/*bool isInside(ConnectedComponent<EdgeWithConnComp>* cc);*/
	~EdgeWithConnComp();
};

typedef std::unordered_set<EdgeWithConnComp*> oe_t;

class SteinerTreeGraph : public UndirectedGraph<EdgeWithConnComp>
{
	typedef std::map<_index, ConnectedComponent<EdgeWithConnComp>*> cc_t;
	typedef PriorityQueue<EdgeWithConnComp>* edge_sampler_t;
	typedef std::map<ConnectedComponent<EdgeWithConnComp>*, std::set<EdgeWithConnComp*>> cc_edge_t;
	typedef std::map<_index, std::set<EdgeWithConnComp*>> n_edge_t;
protected:
	cc_t connComp;
	edge_sampler_t edgeSampler;
	cc_edge_t connCompToEdges;
	n_edge_t nodeToEdges;
	integer numberOfConnComp;
	oe_t orderOfEdges;
	_index countEdges;
public:
	SteinerTreeGraph() {}
	SteinerTreeGraph(GST_instance& SteinerTreeInstance);
	bool isConnected();
	EdgeWithConnComp* sampleEdge();
	dreal sampleAllTightEdges();
	void mergeComponents(EdgeWithConnComp* addedEdge);
	void blowUpEdges(dreal delta);
	void mergeWithSampledEdges();
	void refreshSampler();
	integer get_number_of_cc() const;
};



class GST_solver
{
	/*
		GST approximation algoritm solver.

	*/
	typedef std::map<WeightedEdge*, dreal> edge_var_t;
	typedef std::map< _index, map<_index, dreal > > flow_var_t;
	typedef UndirectedGraph<EdgeWithConnComp> gst_tree;
	typedef enum { ILP, LP, MILP } sol_t;


public:
	GST_solver();
	GST_solver(GST_instance& instance);

	// IO methods		
	status_t write_sol(filename sol_file);
	void     deleteSteinerLeaves();
	dreal    weightApprox();

	// presolve methods
	// check feasibility of data and etc.
	bool check_feasibility(gst_tree &T) const;
	bool check_terminals() const;


	// solvers
	status_t	solve_opt(model_t model); // Domagoj - Gurobi implementacija po ILP kompaktna mrežna formulacija, disertacija
	status_t	solve_approx(model_t model); // Slobodan - implementira algoritam [3]
	status_t	solve_heuristic(string method);
	status_t	solve_Steiner_approx();



protected:
	GST_instance& I; // abstract data structure, need to be implemented		
	// use Gurobi to solve LP formulation
	status_t GRB_solve_netflow_GR(_index root_group,sol_t solve_type);
	status_t GRB_solve_netflow_NR(_index r,sol_t solve_type);
	gst_tree* approximate_solution;

	// find roots
	_index find_root_group();
	// check solution integrality
	status_t check_integrality(edge_var_t *x, flow_var_t *f) const;

	// log files
	std::ofstream opt_log;
	std::ofstream approx_log;
	std::ofstream heuristic_log;


	
public:
	// netflow solution
	// fractional solution
	edge_var_t x;  
	flow_var_t f;
	dreal      frac_OPT;

	// optimal solution
	edge_var_t x_OPT;
	flow_var_t f_OPT;
	dreal      OPT;


	// tentative values
	dreal temp_obj_val;

	// Steiner Tree root
	_index	 root; 
	_index    root_group;
	
	// root processing
	NodeList root_candidates, non_root_groups;
	

	// terminal processing
	status_t define_terminals();
};



#endif