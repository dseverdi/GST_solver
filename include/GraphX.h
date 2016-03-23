#ifndef GRAPHX_H
#define GRAPHX_H


// libraries
#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <cmath>

// user-defined libraries
#include "PQ.h"

// namespaces
using namespace std;

// typedefs


#define EdgeList        std::set< Edge* >
#define Node	        _index 
#define NodeList        std::set<Node>
#define GroupList	    std::map<_index,NodeList>
#define DistanceList    std::vector<sreal>
#define PredAccesList   std::vector<Node>
#define DistanceMatrix  std::vector<DistanceList>
#define PredAccesMatrix std::vector<PredAccesList>


namespace GraphX
{
	// basic IO of a UndirectedGraph
	template <class Graph, class Edge>
	class IO
	{
	public:
		static void print(Graph &);
		static void read(Graph &, const char*);

		// slucajni generatori grafa
		static void randE(Graph &G, int E);
		static void randG(Graph &G, int E);
	};

	class WeightedEdge
	{
		/*
		class that implements weighted edge in an UndirectedGraph
		*/
	protected:
		_index v_id, w_id;
		sreal wght;
	public:
		// constructors
		WeightedEdge() {};
		WeightedEdge(_index v, _index w, sreal wt); 

		// access interface
		_index  v() const;
		_index  w() const;
		sreal  wt() const;
		bool   from(int v) const;
		_index  other(int w) const;

		// modifiers
		_index&  v();
		_index&  w();

		// specific operators
		bool operator<(const WeightedEdge& e) const;
		friend inline std::ostream& operator<<(std::ostream& buffer, const WeightedEdge& e)
		{
			buffer << e.v() << "-" << e.wt() << "-" << e.w() << " ";
			return buffer;
		}

	};

	template <class Edge>
	class UndirectedGraph
	{
		typedef std::map< _index, std::map<_index, Edge*> > AdjDict;
	
	protected:
		// Adjacency list
		AdjDict adj;

		// size of graph
		_index Vcnt, Ecnt;

	public:
		// iterators

		UndirectedGraph():Vcnt(0),Ecnt(0) {}
		UndirectedGraph(NodeList set_of_nodes);

		~UndirectedGraph();

		int V() const;
		int E() const;


		// modifiers
		  // add nodes
		void insert_node(Node v) { if (adj.count(v) == 0) { adj[v]; Vcnt++; } return; };
		EdgeList remove_node(Node v); // ToDo
		  // add edges
		void insert_edge(Edge* e);
		void remove_edge(Edge* e);

		integer degree_of_node(Node v);
		bool isLeaf(Node v);

		// access interface
		Edge* edge(_index v, _index w);
		Edge* edge(_index v, _index w) const;

		// containers
		EdgeList edges();
		NodeList nodes();
		NodeList neighbourhood(_index v);
		NodeList neighbourhood(NodeList& list);
		sreal total_weight() const;
		// graph properties ...
		bool isConnected();
		sreal shorthest_paths(Node s, Node t);
		void  single_source_shorthest_paths(DistanceList& dist, PredAccesList& predac, Node s);
		void  all_pairs_shorthest_paths(DistanceMatrix& dist, PredAccesMatrix& predac);
		sreal  mst_prim(UndirectedGraph<Edge>& mst);
				
		class adjIterator;
		class nodeIterator;
		friend class adjIterator;
		friend class nodeIterator;

	};


	template<class Edge>
	class ConnectedComponent
	{
	protected:
		NodeList nodes;
		//bool CandidateToBlowOn;
		//bool RootComponent;
	public:
		ConnectedComponent() {}
		ConnectedComponent(NodeList& nodes);
		_index size() const;
		void delete_node(_index node);
		void add_node(_index node);
		void add_component(ConnectedComponent<Edge>* cc);
		bool isOnCut(Edge* e)const; // returns true if edge e 
		bool isElement(Edge* e)const;
		bool isElement(Node v) const;
		//bool isCandidateToBlowOn() const;
		//bool isRootComponent() const;
	};

	


	
};


// template definition

// namespaces
using namespace GraphX;
using namespace std;


// WeightedEdge implementation
inline WeightedEdge::WeightedEdge(_index v, _index w, sreal wt) :v_id(v), w_id(w), wght(wt) {};
inline _index    WeightedEdge::v() const { return v_id; }
inline _index    WeightedEdge::w() const { return w_id; }
inline _index&   WeightedEdge::v() { return v_id; }
inline _index&   WeightedEdge::w() { return w_id; }
inline sreal    WeightedEdge::wt() const { return wght; }
inline bool     WeightedEdge::from(int v) const { return (v_id == v); }
inline _index    WeightedEdge::other(int v) const { if (from(v)) return w_id; else return v_id; }
inline bool     WeightedEdge::operator<(const WeightedEdge& e) const { return (this->wt() < e.wt()); }







template<class Edge> UndirectedGraph<Edge>::~UndirectedGraph()
{
	// remove dynamically stored edges
	
}



template<class Edge> inline int   UndirectedGraph<Edge>::V() const { return Vcnt; }
template<class Edge> inline int   UndirectedGraph<Edge>::E() const { return Ecnt; }



template<class Edge> inline void UndirectedGraph<Edge>::insert_edge(Edge *e)
{
	_index v = e->v(), w = e->w();
	 if (adj[v][w] == 0) Ecnt++;
	 adj[v][w] = adj[w][v] = e;
	
	 // update the size of V
	 Vcnt = (_index)adj.size();

	

}

template<class Edge> inline void UndirectedGraph<Edge>::remove_edge(Edge *e)
{
	_index v = e->v(), w = e->w();
	if (adj[v][w] != 0) Ecnt--;
	adj[v].erase(w);
	adj[w].erase(v);



}

template<class Edge> inline EdgeList UndirectedGraph<Edge>::remove_node(Node v)
{
	EdgeList remEdges;
	for (auto x : adj[v])
	{
		remEdges.insert(x.second);
		adj[x.first].erase(v);
		Ecnt--;
	}
	adj.erase(v);
	Vcnt--;
	return remEdges;
}

template<class Edge> integer UndirectedGraph<Edge>::degree_of_node(Node v)
{
	if (adj.count(v)) return adj[v].size();
	else
	{
		//cout << "Node with _index " << v << " does not exist. Return: -1";
		return -1;
	}
	
}

template<class Edge> bool UndirectedGraph<Edge>::isLeaf(Node v)
{
	return degree_of_node(v)==1;
}



template<class Edge> inline Edge* UndirectedGraph<Edge>::edge(_index v, _index w) 
{
	if(adj.count(v) && adj[v].count(w)) return adj[v][w];
	else {
		cout << "No such edge..." <<v<<" "<<w<< endl; return NULL;
	}
}

template<class Edge> inline Edge* UndirectedGraph<Edge>::edge(_index v, _index w) const
{
	if (adj.count(v) && adj[v].count(w)) return adj[v][w];
	else {
		cout << "No such edge..." << v << " " << w << endl; return NULL;
	}
}



template<class Edge> inline EdgeList UndirectedGraph<Edge>::edges() 
{
	EdgeList edge_list;
	
	for (auto& v: adj)
	{
		for (auto& e : v.second)
			if(e.second) edge_list.insert(e.second);
	}
	return edge_list;
}


template<class Edge> inline NodeList UndirectedGraph<Edge>::nodes()
{
	NodeList _nodes;
	for (auto& v : adj)  _nodes.insert(v.first);	

	return _nodes;
}


template<class Edge> inline NodeList UndirectedGraph<Edge>::neighbourhood(_index v) // ovo treba oiptimizirati takod da se susjedstvo vrha obilazi pomo�u iteratora...
{
	NodeList adj_vertices;
	
	for (auto e : adj[v])
	{
		if (e.second) adj_vertices.insert(e.second->other(v));
	}

	return adj_vertices;
}


template<class Edge> inline NodeList UndirectedGraph<Edge>::neighbourhood(NodeList& list)
{
	NodeList nbh;
	for (auto& v : list)
	{
		for (auto& e : adj[v])
		{
			if (!list.count(e.first)) nbh.insert(e.first);
		}
	}
	return nbh;

}
template<class Edge>
sreal UndirectedGraph<Edge>::total_weight() const
{
	sreal twt = 0;
	for (auto& v : adj)
	{
		for (auto& e : v.second)
		{
			twt += e.second->wt();
		}
	}
	return twt / 2;
}


template<class Edge> inline bool UndirectedGraph<Edge>::isConnected() 
{
	_index numOfVisited = 0;
	_index currentNode;
	// get list of nodes
	NodeList _nodes = nodes();
	std::queue<int> nodes_q;

	std::map<_index, color> colors;

	for (auto v : _nodes) colors[v] = WHITE;


	// first node
	_index _first_node = *(_nodes.begin());
	nodes_q.push(_first_node);

	// BFS test of connectivity
	while (!nodes_q.empty())
	{
		currentNode = nodes_q.front();
		nodes_q.pop();
		colors[currentNode] = BLACK;
		numOfVisited++;

		for (auto& adj_iter : adj[currentNode]) 
		{ 
			auto v = adj_iter.first;
			if (colors[v] == WHITE)
			{ 
				nodes_q.push(adj[currentNode][v]->other(currentNode));
				colors[v] = GRAY;
			}
		}
	}

	if (numOfVisited < Vcnt) return false;
	return true;
}



template <class Edge>
class UndirectedGraph<Edge>::adjIterator
{
	typedef std::map<_index, Edge* >						    adj_t;
	typedef typename adj_t::iterator 			     adj_iterator;

	
	UndirectedGraph<Edge> &G;
	int i, v;
	adj_iterator e_iter;

	
	
public:

	
	
	adjIterator(UndirectedGraph<Edge> &G, int v) : G(G), v(v), i(0) {  }

	Edge* begin()
	{
		e_iter = G.adj[v].begin();
		return e_iter->second;
		
	}
	Edge* next()
	{
		e_iter++;
		if (!end()) return e_iter->second;
		else return NULL;
	}
	bool end() 
	{

		if (e_iter!=G.adj[v].end())
		return false;
		else return true;
	}
};

template <class Edge>
class UndirectedGraph<Edge>::nodeIterator
{
	typedef std::map<Node, map<Node,Edge*> > adj_t;
	typedef typename adj_t::iterator node_iter;

	UndirectedGraph<Edge> &G;
	node_iter nit;



public:



	nodeIterator(UndirectedGraph<Edge> &G) : G(G) {  };

	Node begin()
	{
		nit = G.adj.begin();
		return nit->first;

	}

	Node next()
	{
		nit++;
		if (!end()) return nit->first;
		else return NULL;
	}
	bool end()
	{

		if (nit != G.adj.end())
			return false;
		else return true;
	}
};






template <class Graph, class Edge>
void GraphX::IO<Graph, Edge>::print(Graph &G)
{
	std::cout << "\n\r |V| = " << G.V() << ", |E| = " << G.E() << std::endl;
	typedef UndirectedGraph < WeightedEdge > WeightedGraph;
	WeightedGraph::nodeIterator itv(G);

	for (Node v = itv.begin(); v != NULL; v = itv.next())
	{
		std::cout.width(2); std::cout << v << ":";
		//if (G.E() == 0) return;
		typename Graph::adjIterator A(G, v);
		for (Edge* t = A.begin(); !A.end(); t = A.next())
		{
			std::cout.width(2); std::cout << std::setw(3) << t->other(v) << "|" << std::setw(3) << std::setprecision(3) << t->wt() << std::setw(3) << "     ";
		}
		std::cout << std::endl;
	}

};


// graph processing algorithms

template<class Edge> ConnectedComponent<Edge>::ConnectedComponent(NodeList& nodes)
{
	this->nodes = nodes;
	//if (nodes.count(root)) RootComponent = true;
	//else RootComponent = false;
	//CandidateToBlowOn = true;
}

template<class Edge> _index ConnectedComponent<Edge>::size() const
{
	return (_index)nodes.size();
}

template<class Edge> void ConnectedComponent<Edge>::add_node(_index node)
{
	nodes.insert(node);
}


template<class Edge> void ConnectedComponent<Edge>::delete_node(_index node)
{
	nodes.erase(node);
}

template<class Edge> bool ConnectedComponent<Edge>::isOnCut(Edge* e) const
{
	NodeList::iterator vis = nodes.find(e->v());
	NodeList::iterator wis = nodes.find(e->w());

	if (vis != nodes.end() && wis == nodes.end())
	{
		return true;
	}
	if (vis == nodes.end() && wis != nodes.end())
	{
		return true;
	}
	return false;
}


template<class Edge> bool ConnectedComponent<Edge>::isElement(Edge* e) const
{
	NodeList::iterator vis = nodes.find(e->v());
	NodeList::iterator wis = nodes.find(e->w());

	if (vis != nodes.end() && wis != nodes.end()) return true;
}

template<class Edge> bool ConnectedComponent<Edge>::isElement(Node v) const
{
	return (bool)nodes.count(v);
}

template<class Edge> void ConnectedComponent<Edge>::add_component(ConnectedComponent<Edge>* cc)
{
	nodes.insert(cc->nodes.begin(),cc->nodes.end());
}

/*template<class Edge> bool ConnectedComponent<Edge>::isRootComponent() const
{
	return RootComponent;
}

template<class Edge> bool ConnectedComponent<Edge>::isCandidateToBlowOn() const
{
	return CandidateToBlowOn;
}

template<class Edge> bool ConnectedComponent<Edge>::isCandidateToBlowOn()
{
	return CandidateToBlowOn;
}*/

#endif
template<class Object>
class Pair
{
public:
	Object v;
	sreal key;
	Pair(Object obj, sreal key) : v(obj), key(key) {};
	sreal get_price() const { return key; }
	void set_price(sreal k) { key = k; return; }
	template<class U> friend bool operator<=(Pair<U>& p1, Pair<U>& p2) { return p1.key <= p2.key; }
	template<class U> friend bool operator<(Pair<U>& p1, Pair<U>& p2) { return p1.key < p2.key; }
	~Pair() {};
};


/*template<class Object, class Priority>
class ObjectWithPriority : class Object
{
public:
	Object v;
	sreal key;
	Pair(Object obj, sreal key) : v(obj), key(key) {};
	sreal get_price() const { return key; }
	void set_price(sreal k) { key = k; return; }
	template<class U> friend bool operator<=(Pair<U>& p1, Pair<U>& p2) { return p1.key <= p2.key; }
	template<class U> friend bool operator<(Pair<U>& p1, Pair<U>& p2) { return p1.key < p2.key; }
	~Pair() {};
};*/



template<class Edge>
void UndirectedGraph<Edge>::single_source_shorthest_paths(DistanceList& dist, PredAccesList& plist, Node s)
{
	size_t n = dist.size();
	if (n != this->Vcnt) return;
	PriorityQueue<Pair<Node>> Q(n);
	vector<Pair<Node>*> pairs;
	for (Node i = 0; i < n; i++) { pairs.push_back(new Pair<Node>(i, dist[i])); Q.push(pairs[i]); }
	while (!Q.isEmpty())
	{
		Pair<Node>* u = Q.top();
		Q.pop();
		//onQ[u->v] = false; 
		for (auto it = adj[u->v].begin(); it != adj[u->v].end(); ++it)
		{
			if (dist[it->first] > it->second->wt() + dist[u->v])
			{
				dist[it->first] = it->second->wt() + dist[u->v];
				plist[it->first] = u->v;
				Q.changeKey(pairs[it->first], dist[it->first]);
			}
		}
	}
	
	return;

}

template<class Edge>
sreal UndirectedGraph<Edge>::shorthest_paths(Node s, Node t)
{
	DistanceList dist(this->Vcnt, INFINITY);
	PredAccesList prlist(this->Vcnt, 0);
	dist[s] = 0;
	this->single_source_shorthest_paths(dist, prlist, s);
	return dist[t];
		}



template<class Edge>
void UndirectedGraph<Edge>::all_pairs_shorthest_paths(DistanceMatrix& dist, PredAccesMatrix& plist)
{
	size_t n = dist.size();
	if (n != this->Vcnt) return; 
	for (int i = 0; i < n; i++)
	{
		dist[i][i] = 0;
		single_source_shorthest_paths(dist[i], plist[i], i);
	}
	
	return;
}

template<class Edge>
sreal UndirectedGraph<Edge>::mst_prim(UndirectedGraph<Edge>& mst)
{

	size_t n = this->Vcnt;
	sreal mst_weight = 0;
	PriorityQueue<Pair<Node>> Q(n);
	map<Node, Node> predac;
	map<Node, bool> included;
	for (auto& x : adj) included[x.first] = false;
	_index root = (adj.begin())->first;
	included[root] = true;
	map<Node, Pair<Node>*> pairs;
	for (auto& x : adj) { pairs[x.first] = new Pair<Node>(x.first, INFINITY); Q.push(pairs[x.first]); }
	Q.deleteElement(pairs[root]);
	_index u = root;
	predac[u] = u;


	while (!Q.isEmpty())
	{
		for (
				typename map<_index,Edge*>::iterator it = adj[u].begin();
				it != adj[u].end();
				++it
				)
		{
			sreal wt = it->second->wt();
			if (!included[it->first] && wt < pairs[it->first]->get_price())
			{
				Q.changeKey(pairs[it->first], wt);
				predac[it->first] = u;
			}
		}

		Pair<Node>* wnode = Q.top();
		u = wnode->v;
		included[u] = true;
		Q.pop();
		mst_weight += (adj[u][predac[u]])->wt();
		mst.insert_edge(adj[u][predac[u]]);
	}

	return mst_weight;
}
