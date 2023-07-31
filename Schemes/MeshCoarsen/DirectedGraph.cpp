#include "DirectedGraph.hpp"

// ----------------------------------------------------------------------------
//                            DirectedEdge
// ----------------------------------------------------------------------------

DirectedEdge::DirectedEdge() : a(0), b(0) {}
DirectedEdge::~DirectedEdge() {}
DirectedEdge::DirectedEdge(std::size_t _a, std::size_t _b) : a(_a), b(_b) {}

DirectedEdge DirectedEdge::anti_edge()
{
    return DirectedEdge(b, a);
}

bool operator==(const DirectedEdge &lhs, const DirectedEdge &rhs)
    {
        return ((lhs.a == rhs.a) && (lhs.b == rhs.b));
    }

// ----------------------------------------------------------------------------
//                            DirectedCell
// ----------------------------------------------------------------------------

DirectedCell::DirectedCell() {}
DirectedCell::DirectedCell(std::size_t cell_num) {part.push_back(cell_num);}
DirectedCell::~DirectedCell() {}
DirectedCell::DirectedCell(std::vector<DirectedEdge> &_T, std::size_t cell_num) : T(_T) 
{
	part.push_back(cell_num);
}

void DirectedCell::add_edge(DirectedEdge &e)
{
    T.push_back(e);
}

void DirectedCell::remove_edge(std::size_t pos)
{
    T.erase(std::begin(T) + pos);
}

void DirectedCell::append_cell(DirectedCell &K)
{
    T.insert(std::end(T), std::begin(K.T), std::end(K.T));
    part.insert(std::end(part), std::begin(K.part), std::end(K.part));
}

void DirectedCell::remove_duplicate_edges()
{    
    std::size_t iT = 0;
    while (iT < T.size() - 1)
    {
        for (std::size_t jT = iT + 1; jT < T.size(); ++jT)
        {
            if (T[iT] == T[jT].anti_edge())
            {
                remove_edge(jT);
                remove_edge(iT);
                iT--;
                break;
            }
        }
        iT++;
    }
}

bool DirectedCell::is_ordered()
{
    for (std::size_t it = 0; it < T.size() - 1; ++it)
    {
        if (T[it].b != T[it + 1].a)
        {
            return false;
        }
    }
    if (T[T.size() - 1].b != T[0].a)
    {
        return false;
    }
    return true;
}

std::vector<DirectedEdge> DirectedCell::order_edges()
{
	std::vector<DirectedEdge> edges_to_remove;    
    for (std::size_t iT = 0; iT < T.size() - 1; iT++)
    {
    	if (T[iT].b != T[iT + 1].a)
        {
		    for (std::size_t jT = iT + 2; jT < T.size(); jT++)
		    {
		        if (T[iT].b == T[jT].a)
		        {
		        	std::swap(T[iT+1], T[jT]);
		            break;
		        }
		    }
        }
    }

    if (!is_ordered())
    {
    	/* If edges are not ordered - then cell is not simply connected.
    	 * Cell will consist of disconnected loops.
    	 * Find all loops, set the cell to the longest loop and add
    	 * other loops to be removed. This will only fail if the outer 
    	 * most loop (the new cell) is not the longest loop. This is
    	 * very unlikely. In this case the whole graph will not be 
    	 * simply connected, and a loop region will be missing.
    	 */
    	std::cout << "Located non-simply connected cell\n";
    	std::vector<std::vector<DirectedEdge>> sub_vectors;
    	std::size_t pos = 0;    	
    	for(std::size_t i = 0; i < T.size() - 1; i++)
    	{
    		if(T[i].b != T[i+1].a)
    		{
    			std::vector<DirectedEdge> sub_vec(std::begin(T) + pos, std::begin(T) + i + 1);
    			pos = i+1;
    			sub_vectors.push_back(sub_vec);
    		}
    	}
    	std::vector<DirectedEdge> sub_vec(std::begin(T) + pos, std::end(T));
    	sub_vectors.push_back(sub_vec);
    	
    	std::vector<DirectedEdge> max_vec = sub_vectors[0];
    	for(std::size_t i = 1; i < sub_vectors.size(); i++)
    	{
    		if(sub_vectors[i].size() > max_vec.size())
    		{
    			max_vec = sub_vectors[i];
    		}
    	}
    	T = max_vec;
    	for(std::size_t i = 0; i < sub_vectors.size(); i++)
    	{
    		if(sub_vectors[i] != max_vec)
    		{
    			edges_to_remove.insert(std::end(edges_to_remove), std::begin(sub_vectors[i]), std::end(sub_vectors[i]));
    		}
    	}
    }
    
    return edges_to_remove;
}

bool DirectedCell::check_nodes()
{
	std::vector<size_t> as;
    std::vector<size_t> bs;
    for(size_t i = 0; i < T.size(); i++)
    {
    	if((find(begin(as), end(as), T[i].a) != end(as)) || (find(begin(bs), end(bs), T[i].b) != end(bs)))
    	{
    		return false;;
    	}
    	as.push_back(T[i].a);
    	bs.push_back(T[i].b);
    }
    for(std::size_t i = 0; i < T.size(); i++)
    {
        for (std::size_t j = i + 1; j < T.size(); ++j)
        {
            if (as[i] > as[j])
            {
                std::swap(as[i] , as[j]);
            }
            if (bs[i] > bs[j])
            {
                std::swap(bs[i] , bs[j]);
            }
        }
    }
    if(as != bs)
    {
    	return false;
    }
    return true;
}

bool DirectedCell::has_edge(DirectedEdge e)
{
    for (std::size_t it = 0; it < T.size(); it++)
    {
        if(T[it] == e)
        {
        	return true;
        }
    }
    return false;
}

void DirectedCell::print()
{
    for (std::size_t it = 0; it < T.size(); it++)
    {
        std::cout << T[it].a << "," << T[it].b << " ";
    }
    std::cout << "\n";
}

// ----------------------------------------------------------------------------
//                            NodeArray
// ----------------------------------------------------------------------------

NodeArray::NodeArray() {}
NodeArray::~NodeArray() {}
NodeArray::NodeArray(Array &_A) : A(_A) {}
bool NodeArray::node_exists(std::size_t node)
{
    for (std::size_t i = 0; i < A.size(); i++)
    {
        if (find(begin(A[i]), end(A[i]), node) != end(A[i]))
        {
            return true;
        }
    }
    return false;
}

void NodeArray::renum_nodes(std::size_t node)
{
    for (std::size_t i = 0; i < A.size(); i++)
    {
        for (std::size_t j = 0; j < A[i].size(); j++)
        {
            if (A[i][j] > node)
            {
                A[i][j]--;
            }
        }
    }
}

void NodeArray::print(std::ofstream *out, std::size_t width)
{
    *out << "cells"
         << "\n";
    *out << A.size() << "\n";;
    for (std::size_t cell_i = 0; cell_i < A.size(); cell_i++)
    {
        *out << std::setw(width) << A[cell_i].size();
        for (std::size_t node_j = 0; node_j < A[cell_i].size(); node_j++)
        {
            *out << std::setw(width) << A[cell_i][node_j];
        }        
        *out << "\n";
    }
}

// ----------------------------------------------------------------------------
//                            DirectedGraph
// ----------------------------------------------------------------------------

DirectedGraph::DirectedGraph() {}
DirectedGraph::~DirectedGraph() {}
DirectedGraph::DirectedGraph(std::vector<DirectedCell> &_G) : G(_G) {}

void DirectedGraph::add_cell(DirectedCell &T)
{
    G.push_back(T);
}

void DirectedGraph::remove_cell(std::size_t pos)
{
    G.erase(std::begin(G) + pos);
}

bool DirectedGraph::test_graph()
{
    for (std::size_t i = 0; i < G.size(); ++i)
    {
        if (!(G[i].check_nodes()))
        {
            std::cerr << "Cell " + std::to_string(i) + " has repeated or lone nodes\n";
            G[i].print();
            return false;
        }
        
        if (!(G[i].is_ordered()))
        {       	
            std::cerr << "Cell " + std::to_string(i) + " is not ordered\n";
            G[i].print();
            return false;
        }
    }
    return true;
}

void DirectedGraph::randomise()
{
    // std::srand(std::time(nullptr));
    // std::random_shuffle(begin(G), end(G));
    std::random_device rd;
    std::mt19937 g(rd());
 
    std::shuffle(G.begin(), G.end(), g);
}

void DirectedGraph::order()
{
    for(std::size_t i = 0; i < G.size(); i++)
    {
        for (std::size_t j = i + 1; j < G.size(); ++j)
        {
            if (G[i].T.size() > G[j].T.size())
            {
                std::swap(G[i] , G[j]);
            }
        }
    }
}

void DirectedGraph::coarsen()
{
    for (std::size_t it = 0; it < G.size() - 1; it++)
    {
    	// Merge for any number of edges shared
//        for (std::size_t jt = it + 1; jt < G.size(); jt++)
//        {
//        	
//        	std::vector<size_t> edges_shared;
//            if (can_merge(it, jt, edges_shared))
//            {
//                merge_cells(it, jt, edges_shared);
//                break;
//            }
//        }

		// Routine to merge largest number of shared edges
		bool flag = true;
    	std::vector<std::vector<size_t>> shared_edges;
        for (std::size_t jt = it + 1; jt < G.size(); jt++)
        {
        	std::vector<size_t> edges_shared;
            if (can_merge(it, jt, edges_shared))
            {
                if(edges_shared.size() >= G[it].T.size() - 2)
                {
                	merge_cells(it, jt, edges_shared);
                	flag = false;
                	break;
                } 
                else
                {
                	shared_edges.push_back(edges_shared);
                }
            }
            else
            {
            	std::vector<size_t> null_vec;
            	shared_edges.push_back(null_vec);
            }
        }
        if(flag)
        {
        	std::size_t pos = 0;
        	for(std::size_t jt = 1; jt < shared_edges.size(); ++jt)
        	{
        		if(shared_edges[jt].size() > shared_edges[pos].size())
        		{
        			pos = jt;
        		}
        	}
        	if(shared_edges[pos].size() > 0)
        	{
        		merge_cells(it, pos + it + 1, shared_edges[pos]);        	
        	}
        }
    }
}

NodeArray DirectedGraph::graph_to_array()
{
    Array array(G.size());
    for (size_t it = 0; it < G.size(); ++it)
    {
        for (size_t jt = 0; jt < G[it].T.size(); ++jt)
        {
            array[it].push_back(G[it].T[jt].a);
        }
    }

    return NodeArray(array);
}

std::string DirectedGraph::get_partition()
{
	std::string partition = "partition\n";
	for(size_t i = 0; i < G.size(); i++)
	{
		for(size_t j = 0; j < G[i].part.size(); j++)
		{
			partition += std::to_string(G[i].part[j]);
			partition += " ";
		}
		partition += "\n";
	}
	return partition;
}

void DirectedGraph::plotfile(std::ofstream *out, std::vector<Eigen::VectorXd> & vertices)
{
	std::vector<DirectedEdge> edges;
	for(auto & cell : G)
	{
		for(auto & edge : cell.T)
		{
			if( find(std::begin(edges), std::end(edges), edge.anti_edge()) == std::end(edges) )
			{
				edges.push_back(edge);
			}		
		}
	}
	for(auto & edge : edges)
	{
		*out << vertices[edge.a - 1](0) << " " << vertices[edge.a - 1](1) << std::endl;
		*out << vertices[edge.b - 1](0) << " " << vertices[edge.b - 1](1) << std::endl;
		*out << std::endl;
	}
}

bool DirectedGraph::can_merge(std::size_t pos1, std::size_t pos2, std::vector<size_t> &edges_shared)
{
    std::vector<DirectedEdge> cell1 = G[pos1].T;
    std::vector<DirectedEdge> cell2 = G[pos2].T;
    bool vert_flag = true;
    for (std::size_t i = 0; (i < cell1.size()) && vert_flag; i++)
//    for (std::size_t i = 0; i < cell1.size(); i++)
    {
        for (std::size_t j = 0; j < cell2.size(); j++)
        {
        	if(cell1[i].a == cell2[j].b) //tail of edge1 is head of edge2
        	{
        		if(cell1[i].b == cell2[j].a) //are they same edge?
        		{
        			edges_shared.push_back(i);
        			edges_shared.push_back(j + cell1.size());
        		} 
        		else 
        		{
        			std::size_t ipos = (i != 0 ? i - 1 :  cell1.size() - 1);
			        std::size_t jpos = (j != cell2.size() - 1 ? j + 1 : 0);
			        if( !( cell1[ipos].a == cell2[jpos].b ) ) //is the edge prior to edge1 the edge past edge2?
			        {
						vert_flag = false;
						break;
					}
        		}
        	}
        }
    }
    bool edge_flag = (edges_shared.size() > 0);
//    if(edge_flag && (!vert_flag))
//    {
//    	std::cout << "Connection at non-edge avoided" << "\n";
//    }
    return vert_flag && edge_flag;
}

void DirectedGraph::merge_cells(std::size_t pos1, std::size_t pos2, std::vector<size_t> &edges_shared)
{
    G[pos1].append_cell(G[pos2]);
    remove_cell(pos2);
    std::sort(std::begin(edges_shared), std::end(edges_shared), std::greater<size_t>());
    for(std::size_t ie = 0; ie < edges_shared.size(); ie++)
    {
    	G[pos1].remove_edge(edges_shared[ie]);
    }
    
    std::vector<DirectedEdge> to_remove = G[pos1].order_edges();    
    if(to_remove.size() != 0)
    {
    	for(size_t i = 0; i < to_remove.size(); i++)
    	{
    		size_t j = 0;
    		while(j < G.size())
    		{
    			if(j != pos1 && G[j].has_edge(to_remove[i].anti_edge()))
    			{
    				G[pos1].part.insert(std::end(G[pos1].part), std::begin(G[j].part), std::end(G[j].part));
    				remove_cell(j);
    				if(j < pos1)
    				{
    					pos1--;
    				}
    				j--;
    			}
    			j++;
    		}    		
    	}
    	test_graph();
    	std::cout << "Successfully connected non-simply connected cell\n";
    }
}

bool DirectedGraph::random_bool(double prob)
{
    std::srand(std::time(nullptr));
    std::bernoulli_distribution d(prob);
    return d(rand_engine);
}
