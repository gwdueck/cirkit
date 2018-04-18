/* CirKit: A circuit toolkit
 * Copyright (C) 2009-2015  University of Bremen
 * Copyright (C) 2015-2017  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "IBMgraph.hpp"
#include "functions/move_qubit.hpp"


namespace cirkit
{
    
    bool read_graph( const std::string& filename )
    {
        int v,w;
        std::fstream graphfile;
        graphfile.open ( filename, std::fstream::in );
        if ( !graphfile.is_open() )
        {
            return false;
        }
        graphfile >> graph_size;
        graph_adjacency = new bool*[graph_size];
        for(int i = 0; i < graph_size ; i++ )
        {
            graph_adjacency[i] = new bool[graph_size];
            for(int j = 0; j < graph_size; j++)
            {
                graph_adjacency[i][j] = false;
            }
        }
        while ( !graphfile.eof() )
        {
            graphfile >> v >> w;
            graph_adjacency[v][w] = true;
        }
        return true;
    }
    
    void print_graph( ){
        for( int i = 0; i < graph_size; i++ ){
            for( int j = 0; j < graph_size; j++ ){
                std::cout << (graph_adjacency[i][j] ? "X " : "- ");
            }
            std::cout << std::endl;
        }
    }
    
    void delete_graph( ){
        for( int i = 0; i < graph_size; i++ ){
            delete [] graph_adjacency[i];
            delete [] trans_path[i];
            delete [] trans_cost[i];
            
        }
        delete [] graph_adjacency;
        delete [] trans_path;
        delete [] trans_cost;
        graph_size = 0;
    }
    
    // When a path is found add it to global vector of paths.
    void find_all_paths( int v, int w, TransPath &tp, bool *visited )
    {
        bool done = false;
        for( int i = 0; i < graph_size; i++ )
        {
            if( graph_adjacency[v][w] )
            {
                tp.add( MoveQubit( nop, v, w ));
                path_list.push_back( tp );
                tp.remove_last();
                done = true;
            }
            if( graph_adjacency[w][v] )
            {
                tp.add( MoveQubit( flip, v, w ));
                path_list.push_back( tp );
                tp.remove_last();
                done = true;
            }
            if( done ) return; // no need to go further
            if( graph_adjacency[v][i] && !visited[i] )
            {
                tp.add( MoveQubit( cab, v, i ));
                visited[i] = true;
                find_all_paths( i,  w, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[i][v] && !visited[i] )
            {
                tp.add( MoveQubit( cba, v, i ));
                visited[i] = true;
                find_all_paths( i,  w, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[w][i] && !visited[i] )
            {
                tp.add( MoveQubit( tab, w, i ));
                visited[i] = true;
                find_all_paths( v,  i, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
            if( graph_adjacency[i][w] && !visited[i] )
            {
                tp.add( MoveQubit( tba, w, i ));
                visited[i] = true;
                find_all_paths( v,  i, tp, visited );
                tp.remove_last();
                visited[i] = false;
            }
        }
    }
    
    /* Precondition: the path_list contains all paths from v to w
     Postcondition: the best path and its cost will be stored in:
        - trans_cost[v][w]
        - trans_path[v][w]
     */
    
    void set_best_path(int v, int w )
    {
        TransPath tp, best_tp;
        int best_cost;
        best_tp = path_list[0];
        best_cost = best_tp.costPlus();
        for ( auto &p : path_list )
        {
            if(p.costPlus() < best_cost )
            {
                best_cost = p.costPlus();
                best_tp = p;
            }
        }
        trans_cost[v][w] = best_cost;
        trans_path[v][w] = best_tp;
    }
    
    void allocate_data_stuctures(){
        trans_cost = new int*[graph_size];
        trans_path = new TransPath*[graph_size];
        for(int i = 0; i < graph_size ; i++ )
        {
            trans_cost[i] = new int[graph_size];
            trans_path[i] = new TransPath[graph_size];
        }
    }
    
    void create_trans()
    {
        allocate_data_stuctures();
        path_list.clear();
        bool visited[graph_size];
        for( int i = 0; i < graph_size ; i++ )
        {
            visited[i] = false;
        }
        
        TransPath tp;
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
            {
                if( v == w )
                {
                    trans_cost[v][w] = 0;
                }
                else{
                    for( int i = 0; i < graph_size ; i++ ) visited[i] = false;
                    visited[v] = true;
                    visited[w] = true;
                    path_list.clear();
                    tp.clear();
                    find_all_paths( v,  w, tp, visited );
                    set_best_path( v,  w );
                }
            }
        }
        for( int v = 0; v < graph_size; v++)
        {
            for( int w = 0; w < graph_size; w++)
            {
                std::cout << trans_cost[v][w] << " ";
            }
            std::cout << std::endl;
        }
    }
}
