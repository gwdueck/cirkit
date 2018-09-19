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

/**
 * @file IBMgraph.hpp
 *
 * @brief Data structure to represent IBM quantum computer architecture
 *
 * @author Gerhard Dueck
 * @since  2.3
 */


#ifndef IBMgraph_hpp
#define IBMgraph_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <reversible/functions/trans_path.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/circuit.hpp>

namespace cirkit
{

	extern bool **graph_adjacency;  // graph structure
	extern int **trans_cost;        // the cost of each cnot
	extern TransPath **trans_path;  // the transformation path or a given cnot

	int graph_size = 0;
	
extern std::vector<TransPath> path_list; //-- it was giving memory error (now it is in cpp file)

	void allocate_data_stuctures();
	bool read_graph( const std::string& filename );
	bool write_to_file( const std::string& filename );
	bool read_from_file( const std::string& filename );
	void print_graph( );
	void print_matrix( );
	void delete_graph( );
	void create_trans( bool verbose );
    void the_mapping( circuit& circ_out, const circuit& circ_in );
    void expand_cnots( circuit& circ_out, const circuit& circ_in );

}
#endif /* IBMgraph_hpp */
