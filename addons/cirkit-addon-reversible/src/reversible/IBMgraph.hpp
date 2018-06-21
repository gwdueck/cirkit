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
#include "functions/trans_path.hpp"

namespace cirkit
{

	bool **graph_adjacency = NULL;  // graph structure
	int **trans_cost = NULL;        // the cost of each cnot
	TransPath **trans_path = NULL;  // the transformation path or a given cnot

	int graph_size = 0;
	// std::vectosr<TransPath> path_list;

	bool read_graph( const std::string& filename );
	void print_graph( );
	void delete_graph( );
	void create_trans( bool verbose );

}
#endif /* IBMgraph_hpp */
