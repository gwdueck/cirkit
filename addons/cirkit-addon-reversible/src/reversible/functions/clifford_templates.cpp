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

#include "clifford_templates.hpp"

#include <iostream>
#include <fstream>
#include <string>

namespace cirkit
{

bool is_CNOT_gate( const gate& g )
{
    return is_toffoli( g );
}

void Clifford_Template::read( std::ifstream &infile )
{
    Cliff_Gate cliffg;
    char gcode;
    int ngates_match, ngates_repl;
    infile >> num_qubits >> ngates_match >> ngates_repl;  
    for(int i = 0; i < ngates_match + ngates_repl; i++){
        infile >> gcode;
        while( gcode == ' ' ) infile >> gcode;
        cliffg.gtype = cliff_map.find(gcode)->second;
        infile >> cliffg.target;
        if( cliffg.gtype == CNOT )
        {
            cliffg.control = cliffg.target;
            infile >> cliffg.target;
        }
        if( i < ngates_match )
        {
            gates_matched.push_back( cliffg );
        }
        else
        {
            gates_replaced.push_back( cliffg );
        }
    }
}

void Clifford_Template::print( )
{
    std::cout << "num_qubits = " << num_qubits << " ";
    for (std::vector<Cliff_Gate>::iterator it = gates_matched.begin() ; it != gates_matched.end(); ++it)
    {
        std::cout << gate_name[it->gtype] << " " << it->target << " ";
    }
    std::cout << " ==> ";
    for (std::vector<Cliff_Gate>::iterator it = gates_replaced.begin() ; it != gates_replaced.end(); ++it)
    {
        std::cout << gate_name[it->gtype] << " " << it->target << " ";
    }
    std::cout << std::endl;
}

void Clifford_Template::clear( )
{
    gates_matched.clear();
    gates_replaced.clear();
}
}



// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
