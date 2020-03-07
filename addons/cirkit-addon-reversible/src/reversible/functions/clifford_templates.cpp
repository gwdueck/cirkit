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
#include <reversible/functions/add_gates.hpp>

namespace cirkit
{

bool is_CNOT_gate( const gate& g )
{
    return is_toffoli( g ) && ( g.controls().size() == 1 );
}

void Cliff_Gate::convert_gate( const gate &g )
{
    control = -1;
    target = g.targets().front();
    if( is_hadamard( g ) )
    {
        gtype = H;
    }
    else if( is_T_gate( g ) )
    {
        gtype = T;
    } 
    else if( is_T_star_gate( g ) )
    {
        gtype = Ts;
    } 
    else if( is_S_gate( g ) )
    {
        gtype = S;
    } 
    else if( is_S_star_gate( g ) )
    {
        gtype = Ss;
    } 
    else if( is_Z_gate( g ) )
    {
        gtype = Z;
    } 
    else if( is_Y_gate( g ) )
    {
        gtype = Y;
    } 
    else if( is_X_gate( g ) )
    {
        gtype = X;
    } 
    else if( is_CNOT_gate ( g ) )
    {
        gtype = CNOT;
        control = g.controls().front().line();

    }
    else
    {
        std::cout << "ERROR gate not supported! in Cliff_Gate::convert_gate\n";
    }
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

void Clifford_Template::convert_circ( circuit &circ )
{
    num_qubits = circ.lines();
    Cliff_Gate cliffg;
    int n = circ.num_gates(), i = 0;
    for ( auto& gate : circ )
    {
        cliffg.convert_gate( gate );
        if( i <= n/2 )
        {
            gates_matched.push_back( cliffg );
        }
        else
        {
            gates_replaced.insert( gates_replaced.begin(), cliffg );
        }
        i++;
    }

}

void Clifford_Template::print( )
{
    std::cout  << num_qubits << " " << gates_matched.size() << " " << gates_replaced.size() << " ";
    std::string control_str;
    for (std::vector<Cliff_Gate>::iterator it = gates_matched.begin() ; it != gates_matched.end(); ++it)
    {
        control_str = ( it->gtype !=  CNOT )? "": std::to_string( it->control ) +  " ";
        std::cout << gate_name[it->gtype] << " " << control_str << it->target << " " ;
    }
    for (std::vector<Cliff_Gate>::iterator it = gates_replaced.begin() ; it != gates_replaced.end(); ++it)
    {
        control_str = (  it->gtype !=  CNOT )? "": std::to_string( it->control ) +  " ";
        std::cout << gate_name[it->gtype] << " " << control_str << it->target << " " ;
    }
    std::cout << std::endl;
}

void Clifford_Template::clear( )
{
    gates_matched.clear();
    gates_replaced.clear();
}

void append_cliff_gate( circuit &circ, Cliff_Gate gate ){
    std::vector<unsigned> control;
    switch( gate.gtype )
    {
        case H :
            append_hadamard( circ, gate.target );
            break;
        case T :
            append_pauli( circ, gate.target, pauli_axis::Z, 4u, true );
            break;
        case Ts :
            append_pauli( circ, gate.target, pauli_axis::Z, 4u, false );
            break;
        case S :
            append_pauli( circ, gate.target, pauli_axis::Z, 2u, true );
            break;
        case Ss :
            append_pauli( circ, gate.target, pauli_axis::Z, 2u, false );
            break;
        case Z :
            append_pauli( circ, gate.target, pauli_axis::Z );
            break;
        case Y :
            append_pauli( circ, gate.target, pauli_axis::Y );
            break;
        case X :
            append_not( circ, gate.target );
            break;
        case CNOT :
            control.push_back( gate.control );
            append_toffoli( circ, control, gate.target );
            break;
        default :
            std::cout << "error in append_cliff_gate " << gate.gtype << "\n";
    }
}
// convert the given template to an identity circuit
circuit Clifford_Template::convert_to_circ(  )
{
    circuit circ;
    circ.set_lines( num_qubits );
    std::vector<std::string> inputs;
    std::vector<constant> constants;
    for( int i = 0; i < num_qubits; i++ )
    {
        inputs.push_back( std::to_string(i) );
        constants.push_back( constant( false ));
    }
    circ.set_inputs( inputs );
    circ.set_constants( constants );
    for( auto & gate : gates_matched )
    {
        append_cliff_gate( circ, gate );
    }
    for (auto it = gates_replaced.rbegin(); it != gates_replaced.rend(); ++it)
    {
        append_cliff_gate( circ, *it);
    }
    return circ;
}

}



// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
