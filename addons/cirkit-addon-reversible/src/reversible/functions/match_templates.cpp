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

#include "match_templates.hpp"
#include <reversible/functions/clifford_templates.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/io/print_circuit.hpp>

namespace cirkit
{

bool gate_matches_template( const gate& g , const Cliff_Gate& g_temp, int qubits[] )
{
	bool match = is_Gate[ g_temp.gtype ]( g );
	if( match )
	{
		if( qubits[ g_temp.target ] == -1){
			qubits[ g_temp.target ] = g.targets().front();
		}
		else{
			match = qubits[ g_temp.target ] == g.targets().front();
		}
	}
	if( match && g_temp.gtype == CNOT ){
		if( qubits[ g_temp.control ] == -1){
			qubits[ g_temp.control ] = g.controls().front().line();
		}
		else{
			match = qubits[ g_temp.control ] == g.controls().front().line();
		}
	}
	return match;
}

/*
	Template has been matched at start in circ.
	Gates will be replace according to the template.
*/
void replace_matched_template( circuit& circ, Clifford_Template &ctempl, int qubit_map[], int start )
{
	for( int i = 0; i < ctempl.gates_matched.size(); i++ )
	{
		circ.remove_gate_at( start );
	}
	for( int i = 0; i < ctempl.gates_replaced.size(); i++ )
	{
		gate g;
		g.add_target( qubit_map[ ctempl.gates_replaced[i].target ] );
		switch ( ctempl.gates_replaced[i].gtype )
		{
			case H:
				g.set_type( hadamard_tag() );
				break;
			case T:
				g.set_type( pauli_tag( pauli_axis::Z, 4u, false ) );
				break;
			case Ts:
				g.set_type( pauli_tag( pauli_axis::Z, 4u, true ) );
				break;
			case S:
				g.set_type( pauli_tag( pauli_axis::Z, 2u, false ) );
				break;
			case Ss:
				g.set_type( pauli_tag( pauli_axis::Z, 2u, true ) );
				break;
			case Z:
				g.set_type( pauli_tag( pauli_axis::Z, 1u, false ) );
				break;
			case Y:
				g.set_type( pauli_tag( pauli_axis::Y, 1u, false ) );
				break;
			case RZ:
				std::cout << "ERROR in replace_matched_template(): RZ gate not implemented!\n";
				break;
//			case V: 	// not Clifford gates 
//			case Vs:	// not needed for now
			case X:
				g.set_type( toffoli_tag() );
				break;
			case CNOT: 
				g.set_type( toffoli_tag() );
				g.add_control( make_var( qubit_map[ ctempl.gates_replaced[i].control ], true ) );
				break;
			otherwise:
				std::cout << "ERROR in replace_matched_template(): gate not implemented!\n";
		}
		circ.insert_gate( start ) = g;
		start++;
	}
}

/*
	Check if there is a match of the template in the circuit.
	If there is, replace it.
*/
bool match_template( circuit& circ, Clifford_Template &ctempl )
{
	bool match; 
	int qubits[ ctempl.num_qubits ];

	int start = 0, len = ctempl.gates_matched.size();
	while( start + len  <= circ.num_gates() )
	{
		int i = 0;
		std::fill_n(qubits, ctempl.num_qubits, -1);
		match = gate_matches_template( circ[ start ], ctempl.gates_matched[i], qubits );
		i++;
		while( match && i < len )
		{
			int j = start + i ;
			bool next_match = false;
			while( j < circ.num_gates() && !next_match ){
				next_match = gate_matches_template( circ[ j ], ctempl.gates_matched[i], qubits );
				j++;
			}
			match = false;
			if( next_match )
			{
				j--;
				match = move_gate( circ,  start + i - 1, j);
			}
			i++;
		}
		if( match )
		{
			std::cout << "Matched! In match_template\n";
			ctempl.print();
			std::cout << "start = " << start << " ctempl.gates_matched.size() = " << ctempl.gates_matched.size() << "\n";
//			std::cout << circ ;
			replace_matched_template( circ, ctempl, qubits, start );
			return true;
		}
		start++;
	}
	return false;
}

/*
	Check if there is a match with any of the templates given.
	If there is one, then apply it and return true.
*/

bool match_any_template( circuit& circ, std::vector<Clifford_Template> &ctempls )
{
	int i = 0;
	for (std::vector<Clifford_Template>::iterator it = ctempls.begin() ; it != ctempls.end(); ++it)
	{
		if( match_template( circ, *it ) )
		{
			return true;
		}
	}
	return false;
}
}



// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
