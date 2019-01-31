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

#include "alex.hpp"

#include <cmath>
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <iostream>
#include <core/utils/program_options.hpp>
#include <boost/program_options.hpp>
#include <reversible/functions/ibm_helper.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <reversible/functions/add_circuit.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/utils/matrix_utils.hpp>
#include <alice/rules.hpp>
#include <reversible/functions/clear_circuit.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/target_tags.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/find_lines.hpp>
#include <reversible/mapping/nct_mapping.hpp>
#include <reversible/functions/remove_dup_gates.hpp>

namespace cirkit
{

using boost::program_options::value;
using matrix = std::vector<std::vector<unsigned>>;

alex_command::alex_command( const environment::ptr& env )
	: cirkit_command( env, "Random projects" )
{
	opts.add_options()
    ( "interchange,i",  "interchange controls" )
    ( "clifford,c",  "transform to Clifford+T" )
    ( "direct,d",  "transform using Toffoli with 17 Clifford+T gates" )
    ( "permutations,p",  "try_all_permutations permutations (only 5 qubits)" )
    ( "single_Toffoli,t",  "try different V gates for a single Toffoli" )
    ( "remove,r",  "remove_dup_gates" )
	;

}

bool interchange = false;
bool remove = false;
bool direct = false;

command::rules_t alex_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

void toffoli_to_v_lower_cost(circuit& circ_out, gate g, matrix& cnot_costs)
{
    std::vector<unsigned int> controls_a, controls_b;
    controls_a.clear(); controls_b.clear();
	
    unsigned a, b, c, cost_a = 0, cost_b = 0;
    unsigned polarity_a, polarity_b;
    polarity_a = g.controls().front().polarity();
    polarity_b = g.controls().back().polarity();

    a = g.controls().front().line();
    b = g.controls().back().line();
    c = g.targets().front();

	cost_a += cnot_costs[a][b] * 2;
   	if( cnot_costs[a][c] < cnot_costs[c][a] )
		cost_a += cnot_costs[a][c];
	else
		cost_a += cnot_costs[c][a];
	if( cnot_costs[b][c] < cnot_costs[c][b] )
		cost_a += cnot_costs[b][c] * 4;
	else
		cost_a += cnot_costs[c][b] * 4;

	cost_b += cnot_costs[b][a] * 2;
   	if( cnot_costs[b][c] < cnot_costs[c][b] )
		cost_b += cnot_costs[b][c];
	else
		cost_b += cnot_costs[c][b];
	if( cnot_costs[a][c] < cnot_costs[c][a] )
		cost_b += cnot_costs[a][c] * 4;
	else
		cost_b += cnot_costs[c][a] * 4;

	// std::cout << "cost_a: " << cost_a << " cost_b: " << cost_b << std::endl;
  	if ( cost_b < cost_a )
	{
    	controls_a.push_back(g.controls().back().line());
    	controls_b.push_back(g.controls().front().line());
	}
	else
	{
		controls_a.push_back(g.controls().front().line());
    	controls_b.push_back(g.controls().back().line());	
	}
	
	if( polarity_a == 1 && polarity_b == 1 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), false );
	}
	else if( polarity_a == 0 && polarity_b == 0 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), false );
		append_not( circ_out, g.targets().front() );
	}
	else if( polarity_a == 1 && polarity_b == 0 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), cost_b < cost_a ? false : true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), cost_b < cost_a ? true : false );
	}
	else if( polarity_a == 0 && polarity_b == 1 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), cost_b < cost_a ? true : false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), cost_b < cost_a ? false : true );
	}
}

void toffoli_to_v(circuit& circ_out, gate g)
{
    std::vector<unsigned int> controls_a, controls_b;
    std::vector<unsigned int> polarities;
    polarities.clear(); controls_a.clear(); controls_b.clear();
	polarities.push_back(g.controls().front().polarity());
    polarities.push_back(g.controls().back().polarity());

  	if ( interchange && polarities.front() == polarities.back() )
	{
    	controls_a.push_back(g.controls().back().line());
    	controls_b.push_back(g.controls().front().line());
	}
	else
	{
		controls_a.push_back(g.controls().front().line());
    	controls_b.push_back(g.controls().back().line());	
	}
	
	if( polarities.front() == 1 && polarities.back() == 1 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), false );
	}
	else if( polarities.front() == 0 && polarities.back() == 0 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), false );
		append_not( circ_out, g.targets().front() );
	}
	else if( polarities.front() == 1 && polarities.back() == 0 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), false );
	}
	else if( polarities.front() == 0 && polarities.back() == 1 )
	{
		append_v( circ_out,  controls_a, g.targets().front(), false );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), true );
		append_toffoli( circ_out, controls_a, controls_b.front() );
		append_v( circ_out,  controls_b, g.targets().front(), true );
	}
}

void v_to_clifford(circuit& circ_in, circuit& circ_out, gate g)
{
    const auto& tag = boost::any_cast<v_tag>( g.type() );
	if ( tag.adjoint )
	{
		append_hadamard( circ_out, g.targets().front() );
		append_toffoli( circ_out, g.controls(), g.targets().front() );
		append_pauli( circ_out,  g.targets().front(), pauli_axis::Z, 4u, false );
		append_toffoli( circ_out, g.controls(), g.targets().front() );
		append_pauli( circ_out,  g.controls()[0].line(), pauli_axis::Z, 4u, true );
		append_pauli( circ_out,  g.targets().front(), pauli_axis::Z, 4u, true );
		append_hadamard( circ_out, g.targets().front() );
	}
	else
	{
		append_hadamard( circ_out, g.targets().front() );
		append_toffoli( circ_out, g.controls(), g.targets().front() );
		append_pauli( circ_out,  g.targets().front(), pauli_axis::Z, 4u, true );
		append_toffoli( circ_out, g.controls(), g.targets().front() );
		append_pauli( circ_out,  g.controls()[0].line(), pauli_axis::Z, 4u, false );
		append_pauli( circ_out,  g.targets().front(), pauli_axis::Z, 4u, false );
		append_hadamard( circ_out, g.targets().front() );
	}
}

void toffoli_to_clifford(circuit& circ, gate g)
{
	std::vector<unsigned> ga, gb;
	if( g.controls().front().line() < g.controls().back().line() )
	{
		ga.push_back(g.controls().front().line());
		gb.push_back(g.controls().back().line());
	}
	else
	{
		ga.push_back(g.controls().back().line());
		gb.push_back(g.controls().front().line());
	}

	append_hadamard( circ, g.targets().front() );
	append_toffoli( circ, gb, g.targets().front() );
	append_pauli( circ,  g.targets().front(), pauli_axis::Z, 4u, true );
	append_toffoli( circ, gb, g.targets().front() );
	append_pauli( circ,  gb.front(), pauli_axis::Z, 4u, false );
	append_toffoli( circ, ga, gb.front() );
	append_toffoli( circ, gb, g.targets().front() );
	append_pauli( circ,  g.targets().front(), pauli_axis::Z, 4u, false );
	append_toffoli( circ, gb, g.targets().front() );
	append_pauli( circ,  gb.front(), pauli_axis::Z, 4u, true );
	append_toffoli( circ, ga, gb.front() );
	append_toffoli( circ, ga, g.targets().front() );
	append_pauli( circ,  g.targets().front(), pauli_axis::Z, 4u, true );
	append_toffoli( circ, ga, g.targets().front() );
	append_pauli( circ,  ga.front(), pauli_axis::Z, 4u, false );
	append_pauli( circ,  g.targets().front(), pauli_axis::Z, 4u, false );
	append_hadamard( circ, g.targets().front() );
}

void transform_to_v(circuit& circ_in, circuit& circ_out, matrix& cnot_costs)
{
	copy_metadata(circ_in, circ_out);
	for ( const auto& gate : circ_in )
    {
		if ( is_toffoli( gate ) )
        {
            if(direct)
            	toffoli_to_clifford(circ_out, gate);
            else if( gate.controls().size() == 2 )
          	 	toffoli_to_v_lower_cost(circ_out, gate, cnot_costs);
          	 	// toffoli_to_v(circ_out, gate, cnot_costs);
            else
            	circ_out.append_gate() = gate;
		}
		else
			circ_out.append_gate() = gate;
	}
}

circuit Transform_to_v(circuit& circ_in, matrix& cnot_costs)
{
	circuit circ_out;
	copy_metadata(circ_in, circ_out);
	for ( const auto& gate : circ_in )
    {
		if ( is_toffoli( gate ) )
        {
            if( gate.controls().size() == 2 )
          	 	toffoli_to_v_lower_cost(circ_out, gate, cnot_costs);
            else
            	circ_out.append_gate() = gate;
		}
		else
			circ_out.append_gate() = gate;
	}
	return circ_out;
}

void transform_to_clifford(circuit& circ_in, circuit& circ_out)
{
	copy_metadata(circ_in, circ_out);
	for ( const auto& gate : circ_in )
    {
		if( is_v( gate ) )
	        v_to_clifford(circ_in, circ_out, gate);
		else
			circ_out.append_gate() = gate;
	}
}

void try_all_permutations(circuit circ)
{
    int perm[5] = {0, 1, 2, 3, 4}, inv_perm[5], best_perm[5] = {0, 1, 2, 3, 4};
    circuit circ_test;
    do
    {
    	clear_circuit( circ_test );
    	copy_circuit(circ, circ_test);
	    permute_lines( circ_test , perm );
        circ_test = remove_dup_gates( circ_test );
        for( int i = 0; i < 5; i++ )
            std::cout << perm[i] << " ";
        std::cout << "gates = " << circ_test.num_gates() << std::endl;
   	} while ( std::next_permutation(perm,perm+5) );
}

void single_toffoli_different_v_gates(circuit circ_in)
{
	circuit aux;
    std::cout << "Original" << std::endl;
    std::cout << circ_in << std::endl;
	if( circ_in.num_gates() == 1 )
	{
		for ( const auto& gate : circ_in )
	    {
			if ( is_toffoli( gate ) )
	        {
	            if( gate.controls().size() == 2 )
	            {
	            	clear_circuit(aux);
	            	copy_metadata(circ_in, aux);
	          	 	toffoli_to_v(aux, gate);
	          	 	std::cout << aux << std::endl;
	          	 	aux = transform_to_IBMQ( aux, map_method_qx4, true );
	          	 	if(remove)
	          	 		aux = remove_dup_gates( aux );
	          	 	std::cout << aux << std::endl;
	          	 	std::cout << aux.num_gates() << std::endl;
	          	 	
    		         	 	
    				std::cout << "\nNext" << std::endl;
	          	 	clear_circuit(aux);
	            	copy_metadata(circ_in, aux);
	          	 	interchange = true;
	          	 	toffoli_to_v(aux, gate);
	          	 	std::cout << aux << std::endl;
	          	 	aux = transform_to_IBMQ( aux, map_method_qx4, true );
	          	 	if(remove)
	          	 		aux = remove_dup_gates( aux );
	          	 	std::cout << aux << std::endl;
	          	 	std::cout << aux.num_gates() << std::endl;
	            }
	            else
	            	std::cout << "Not a toffoli" << std::endl;
			}
			else
	           	std::cout << "Not a toffoli" << std::endl;
		}	
	}	
	
}

bool alex_command::execute()
{
	matrix qx4 ={{0,4,4,7,7},{0,0,4,7,7},{0,0,0,4,4},{3,3,0,0,0},{3,3,0,4,0}};
	is_set("interchange") ? interchange = true : interchange = false;
	is_set("direct") ? direct = true : direct = false;
	auto& circuits = env->store<circuit>();
	circuit circ = circuits.current();
	// std::cout << "	" << circ.num_gates() << std::endl;

 	if( is_set("single_Toffoli") )
 	{
 		if( is_set("remove") )
 			remove = true;
 		else
 			remove = false;
		single_toffoli_different_v_gates(circ);
		return true;
 	}

 	unsigned target, control;
    std::vector<unsigned int> controls_a, controls_b;

    circuit vCircuit;
   	circuit cliffordCircuit;

   	auto settings = make_settings();
 	settings->set( "controls_threshold", 2u );
 	circ = nct_mapping( circ, settings, statistics );
   	transform_to_v(circ, vCircuit, qx4);

	if ( is_set( "clifford" ) )
    	transform_to_clifford(vCircuit, cliffordCircuit);

	circuits.extend();
	if ( is_set( "clifford" ) )
	{
		if ( is_set( "permutations" ) )
			try_all_permutations(cliffordCircuit);
		circuits.current() = cliffordCircuit;
	}
	else
	{
		if ( is_set( "permutations" ) )
			try_all_permutations(vCircuit);
		circuits.current() = vCircuit;
	}
	
	return true;
}

command::log_opt_t alex_command::log() const
{
  return log_opt_t({
	  {"runtime",       statistics->get<double>( "runtime" )}
	});
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
