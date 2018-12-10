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

#include <boost/program_options.hpp>

#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include <cli/reversible_stores.hpp>
#include <reversible/utils/matrix_utils.hpp>
#include <alice/rules.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/add_line_to_circuit.hpp>
#include <reversible/circuit.hpp>
#include <reversible/gate.hpp>
#include <reversible/pauli_tags.hpp>
#include <reversible/rotation_tags.hpp>
#include <reversible/target_tags.hpp>

namespace cirkit
{

using boost::program_options::value;

alex_command::alex_command( const environment::ptr& env )
	: cirkit_command( env, "Random projects" )
{
	opts.add_options()
    ( "interchange,c",  "interchange controls" )
    ( "clifford,t",  "transform to Clifford+T" )
	;

}

bool interchange = false;

command::rules_t alex_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

void toffoli_to_v(circuit& circ_in, circuit& circ_out, gate g)
{
    std::vector<unsigned int> controls_a, controls_b;
    std::vector<unsigned int> polarities;
    polarities.clear();
    controls_a.clear();
	controls_b.clear();
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

void multiControlToffoli( circuit& circ_in, gate g )
{
	circuit c;
	copy_metadata(circ_in, c);
	// append_toffoli( circ_out, g.controls(), g.targets().front() );

}

void transform_to_toffoli(circuit& circ_in, circuit& circ_out)
{
	copy_metadata(circ_in, circ_out);
	for (int i = circ_in.lines() - 1; i > 3; --i)
	{
		for ( const auto& gate : circ_in )
    	{
			if ( is_toffoli( gate ) )
	        {
	            if( gate.controls().size() == circ_in.lines() - 1 )
					add_line_to_circuit( circ_out, "i" + boost::lexical_cast<std::string>(circ_in.lines()) , "o" + boost::lexical_cast<std::string>(circ_in.lines()));
	            else if( gate.controls().size() > 2 )
	          	 	append_circuit( multiControlToffoli(), circ_out );
	            else
	            	circ_out.append_gate() = gate;
			}
			else
				circ_out.append_gate() = gate;
		}
	}
}

void transform_to_v(circuit& circ_in, circuit& circ_out)
{
	copy_metadata(circ_in, circ_out);
	for ( const auto& gate : circ_in )
    {
		if ( is_toffoli( gate ) )
        {
            if( gate.controls().size() == 2 )
          	 	toffoli_to_v(circ_in, circ_out, gate);
            else
            {
            	multiControlToffoli(gate, circ_in.lines());
            	circ_out.append_gate() = gate;
            }
		}
		else
			circ_out.append_gate() = gate;
	}
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

bool alex_command::execute()
{
	if ( is_set("interchange") )
		interchange = true;
	else
		interchange = false;

	auto& circuits = env->store<circuit>();
	circuit circ = circuits.current();
	// std::cout << "	" << circ.num_gates() << std::endl;
 	
 	unsigned target, control;
    std::vector<unsigned int> controls_a, controls_b;

    circuit vCircuit;
   	circuit cliffordCircuit;
    
    transform_to_v(circ, vCircuit);
    transform_to_toffoli(circ, vCircuit);

	if ( is_set( "clifford" ) )
    	transform_to_clifford(vCircuit, cliffordCircuit);
	
	circuits.extend();
	if ( is_set( "clifford" ) )
		circuits.current() = cliffordCircuit;
	else
		circuits.current() = vCircuit;
	
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
