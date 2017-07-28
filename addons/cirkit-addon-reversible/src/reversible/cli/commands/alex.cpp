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

#include <iostream>

#include <reversible/cli/stores.hpp>
#include <reversible/circuit.hpp>
#include <reversible/io/print_circuit.hpp>


namespace cirkit
{

using boost::program_options::value;

    alex_command::alex_command(const environment::ptr& env )
    : cirkit_command( env, "Alex test" )
    {
    }

                           
                           
bool alex_command::execute()
{
	auto& circuits = env->store<circuit>();
	
	if( env->store<circuit>().current_index() >= 0 )
	{
	 	circuit circ;
	 	circ = circuits.current();
	 	testing ( circ );
	}
	else
	{
		std::cout << "no circuit in store" << std::endl;
	}
	return true;
}

void alex_command::testing ( circuit& circ )
{
	std::cout << circ << std::endl;
	std::cout << "Number of gates: " << circ.num_gates() << std::endl;
	std::cout << "Number of lines: " << circ.lines() << std::endl;
	
	std::cout << std::endl;
	std::cout << "Removing gate 3..." << std::endl;
	circ.remove_gate_at( 2 );
	std::cout << circ << std::endl;
	std::cout << "Number of gates: " << circ.num_gates() << std::endl;
	
	std::cout << std::endl;
	std::cout << "Iterate through the circuit..." << std::endl;
	for ( circuit::const_iterator itGate = circ.begin(); itGate != circ.end(); ++itGate )
  	{
  		unsigned itGateIndex = itGate - circ.begin();
  		std::cout << "Gate: " << itGateIndex + 1 << "\t";
  		
  		for ( const auto& v : itGate->controls() )
		{
			std::cout << "\tControl: " << v.line();
		}
		
		for ( const auto& l : itGate->targets() )
		{
			std::cout << "\tTarget: " << l;
		}
		
		std::cout << std::endl;
  	}
  	
  	std::cout << std::endl;
  	std::cout << "Verifying if two adjacents gates are equal..." << std::endl;
	for ( circuit::const_iterator itGate = circ.begin(), nextGate = ++circ.begin(); nextGate != circ.end(); ++itGate, ++nextGate )
	{
		unsigned itGateIndex = itGate - circ.begin();
		unsigned nextGateIndex = nextGate - circ.begin();
		
		if ( itGate->controls() == nextGate->controls() && itGate->targets() == nextGate->targets() )
		{
			std::cout << itGateIndex + 1 << " - " << nextGateIndex + 1 << "\tEQUAL!" << std::endl;
		}
		else
		{
			std::cout << itGateIndex + 1 << " - " << nextGateIndex + 1 << "\tDIFFERENT!" << std::endl;
		}
	}
	
	std::cout << std::endl;
	std::cout << "Changing all the positive controls to negative..." << std::endl;
	for ( circuit::const_iterator itGate = circ.begin(); itGate != circ.end(); ++itGate )
	{
		for ( auto& v : itGate->controls() )
		{
			v.set_polarity( false );
		}
	}
	std::cout << circ << std::endl;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
