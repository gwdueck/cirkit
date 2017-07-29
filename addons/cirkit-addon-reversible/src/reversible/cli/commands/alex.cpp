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
#include <reversible/functions/add_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/is_identity.hpp>
#include <reversible/functions/reverse_circuit.hpp>


namespace cirkit
{

using boost::program_options::value;

    alex_command::alex_command(const environment::ptr& env )
    : cirkit_command( env, "Alex test" )
    {
    }

//Remove control of the line with target
void remove_line_control_target ( gate& control, const gate& target )
{
	for ( auto& v : control.controls() )
		if ( v.line() == target.targets().at(0) ) 
			control.remove_control ( v ); 

}

//Control is in the same line of the target?
bool line_control_target ( const gate& control, const gate& target )
{
	for ( auto& v : control.controls() )
		if ( v.line() == target.targets().at(0) )
			return true;
	return false;
}

//Different polarity in the same line?
bool different_polarity_controls ( const gate& ga, const gate& gb )
{
	for ( const auto& v : ga.controls() )
		for ( const auto& z : gb.controls() )
			if ( v.line() == z.line() && v.polarity() != z.polarity() )
				return true;
	return false;
}

//Controls in the same line?
bool controls_same_line ( const gate& ga, const gate& gb )
{
	for ( const auto& v : ga.controls() )
		for ( const auto& z : gb.controls() )
			if ( v.line() == z.line() )
				return true;
	return false;
}

//Targets in the same line?
bool targets_same_line ( const gate& ga, const gate& gb )
{
	if ( ga.targets() == gb.targets() )
		return true;
	return false;
}

//Change the order of two gates
void swap_gates ( circuit::const_iterator itGate, circuit::const_iterator nextGate )
{
	gate g;
	g.controls() = itGate->controls();
	g.targets() = itGate->targets();
	itGate->controls() = nextGate->controls();
	itGate->targets() = nextGate->targets();
	nextGate->controls() = g.controls();
	nextGate->targets() = g.targets();
}

//Verifying if two adjacents gates are equal
bool verify_rule_Done ( const gate& ga, const gate& gb )
{
	if ( ga.controls() == gb.controls() && ga.targets() == gb.targets() )
		return true;
	return false;
}

//Remove equal adjacent gates
void apply_rule_Done ( circuit& circ, unsigned itGateIndex, unsigned nextGateIndex )
{
	circ.remove_gate_at ( nextGateIndex );
	circ.remove_gate_at ( itGateIndex );
}

//Moving rule ( control positive -> control negative )
bool verify_rule_Dtwo (  gate ga, gate gb )
{
	if ( line_control_target ( ga, gb  ) && !line_control_target ( gb , ga ) )
	{	
		remove_line_control_target ( ga, gb );
		if ( ga.controls() == gb.controls() )
			return true;
	}
	if ( line_control_target ( gb, ga  ) && !line_control_target ( ga , gb ) )
	{	
		remove_line_control_target ( gb, ga );
		if ( ga.controls() == gb.controls() )
			return true;
	}		
	return false;
}

void apply_rule_Dtwo ( circuit::const_iterator itGate, circuit::const_iterator nextGate )
{
	swap_gates ( itGate, nextGate );
	
	for ( auto& v : itGate->controls() )
		if ( v.line() == nextGate->targets().at(0) )
		{
			v.set_polarity ( !v.polarity() ); 
			return;
		}
	
	for ( auto& v : nextGate->controls() )
		if ( v.line() == itGate->targets().at(0) )
			v.set_polarity ( !v.polarity() ); 
}

//Control rule
bool verify_rule_Dthree ( const gate& ga, const gate& gb )
{
	if ( targets_same_line ( ga, gb ) )
		if ( different_polarity_controls ( ga, gb ) )
			if ( ga.controls().size() == 1 && gb.controls().size() == 1 )
				return true;
	return false;
}

void apply_rule_Dthree ( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate )
{
	for ( const auto& v : itGate->controls() )
		for ( const auto& z : nextGate->controls() )
			if ( v.line() == z.line() && v.polarity() != z.polarity() )
				itGate->remove_control ( v );
	circ.remove_gate_at ( nextGate - circ.begin() );
}

//Merge rule
bool verify_rule_Dfour ( const gate& ga, const gate& gb )
{
	gate gc, gd;
	if ( targets_same_line ( ga, gb ) && !different_polarity_controls ( ga, gb ) )
	{
		std::set_difference( ga.controls().begin(), ga.controls().end(), gb.controls().begin(), gb.controls().end(), std::inserter( gc.controls() , gc.controls().end() ) );
		if ( gc.controls().size() == 1 )
			return true;
		std::set_difference( gb.controls().begin(), gb.controls().end(), ga.controls().begin(), ga.controls().end(), std::inserter( gd.controls() , gd.controls().end() ) );
		if ( gd.controls().size() == 1 )
			return true;	
	}
	return false;
}

void apply_rule_Dfour ( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate )
{
	if ( itGate->controls().size() < nextGate->controls().size() )
	{	
		if ( itGate->controls().size() == 0 )
			nextGate->controls().at(0).set_polarity ( !nextGate->controls().at(0).polarity() );
		else
			for ( auto v = itGate->controls().begin(), z = nextGate->controls().begin();  z != nextGate->controls().end(); ++v, ++z )
				if ( v->line() != z->line() )
				{
					z->set_polarity ( !z->polarity() );
					break;
				}
		circ.remove_gate_at ( itGate - circ.begin() );
	}
	else
	{
		if ( nextGate->controls().size() == 0 )
			itGate->controls().at(0).set_polarity ( !itGate->controls().at(0).polarity() );
		else
			for ( auto v = itGate->controls().begin(), z = nextGate->controls().begin();  v != itGate->controls().end(); ++v, ++z )
				if ( v->line() != z->line() )
				{
					v->set_polarity ( !v->polarity() );
					break;
				}
		circ.remove_gate_at ( nextGate - circ.begin() );
	}
}

//Moving rule ( controls with different polarities )
bool verify_rule_Rfour ( const gate& ga, const gate& gb )
{
	if ( different_polarity_controls ( ga, gb ) )
		if ( !targets_same_line ( ga, gb ) )
		return true;
	return false;
}

void apply_rule_Rfive ( gate& ga, gate& gb )
{	
	for ( auto v = ga.controls().begin();  v != ga.controls().end(); ++v )
		for ( auto z = gb.controls().begin();  z != gb.controls().end(); ++z )
			if ( v->line() == z->line() && v->polarity() == z->polarity() )
			{
				ga.remove_control ( ga.controls().at( v - ga.controls().begin() ) );
				gb.remove_control ( gb.controls().at( z - gb.controls().begin() ) );
				--v; --z;
			}
}

//Moving rule
bool verify_rule_Dsix ( const gate& ga, const gate& gb )
{
	if ( !line_control_target ( ga, gb ) && !line_control_target ( gb, ga ) )
		if ( !targets_same_line ( ga, gb ) ) 
			return true;
	
	if ( targets_same_line ( ga, gb ) )
		return true;
		
	return false;
}

//Moving rule with ++cost
bool verify_rule_Dseven ( const gate& ga, const gate& gb ) 
{
	if ( line_control_target ( ga, gb ) && !line_control_target ( gb, ga ) )
		if ( !controls_same_line ( ga, gb ) )
			return true;
	
	if ( line_control_target ( gb, ga ) && !line_control_target ( ga, gb ) )
		if ( !controls_same_line ( ga, gb ) )
			return true;
	
	return false;
}

void apply_rule_Dseven ( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate )
{	
	unsigned n = ( nextGate - circ.begin() ) + 1;
	unsigned target;
	bool it = true;
	gate g;
	
	for ( const auto& v : itGate->controls() )
		if ( v.line() == nextGate->targets().at( 0 ) ) 
			it = false;
		else
			g.add_control ( v );
				
	for ( const auto& v : nextGate->controls() )
		if ( v.line() != itGate->targets().at( 0 ) ) 
			g.add_control ( v );
	
	if ( it )
		target = nextGate->targets().at( 0 );
	else
		target = itGate->targets().at( 0 );
		
	std::sort ( g.controls().begin(), g.controls().end() );
	auto i = std::unique ( g.controls().begin(), g.controls().end() );
	g.controls().erase ( i, g.controls().end() );
	
	swap_gates ( itGate, nextGate );
	insert_toffoli ( circ, n, g.controls(), target );
}

void applying_rules ( circuit& circ, unsigned x[][3], unsigned y )
{
	circuit::const_iterator itGate = circ.begin(), nextGate = circ.begin();
	switch ( x[y][0] )
	{
		case 1:
			apply_rule_Done ( circ, x[y][1], x[y][2] );
			break;
		case 2:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dtwo ( itGate, nextGate );
			break;
		case 3:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dthree ( circ, itGate, nextGate );
			break;
		case 4:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			swap_gates ( itGate, nextGate );
			break;
		case 6:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			swap_gates ( itGate, nextGate );
			break;
		case 7:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dseven ( circ, itGate, nextGate );
			break;
		case 14:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dfour ( circ, itGate, nextGate );
			break;
		default:
			std::cout << "Something got wrong!" << std::endl;
			break;
	}
}

//Circuit is the same truth table?
void circuit_is_the_same ( circuit circ, circuit orig )
{
	reverse_circuit ( orig );
	append_circuit ( circ, orig );
	if ( !is_identity ( circ ) )
		std::cout << "Some rule changed the truth table of the circuit!" << std::endl;
}

//Choose the rule to be applied
void choose_rule ( unsigned& y, const unsigned i )
{
	while ( true )
	{
		std::cout << "\nWich rule will be applied? 0 to exit. ";
		std::cin >> y;
		
		if ( std::cin.fail() || y < 0 || y > i )
		{
			std::cin.clear();
			std::cin.ignore();
			std::cout << "Choose between 1 and " << i << ". 0 to exit. ";
		}
		else
			break;
	}
}

void testing ( circuit& circ )
{
	gate ga, gb;
	unsigned y = 0, i = 0;
	unsigned x[100][3];
	
	while ( true )
	{ 
		i = 0;
		
		//Printing the circuit
		std::cout << std::endl << circ << std::endl;
		
		std::cout << "Number of gates: " << circ.num_gates() << std::endl;
		
		std::cout << "Iterating through the circuit..." << std::endl;
		
		for ( circuit::const_iterator itGate = circ.begin(), nextGate = ++circ.begin();  nextGate != circ.end(); ++itGate, ++nextGate )
		{
			unsigned itGateIndex = itGate - circ.begin();
			unsigned nextGateIndex = nextGate - circ.begin();
			
			std::sort ( itGate->controls().begin(), itGate->controls().end() );
			std::sort ( nextGate->controls().begin(), nextGate->controls().end() );
			
			ga.controls() = itGate->controls();
			gb.controls() = nextGate->controls();
			ga.targets() = itGate->targets();
			gb.targets() = nextGate->targets();
						
			if ( verify_rule_Done ( ga, gb ) )
			{
				++i;
				std::cout << i << ". Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be removed.\t\tCost=-2" << std::endl;
				x[i][0] = 1;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
			}
		
			apply_rule_Rfive ( ga, gb );

			if ( verify_rule_Dtwo ( ga, gb ) ) 
			{
				++i;
				std::cout << i << ". Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
				x[i][0] = 2;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
			}

			if ( verify_rule_Dthree ( ga, gb ) )
			{
				++i;
				std::cout << i << ". Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be merged.\t\t--Cost" << std::endl;
				x[i][0] = 3;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
			}

			if ( verify_rule_Dfour ( ga, gb ) )
			{
				++i;
				std::cout << i << ". Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be merged.\t\t--Cost" << std::endl;
				x[i][0] = 14;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
			}

			if ( verify_rule_Rfour ( ga, gb ) )
			{
				++i;
				std::cout << i << ". Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
				x[i][0] = 4;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
			}

			if ( verify_rule_Dsix ( ga, gb ) )
			{
				++i;
				std::cout << i << ". Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
				x[i][0] = 6;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
				//swap_gates ( itGate, nextGate );  
			}

			if ( verify_rule_Dseven ( ga, gb ) )
			{
				++i;
				std::cout << i << ". Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\t++Cost" << std::endl;
				x[i][0] = 7;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
			}

		}
		
		if ( i == 0 )
		{
			std::cout << "\nNo rule to be applied!" << std::endl;
			break;
		}
		else
		{
			choose_rule ( y, i);
			if ( y == 0 )
				break;
			applying_rules ( circ, x, y );
		}
	}
}

bool alex_command::execute()
{
	auto& circuits = env->store<circuit>();
	
	if( env->store<circuit>().current_index() >= 0 )
	{
	 	circuit circ, orig;
	 	circ = circuits.current();
	 	copy_circuit ( circ, orig );
	 	testing ( circ );
	 	circuit_is_the_same ( circ, orig );
	 	circuits.current() = circ;
	}
	else
	{
		std::cout << "no circuit in store" << std::endl;
	}
	return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
