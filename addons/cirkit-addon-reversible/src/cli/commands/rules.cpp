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

#include "rules.hpp"

#include <iostream>

#include <reversible/circuit.hpp>


namespace cirkit
{


//Remove control of the line with target
void remove_line_control_target( gate& control, const gate& target )
{
	for( auto& v : control.controls() )
		if( v.line() == target.targets().at(0) ) 
			control.remove_control ( v ); 

}

//Control is in the same line of the target?
bool line_control_target( const gate& control, const gate& target )
{
	for( auto& v : control.controls() )
		if( v.line() == target.targets().at(0) )
			return true;
	return false;
}

//Different polarity in the same line?
bool different_polarity_controls( const gate& ga, const gate& gb )
{
	for( const auto& v : ga.controls() )
		for( const auto& z : gb.controls() )
			if( v.line() == z.line() && v.polarity() != z.polarity() )
				return true;
	return false;
}

//Only one control different comparing two adjacents gates
bool single_control( const gate& ga, const gate& gb )
{
	unsigned i = 0;
	for( const auto& v : ga.controls() )
	{
		bool nex = true;
		for( const auto& z : gb.controls() )
			if( v.line() == z.line() )
				nex = false;
		if( nex )
			++i;	
	}
	
	for( const auto& v : gb.controls() )
	{
		bool nex = true;
		for( const auto& z : ga.controls() )
			if( v.line() == z.line() )
				nex = false;
		if( nex )
			++i;	
	}
	
	if( i == 1 )
		return true;
	return false;
}

//Controls in the same line?
bool controls_same_line( const gate& ga, const gate& gb )
{
	for( const auto& v : ga.controls() )
		for( const auto& z : gb.controls() )
			if( v.line() == z.line() )
				return true;
	return false;
}

//Targets in the same line?
bool targets_same_line( const gate& ga, const gate& gb )
{
	if( ga.targets() == gb.targets() )
		return true;
	return false;
}

//Change the order of two gates
void swap_gates( circuit::const_iterator itGate, circuit::const_iterator nextGate )
{
	gate g;
	g.controls() = itGate->controls();
	g.targets() = itGate->targets();
	itGate->controls() = nextGate->controls();
	itGate->targets() = nextGate->targets();
	nextGate->controls() = g.controls();
	nextGate->targets() = g.targets();
}

void apply_rule_Rfive( gate& ga, gate& gb )
{	
	for( auto v = ga.controls().begin();  v != ga.controls().end(); ++v )
		for( auto z = gb.controls().begin();  z != gb.controls().end(); ++z )
			if( v->line() == z->line() && v->polarity() == z->polarity() )
			{
				ga.remove_control ( ga.controls().at( v - ga.controls().begin() ) );
				gb.remove_control ( gb.controls().at( z - gb.controls().begin() ) );
				--v; --z;
			}
}

//Verifying if two adjacents gates are equal
bool verify_rule_Done( const gate& ga, const gate& gb )
{
	if( ga.controls() == gb.controls() && ga.targets() == gb.targets() )
		return true;
	return false;
}

//Remove equal adjacent gates
void apply_rule_Done( circuit& circ, unsigned itGateIndex, unsigned nextGateIndex )
{
	circ.remove_gate_at( nextGateIndex );
	circ.remove_gate_at( itGateIndex );
}

//Control rule
bool verify_rule_Dthree( const gate& ga, const gate& gb )
{
	if( targets_same_line( ga, gb ) )
		if( different_polarity_controls ( ga, gb ) )
			if( ga.controls().size() == 1 && gb.controls().size() == 1 )
				return true;
	return false;
}

void apply_rule_Dthree( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate )
{
	for( const auto& v : itGate->controls() )
		for( const auto& z : nextGate->controls() )
			if( v.line() == z.line() && v.polarity() != z.polarity() )
				itGate->remove_control( v );
	circ.remove_gate_at( nextGate - circ.begin() );
}

//Merge rule
bool verify_rule_Dfour( const gate& ga, const gate& gb )
{	
	if( targets_same_line( ga, gb ) && !different_polarity_controls ( ga, gb ) )
		if( single_control( ga, gb ) )
				return true;	
	return false;
}

void apply_rule_Dfour( circuit& circ, circuit::const_iterator& itGate, circuit::const_iterator& nextGate )
{
	/*for( auto v = itGate->controls().begin();  v != itGate->controls().end(); ++v )
		std::cout << "i: " << v->line() << v->polarity() << std::endl;
	for( auto z = nextGate->controls().begin();  z != nextGate->controls().end(); ++z )
		std::cout << "n: " << z->line() << z->polarity() << std::endl;*/
	if( itGate->controls().size() < nextGate->controls().size() )
	{	
		for( auto& z : nextGate->controls() )
		{
			bool ex = true;
			for( const auto& v : itGate->controls() )
			{
				if( z == v )
					ex = false;
			}
			if( ex )
			{
				z.set_polarity( !z.polarity() );
				break;
			}
		}			
		circ.remove_gate_at( itGate - circ.begin() );
	}
	else
	{
		for( auto& v : itGate->controls() )
		{
			bool ex = true;
			for( const auto& z : nextGate->controls() )
			{
				if( v == z )
					ex = false;
			}
			if( ex )
			{
				v.set_polarity( !v.polarity() );
				break;
			}
		}
		circ.remove_gate_at( nextGate - circ.begin() );
	}
}

/*****************************************************************
			MOVING RULES
*****************************************************************/

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

//Moving rule ( controls with different polarities )
bool verify_rule_Rfour ( const gate& ga, const gate& gb )
{
	if ( different_polarity_controls ( ga, gb ) )
		if ( !targets_same_line ( ga, gb ) )
		return true;
	return false;
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

//Insert control rule
bool verify_rule_Dfive( const gate& ga, const gate& gb )
{	
	if( targets_same_line( ga, gb ) )
		if( ga.controls().size() == 1 && gb.controls().size() == 1 )
			if( !controls_same_line( ga, gb ) )
				return true;
	return false;

}

//inseting because of tabu ... I'll take a closer look 
void apply_rule_Dfive( circuit::const_iterator& itGate, circuit::const_iterator& nextGate )
{
	unsigned l;
	for( auto& v : itGate->controls() )
		if( std::find(nextGate->controls().begin(), nextGate->controls().end(), v ) == nextGate->controls().end() )
		{
			l = v.line();
			v.set_polarity( !v.polarity() ); 
			nextGate->add_control( v );
			v.set_polarity( !v.polarity() );
		}
	for( auto& z : nextGate->controls() )
		if( std::find(itGate->controls().begin(), itGate->controls().end(), z ) == itGate->controls().end() )		
			if( z.line() != l )
			{
				//std::cout << z.line() << l << std::endl;
				z.set_polarity( !z.polarity() ); 
				itGate->add_control( z );
				z.set_polarity( !z.polarity() );
			}
}

//Remove control rule
bool verify_rule_Dfivee( const gate& ga, const gate& gb )
{	
	unsigned l = 0;
	if( targets_same_line( ga, gb ) )
		if( ga.controls().size() == 2 && gb.controls().size() == 2 )
			for( const auto& v : ga.controls() )
				for( const auto& z : gb.controls() )
					if( v.line() == z.line() && v.polarity() != z.polarity() )
						++l;
	if( l == 2 )
		return true;
	else
		return false;
	return false;
}

void apply_rule_Dfivee( circuit::const_iterator& itGate, circuit::const_iterator& nextGate )
{
	bool l = false;
	for( auto& v : itGate->controls() )
		for( auto& z : nextGate->controls() )
			if( v.line() == z.line() && v.polarity() != z.polarity() )
			{
				if( l )
					itGate->remove_control( v );
				else
				{
					nextGate->remove_control( z );
					l = true;
				}
			}	
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
