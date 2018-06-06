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

//#include <reversible/cli/stores.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/circuit.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/functions/add_circuit.hpp>
#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/copy_metadata.hpp>
#include <reversible/functions/clear_circuit.hpp>
#include <reversible/functions/is_identity.hpp>
#include <reversible/functions/reverse_circuit.hpp>

#include <cli/commands/rules.hpp>
#include <cli/commands/perm.hpp>

#include <core/utils/program_options.hpp>
#include <reversible/utils/permutation.hpp>

using boost::program_options::value;

bool selectRandom = false;
//unsigned iteration;

namespace cirkit
{

 program_options opts;
 
alex_command::alex_command( const environment::ptr& env )
    : cirkit_command( env, "Alex test" )
{
	opts.add_options()
    ( "random,r",          					"Select random rules automatically" )
    //( "iterations,t", value(&iteration),	"Select number of iterations" )
    ;

}

//Circuit is the same truth table?
void circuit_is_the_same( circuit circ, circuit orig )
{
	append_circuit( circ, orig );
	if( !is_identity ( circ ) )
	{
		std::cout << "Some rule changed the truth table of the circuit!" << std::endl;
		exit( 1 );
	}
}

void applying_rules( circuit& circ, std::vector<std::vector<int>> x, unsigned y )
{
	circuit::const_iterator itGate = circ.begin(), nextGate = circ.begin();
	switch( x[y][0] )
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

//Choose the rule to be applied
void choose_rule( unsigned& y, const unsigned i )
{
	if ( selectRandom )
	{
		y = 1 + rand() % i;
		std::cout << "Sorteio: " << y << std::endl;
	}
	else
	{
		while( true )
		{
			std::cout << "\nWich rule will be applied? 0 to exit. ";
			std::cin >> y;
		
			if( std::cin.fail() /*|| y < 0 */ || y > i )
			{
				std::cin.clear();
				std::cin.ignore();
				std::cout << "Choose between 1 and " << i << ". 0 to exit. ";
			}
			else
				break;
		}
	}
}

void testing( circuit& circ, const circuit orig )
{
	gate ga, gb;
	unsigned y = 0, i = 0;
	//unsigned x[10000][4];
	unsigned stop = 0;
	unsigned start, min, max;
	circuit mc;

	std::vector <std::vector<int>> x(10000, std::vector<int>(4));
	
	start = min = max = circ.num_gates();
	copy_circuit( circ, mc );
 	//while( true )
	while( stop < 1000 )
	{ 
		i = 0;
		
		if( circ.num_gates() < min )
		{
			min = circ.num_gates();
			clear_circuit( mc );
			copy_circuit( circ, mc );
		}	
		if( circ.num_gates() > max )
			max = circ.num_gates();
		if( circ.num_gates() > min*1.5)
		{
			clear_circuit( circ );	
			copy_circuit( mc, circ );
			
		}
		
		//Printing the circuit
		std::cout << std::endl << circ << std::endl;
		circuit_is_the_same( circ, orig );
		std::cout << "Number of gates: " << circ.num_gates() << std::endl;
		std::cout << "Iteration: " << stop << std::endl;
		std::cout << "Iterating through the circuit..." << std::endl;
		
		for( circuit::const_iterator itGate = circ.begin(), nextGate = ++circ.begin();  nextGate != circ.end(); ++itGate, ++nextGate )
		{
			unsigned itGateIndex = itGate - circ.begin();
			unsigned nextGateIndex = nextGate - circ.begin();
			
			std::sort ( itGate->controls().begin(), itGate->controls().end() );
			std::sort ( nextGate->controls().begin(), nextGate->controls().end() );
			
			ga.controls() = itGate->controls();
			gb.controls() = nextGate->controls();
			ga.targets() = itGate->targets();
			gb.targets() = nextGate->targets();
						
			if( verify_rule_Done ( ga, gb ) )
			{
				++i;
				std::cout << i << ". [D1] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be removed.\t\tCost=-2" << std::endl;
				x[i][0] = 1;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
				x[i][3] = -2;
			}
		
			apply_rule_Rfive ( ga, gb );

			if( verify_rule_Dtwo ( ga, gb ) ) 
			{
				++i;
				std::cout << i << ". [R5] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
				x[i][0] = 2;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
				x[i][3] = 0;
			}

			if( verify_rule_Dthree ( ga, gb ) )
			{
				++i;
				std::cout << i << ". [D3] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be merged.\t\t--Cost" << std::endl;
				x[i][0] = 3;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
				x[i][3] = -1;
			}

			if( verify_rule_Dfour ( ga, gb ) )
			{
				++i;
				std::cout << i << ". [D4] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be merged.\t\t--Cost" << std::endl;
				x[i][0] = 14;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
				x[i][3] = -1;
			}

			if( verify_rule_Rfour ( ga, gb ) )
			{
				++i;
				std::cout << i << ". [R4] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
				x[i][0] = 4;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
				x[i][3] = 0;
			}

			if( verify_rule_Dsix ( ga, gb ) )
			{
				++i;
				std::cout << i << ". [D6] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
				x[i][0] = 6;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;  
				x[i][3] = 0;
			}

			if( verify_rule_Dseven ( ga, gb ) )
			{
				++i;
				std::cout << i << ". [D7] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\t++Cost" << std::endl;
				x[i][0] = 7;
				x[i][1] = itGateIndex;
				x[i][2] = nextGateIndex;
				x[i][3] = 1;
			}

		}
		
		if( i == 0 )
		{
			std::cout << "\nNo rule to be applied!" << std::endl;
			break;
		}
		else
		{
			//x.erase( x.begin()+10 );
			choose_rule( y, i );
			if( y == 0 )
				break;
			applying_rules( circ, x, y );
			++stop;
		}
	}
	if( selectRandom )
		std::cout << "Number of gates-> Start: " << start << " Min: " << min << " Max: " << max << std::endl;
}

void print_perm(std::vector<unsigned>& ppp)
{
	for (int i = 0; i < ppp.size(); ++i)
		std::cout << " " << ppp[i];
	std::cout << std::endl;
}

std::string vector_to_string(std::vector<unsigned>& ppp)
{
	std::string p;
	for (int i = 0; i < ppp.size(); ++i)
		p.append(std::to_string(ppp[i]));
	return p;
}

bool alex_command::execute()
{
	auto& circuits = env->store<circuit>();
	circuit circ, aux;
	circ = circuits.current();
	
	std::vector<unsigned> ppp;
	std::vector<std::string> m;
	std::string p;

	copy_metadata(circ, aux);
	// ppp = circuit_to_permutation(circ);
	// print_perm(ppp);

	unsigned int y = 0;
	unsigned int total = circ.num_gates();
	for(auto g : circ)
	{
		++y;
		std::cout << "\r" << (y*100)/total << "%";
		std::cout.flush();
		aux.append_gate() = g;
		ppp = circuit_to_permutation(aux);
		p = vector_to_string(ppp);
		auto x = std::find(m.begin(), m.end(), p); 
		if(x != m.end())
			std::cout << "\nEncontrou igual! " << y << " " << (x-m.begin())+1 << std::endl;
		m.push_back(p);
	}
	std::cout << std::endl;

	//imprimir o vetor com as permutacoes
	// for (int i = 0; i < m.size(); ++i)
	// 	std::cout << m[i] << std::endl;


	// if ( is_set( "random" ) )
	// 	selectRandom = true;
	// if( env->store<circuit>().current_index() >= 0 )
	// {
	//  	circuit circ, orig;
	//  	circ = circuits.current();
	//  	copy_circuit( circ, orig );
	//  	reverse_circuit( orig );
	//  	testing( circ, orig );
	//  	circuit_is_the_same( circ, orig );
	//  	circuits.current() = circ;
	// }
	// else
	// {
	// 	std::cout << "no circuit in store" << std::endl;
	// }
	// //selectRandom = false;
	return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
