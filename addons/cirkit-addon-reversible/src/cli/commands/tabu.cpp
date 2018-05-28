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

#include "tabu.hpp"

#include <iostream>
#include <time.h>
#include <core/utils/timer.hpp>
#include <alice/rules.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/circuit.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/functions/add_circuit.hpp>
//#include <reversible/functions/add_gates.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/clear_circuit.hpp>
#include <reversible/functions/is_identity.hpp>
#include <reversible/functions/reverse_circuit.hpp>

//#include <reversible/pauli_tags.hpp>
//#include <reversible/target_tags.hpp>

#include <cli/commands/rules.hpp>
#include <core/utils/program_options.hpp>


typedef std::vector<std::vector<int>> matrix;

namespace cirkit
{
 
tabu_command::tabu_command(const environment::ptr& env)
    : cirkit_command(env, "Tabu Search")
{
	opts.add_options()
    //( "random,r",          					"Select random rules automatically" )
    //( "i",value_with_default(&iteration),		"Select number of iterations" )
    //( "p",value_with_default(&penalization),	"Select number of penalization" )
    ( "n",value_with_default(&neighboorhod),	"Select the size of the neighboorhod" )
    ( "o",value_with_default(&overlap),			"Select the size of the overlap (%)" )
    ;
}

command::rules_t tabu_command::validity_rules() const
{
	return {has_store_element<circuit>( env )};
}

//Any rule changed the truth table?
// bool changed_circuit(circuit circ, circuit orig)
// {
// 	append_circuit(circ, orig);
// 	if(!is_identity(circ))
// 	{
// 		std::cout << "Some rule changed the truth table of the circuit!" << std::endl;
// 		return false;
// 	}
// 	return true;
// }

void list_rules( circuit circ, std::vector<std::vector<int>>& x, unsigned begin, unsigned end )
{
	for( circuit::const_iterator itGate = circ.begin() + begin, nextGate = ++( circ.begin() + begin );  nextGate != circ.begin() + end; ++itGate, ++nextGate )
	{
		unsigned itGateIndex = itGate - circ.begin();
		unsigned nextGateIndex = nextGate - circ.begin();
		gate ga, gb;
		
		std::sort( itGate->controls().begin(), itGate->controls().end() );
		std::sort( nextGate->controls().begin(), nextGate->controls().end() );
		
		ga.controls() = itGate->controls();
		gb.controls() = nextGate->controls();
		ga.targets() = itGate->targets();
		gb.targets() = nextGate->targets();
		
		if( verify_rule_Done( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[D1] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be removed.\t\tCost=-2" << std::endl;
			v.push_back(1);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(-2);
			x.push_back( v );
		}
	
		apply_rule_Rfive( ga, gb );

		if( verify_rule_Dtwo( ga, gb ) ) 
		{
			std::vector<int> v;
			//std::cout  << "[R5] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
			v.push_back(2);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			x.push_back( v );
		}

		if( verify_rule_Dthree( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[D3] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be merged.\t\t\t--Cost" << std::endl;
			v.push_back(3);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(-1);
			x.push_back( v );
		}
		
		if( verify_rule_Dfour( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[D4] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be merged.\t\t\t--Cost" << std::endl;
			v.push_back(14);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(-1);
			x.push_back( v );
		}

		if( verify_rule_Rfour( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[R4] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
			v.push_back(4);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			x.push_back( v );
		}
		
		if( verify_rule_Dfive( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[D5] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can insert controls.\t\tNo cost change" << std::endl;
			v.push_back(5);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			x.push_back( v );
		}
		
		if( verify_rule_Dfivee( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[D5] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can remove controls.\t\tNo cost change" << std::endl;
			v.push_back(-5);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			x.push_back( v );
		}

		if( verify_rule_Dsix( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[D6] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\tNo cost change" << std::endl;
			v.push_back(6);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			x.push_back( v );
		}

		if( verify_rule_Dseven( ga, gb ) )
		{
			std::vector<int> v;
			//std::cout  << "[D7] Gates ( " << itGateIndex + 1 << " - " << itGateIndex + 2 << " ) can be interchanged.\t\t++Cost" << std::endl;
			v.push_back(7);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(1);
			x.push_back( v );
		}
	}

}

void apply_rule( circuit& circ, std::vector<std::vector<int>> x, unsigned y )
{
	circuit::const_iterator itGate = circ.begin(), nextGate = circ.begin();
	switch( x[y][0] )
	{
		case 1:
			apply_rule_Done( circ, x[y][1], x[y][2] );
			break;
		case 2:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dtwo( itGate, nextGate );
			break;
		case 3:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dthree( circ, itGate, nextGate );
			break;
		case 4:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			swap_gates( itGate, nextGate );
			break;
		case 5:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dfive( itGate, nextGate );
			break;
		case -5:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dfivee( itGate, nextGate );
			break;
		case 6:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			swap_gates( itGate, nextGate );
			break;
		case 7:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dseven( circ, itGate, nextGate );
			break;
		case 14:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dfour( circ, itGate, nextGate );
			break;
		default:
			std::cout << "Something got wrong!" << std::endl;
			break;
	}
}

//Sort the vector of the rules
void sortrows(matrix& x, int col)
{    
	std::sort(x.begin(),
              x.end(),
              [col](const std::vector<int>& lhs, const std::vector<int>& rhs) {
                  return lhs[col] < rhs[col];
              });
}

//Print vector
void print_matrix(matrix& x)
{
	for(unsigned i = 0; i < x.size(); ++i)
		for(int j = 0; j < 4; ++j)
			std::cout << x[i][j] << ' ';
		std::cout << std::endl;	
}

//Function to decide which rule must be applied
void choose_rule(circuit& circ, matrix x, matrix& tp)
{
	bool insert;
	//std::cout << t.size() << std::endl;
	for(unsigned k = 0; k < x.size(); ++k)
	{
		if(x[k][3] < 0)
		{
			apply_rule(circ, x, k);
			break;
		}
		insert = true;
		for(unsigned i = 0; i < tp.size(); ++i)
		{
			if(x[k][0] == tp[i][0] && x[k][1] == tp[i][1] && x[k][2] == tp[i][2])
			{
				++tp[i][3];
				insert = false;
				break;
			}
		}
		if(insert)
		{
			x[k][3] = 0;
			tp.push_back( x[k] );
			apply_rule( circ, x, k );
			break;
		}
	}
}

//Update the tabu list (penalization)
void update_tabu_list(matrix& tp, int penalization)
{
	for(unsigned i = 0; i < tp.size(); ++i)
	{
		if(tp[i][3] >= penalization)
		{
			tp.erase(tp.begin()+i);
			break;
		}
	}
}

//Save the minimum circuit found
bool update_circuit(circuit& circ, circuit& min)
{
	if(circ.num_gates() > min.num_gates() + 3)
	{
		clear_circuit(circ);
		copy_circuit(min, circ);	
	}
	else if(circ.num_gates() < min.num_gates())
	{
		clear_circuit(min);
		copy_circuit(circ, min);
		//std::cout << " Minimo: " << min.num_gates() << std::endl;
		return true;
	}
	return false;
}

void tabu_search( circuit& circ, unsigned overlap, unsigned neighboorhod, const properties::ptr& statistics )
{
	properties_timer t( statistics );
	unsigned stop = 0, begin, end;
	bool finish = false;
	matrix x; //list of possible rules of the current iteration
	matrix tp; //tabu list with the penalization
	circuit orig, min;

	copy_circuit( circ, orig );
 	copy_circuit( circ, min );
 	reverse_circuit( orig );
 	begin = 0;
	end = neighboorhod;
	overlap = (overlap * neighboorhod) * 0.01;

	while(!finish)
	{ 	
		if(end > circ.num_gates())
		{
			end = circ.num_gates();
			if(end < neighboorhod)
			{
				begin = 0;
			}
			else
			{
				begin = end - neighboorhod;
			}
			finish = true;
		}
		else if(int(begin) < 0)
		{
			begin = 0;
		}
		//std::cout << "begin: " << begin << " end: " << end << std::endl;
	 	stop = 0;
	 	unsigned improvement = circ.num_gates();
	 	while(stop < neighboorhod * 50)
	 	{
	 		list_rules(circ, x, begin, end);
	 		sortrows(x, 3);
	 		choose_rule(circ, x, tp);
	 		//std::cout << "Imprimindo matrix tabu antes de atualizar" << std::endl;
	 		//print_matrix( tp );
	 		update_tabu_list(tp, neighboorhod);
	 		//std::cout << "Imprimindo matrix tabu depois de atualizar" << std::endl;
	 		//print_matrix( tp );
	 		//std::cout << "Imprimindo matrix x" << std::endl;
	 		//print_matrix( x );
	 		if(!update_circuit(circ, min))
	 		{
	 			++stop;
	 		}
	 		else
	 		{
	 			stop = 0;
	 		}
	 		x.clear();
	 		//std::cin.get();
	 		
	 		//++stop;
	 		//std::cout << stop << std::endl;
	 		if(end > circ.num_gates())
				end = circ.num_gates();
			if(circ.num_gates() == 0)
			{
				finish = true;
				break;
			}
	 	}
	 	clear_circuit(circ);
	 	copy_circuit(min, circ);
	 	improvement = improvement - min.num_gates();
	 	//std::cout << "improvement: " << improvement << std::endl;
	 	//std::cin.get();
	 	begin = begin + (neighboorhod - overlap) - improvement;
	 	end = begin + neighboorhod;
	}

		clear_circuit(circ);
	copy_circuit(min, circ);
}


bool tabu_command::execute()
{
	srand (time(NULL));
	clock_t Ticks[2];
  Ticks[0] = clock();
	unsigned initialgates;
	
	auto& circuits = env->store<circuit>();	
	circuit circ;
 	circ = circuits.current();
 	initialgates = circ.num_gates();
 	tabu_search( circ, overlap, neighboorhod, statistics );
 	
 	
 	//std::cout << circ << std::endl;

	
	
 	//std::cout << min << std::endl;
 	

 	circuits.current() = circ;
 	Ticks[1] = clock();
    double Time = (Ticks[1] - Ticks[0]) * 1.0 / CLOCKS_PER_SEC;
    std::cout << initialgates << "	" << circ.num_gates() << "	" << Time << std::endl;
    print_runtime();
	return true;
}

command::log_opt_t tabu_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
