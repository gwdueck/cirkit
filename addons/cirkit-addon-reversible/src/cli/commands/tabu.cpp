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
    ( "neighborhood,g",value_with_default(&neighborhood),	"Select the size of the neighborhood" )
    ( "optimization,t",value_with_default(&opt),				"Select the optimization. 0: gates (default), 1: quantum cost, 2: average of both" )
    ( "overlap,o",value_with_default(&overlap),			"Select the size of the overlap (%)" )
    ( "verbose,v",								"be verbose")
    ( "step,s",									"be verbose and step-by-step")
    ;
    add_new_option();
}

command::rules_t tabu_command::validity_rules() const
{
	return {has_store_element<circuit>( env )};
}

//Any rule changed the truth table?
bool changed_circuit(circuit circ)
{
	circuit aux;
	reverse_circuit(circ, aux);
	append_circuit(circ, aux);
	if(!is_identity(circ))
	{
		std::cout << "Some rule changed the truth table of the circuit!" << std::endl;
		return false;
	}
	return true;
}

int ncv_cost(unsigned int control)
{
	//Assuming negative control with the same cost
	switch (control){
		case 0:
			return 1;
		case 1:
			return 1;
		case 2:
			return 5;
		case 3:
			return 20;
		case 4:
			return 50;
		default:
			return 40*(control-3);		
	}
}

int t_depth(unsigned int control)
{
	//Assuming negative control with the same cost
	switch (control){
		case 0:
			return 0;
		case 1:
			return 0;
		case 2:
			return 3;
		case 3:
			return 12;
		case 4:
			return 30;
		default:
			return 24*(control-3);		
	}
}

unsigned int circuit_ncv_cost(circuit& circ)
{
	unsigned int cost = 0;
	for (auto g : circ)
 		cost += ncv_cost(g.controls().size());
 	return cost;
}

void list_rules( circuit circ, std::vector<std::vector<int>>& x, unsigned begin, unsigned end, bool verbose )
{
	for( circuit::const_iterator itGate = circ.begin() + begin, nextGate = ++( circ.begin() + begin );  nextGate != circ.begin() + end; ++itGate, ++nextGate )
	{
		unsigned int itGateIndex = itGate - circ.begin();
		unsigned int nextGateIndex = nextGate - circ.begin();
		int gaNCVCost, gbNCVCost, gaTdepth, gbTdepth; //Quantum cost of the gates
		unsigned int gaControls, gbControls; //Controls of the gates

		gate ga, gb;
		
		std::sort( itGate->controls().begin(), itGate->controls().end() );
		std::sort( nextGate->controls().begin(), nextGate->controls().end() );
		
		ga.controls() = itGate->controls();
		gb.controls() = nextGate->controls();
		ga.targets() = itGate->targets();
		gb.targets() = nextGate->targets();

		gaControls = ga.controls().size();
		gbControls = gb.controls().size();

		gaNCVCost = ncv_cost(ga.controls().size());
		gbNCVCost = ncv_cost(gb.controls().size());
		gaTdepth = t_depth(ga.controls().size());
		gbTdepth = t_depth(gb.controls().size());

		if( verify_rule_Done( ga, gb ) )
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[D1] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be removed.\t\tCost=-2;\t\tQCost:" << -1*(gaNCVCost + gbNCVCost) << std::endl;
			v.push_back(1);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(-2);
			v.push_back(-1*(gaNCVCost + gbNCVCost));
			x.push_back( v );
		}
	
		apply_rule_Rfive( ga, gb );

		if( verify_rule_Dtwo( ga, gb ) ) 
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[R5] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be interchanged.\t\tNo cost change;" << std::endl;
			v.push_back(2);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			v.push_back(0);
			x.push_back( v );
		}

		if( verify_rule_Dthree( ga, gb ) )
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[D3] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be merged.\t\t\t--Cost;\t\tQCost:" << -1*(gaNCVCost-1) << std::endl;
			v.push_back(3);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(-1);
			v.push_back(-1*(gaNCVCost-1));
			x.push_back( v );
		}
		
		if( verify_rule_Dfour( ga, gb ) )
		{
			std::vector<int> v;
			// if(verbose)
			// 	std::cout  << "[D4] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be merged.\t\t\t--Cost" << std::endl;
			v.push_back(14);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(-1);
			if(gaControls > gbControls)
			{
				v.push_back(-1*gaNCVCost);
				if(verbose)
					std::cout  << "[D4] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be merged.\t\t\t--Cost;\t\tQCost:" << -1*gaNCVCost << std::endl;
			}
			else
			{
				if(verbose)
					std::cout  << "[D4] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be merged.\t\t\t--Cost;\t\tQCost:" << -1*gbNCVCost << std::endl;
				v.push_back(-1*gbNCVCost);
			}
			x.push_back( v );
		}

		if( verify_rule_Rfour( ga, gb ) )
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[R4] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be interchanged.\t\tNo cost change;" << std::endl;
			v.push_back(4);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			v.push_back(0);
			x.push_back( v );
		}
		
		if( verify_rule_Dfive( ga, gb ) )
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[D5] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can insert controls.\t\tNo cost change;\tQCost:" << ncv_cost(gaControls+1) + ncv_cost(gbControls+1) << std::endl;
			v.push_back(5);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			v.push_back(ncv_cost(gaControls+1) + ncv_cost(gbControls+1));
			x.push_back( v );
		}
		
		if( verify_rule_Dfivee( ga, gb ) )
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[D5] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can remove controls.\t\tNo cost change;\tQCost:" << -1*(ncv_cost(gaControls+1) + ncv_cost(gbControls+1)) <<std::endl;
			v.push_back(-5);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			v.push_back(-1*(ncv_cost(gaControls+1) + ncv_cost(gbControls+1)));
			x.push_back( v );
		}

		if( verify_rule_Dsix( ga, gb ) )
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[D6] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be interchanged.\t\tNo cost change;" << std::endl;
			v.push_back(6);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(0);
			v.push_back(0);
			x.push_back( v );
		}

		if( verify_rule_Dseven( ga, gb ) )
		{
			std::vector<int> v;
			if(verbose)
				std::cout  << "[D7] Gates ( " << itGateIndex << " - " << itGateIndex + 1 << " ) can be interchanged.\t\t++Cost;\t\tQCost:" << ncv_cost(gaControls + gbControls) << std::endl;
			v.push_back(7);
			v.push_back(itGateIndex);
			v.push_back(nextGateIndex);
			v.push_back(1);
			v.push_back(ncv_cost(gaControls + gbControls));
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

void apply_rule( circuit& circ, std::vector<std::vector<int>> x, unsigned y, bool verbose )
{
	circuit::const_iterator itGate = circ.begin(), nextGate = circ.begin();
	switch( x[y][0] )
	{
		case 1:
			apply_rule_Done( circ, x[y][1], x[y][2] );
			std::cout << "Rule D1 applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case 2:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dtwo( itGate, nextGate );
			std::cout << "Rule D2 applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case 3:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dthree( circ, itGate, nextGate );
			std::cout << "Rule D3 applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case 4:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			swap_gates( itGate, nextGate );
			std::cout << "Rule swap applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case 5:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dfive( itGate, nextGate );
			std::cout << "Rule D5 applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case -5:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dfivee( itGate, nextGate );
			std::cout << "Rule D5.1 applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case 6:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			swap_gates( itGate, nextGate );
			std::cout << "Rule swap applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case 7:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dseven( circ, itGate, nextGate );
			std::cout << "Rule D7 applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		case 14:
			std::advance( itGate, x[y][1] );
			std::advance( nextGate, x[y][2] );
			apply_rule_Dfour( circ, itGate, nextGate );
			std::cout << "Rule D4 applied " << x[y][1] << " " << x[y][2] << std::endl;
			break;
		default:
			std::cout << "Something got wrong!" << std::endl;
			break;
	}
}

//Sort the vector of the rules
void sortrows(matrix& x, int col)
{    
	std::sort(	x.begin(),
	      		x.end(),
	      		[col](const std::vector<int>& lhs, const std::vector<int>& rhs) {
	          		return lhs[col] < rhs[col];
	      });
}

//Print vector
void print_list(matrix& x)
{
	for(unsigned i = 0; i < x.size(); ++i)
	{
		for(unsigned j = 0; j < x[i].size(); ++j)
			std::cout << x[i][j] << " ";
		std::cout << std::endl;
	}

}

//Function to decide which rule must be applied
void choosing_rule(circuit& circ, matrix x, matrix& tp, bool verbose)
{
	bool insert;
	//std::cout << t.size() << std::endl;
	for(unsigned k = 0; k < x.size(); ++k)
	{
		if(x[k][3] < 0)
		{
			apply_rule(circ, x, k, verbose);
			break;
		}
		insert = true;
		for(unsigned i = 0; i < tp.size(); ++i)
		{
			if(x[k][0] == tp[i][0] && x[k][1] == tp[i][1] && x[k][2] == tp[i][2])
			{
				insert = false;
				break;
			}
		}
		if(insert)
		{
			x[k][3] = 0;
			tp.push_back( x[k] );
			apply_rule( circ, x, k, verbose );
			break;
		}
	}
	for(unsigned i = 0; i < tp.size(); ++i)
		++tp[i][3];
}

//Function to decide which rule must be applied
void choosing_rule(circuit& circ, matrix x, matrix& tp)
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
	for(unsigned i = 0; i < tp.size(); ++i)
		++tp[i][3];	
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
bool update_circuit(circuit& circ, circuit& min, unsigned opt)
{
	switch(opt){
		case 0:
			if(circ.num_gates() < min.num_gates())
			{
				clear_circuit(min);
				copy_circuit(circ, min);
				//std::cout << " Minimo: " << min.num_gates() << std::endl;
				return true;
			}
		default:
			if(circuit_ncv_cost(circ) < circuit_ncv_cost(min))
			{
				clear_circuit(min);
				copy_circuit(circ, min);
				//std::cout << " Minimo: " << min.num_gates() << std::endl;
				return true;
			}
	}
	return false;
}

void tabu_search( circuit& circ, unsigned overlap, unsigned neighborhood, const properties::ptr& statistics, unsigned opt )
{
	properties_timer t( statistics );
	unsigned stop = 0, begin, end;
	bool finish = false;
	matrix x; //list of possible rules of the current iteration
	matrix tp; //tabu list with the penalization
	circuit orig, min;

	copy_circuit( circ, orig );
 	copy_circuit( circ, min );
 	// reverse_circuit( orig );
 	begin = 0;
	end = neighborhood;
	overlap = (overlap * neighborhood) * 0.01;

	while(!finish)
	{ 	
		if(end > circ.num_gates())
		{
			end = circ.num_gates();
			if(end < neighborhood)
			{
				begin = 0;
			}
			else
			{
				begin = end - neighborhood;
			}
			finish = true;
		}
		else if(int(begin) < 0)
		{
			begin = 0;
		}
		//std::cout << "begin: " << begin << " end: " << end << std::endl;
	 	stop = 0;
	 	unsigned int improvement = circ.num_gates();
	 	tp.clear();
	 	while(stop < neighborhood * 50)
	 	{
	 		list_rules(circ, x, begin, end, false);
	 		sortrows(x, opt+3);
	 		choosing_rule(circ, x, tp);
	 		update_tabu_list(tp, neighborhood);
	 		if(!update_circuit(circ, min, opt))
	 			++stop;
	 		else
	 			stop = 0;
	 		x.clear();

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
	 	begin = begin + (neighborhood - overlap) - improvement;
	 	end = begin + neighborhood;
	}
	clear_circuit(circ);
	copy_circuit(min, circ);
}

void tabu_search( circuit& circ, unsigned overlap, unsigned neighborhood, const properties::ptr& statistics, unsigned opt, bool verbose, bool step )
{
	properties_timer t( statistics );
	unsigned stop = 0, begin, end;
	bool finish = false;
	matrix x; //list of possible rules of the current iteration
	matrix tp; //tabu list with the penalization
	circuit orig, min;

	copy_circuit( circ, orig );
 	copy_circuit( circ, min );
 	// reverse_circuit( orig );
 	begin = 0;
	end = neighborhood;
	overlap = (overlap * neighborhood) * 0.01;

	while(!finish)
	{ 	
		if(end > circ.num_gates())
		{
			end = circ.num_gates();
			if(end < neighborhood)
			{
				begin = 0;
			}
			else
			{
				begin = end - neighborhood;
			}
			finish = true;
		}
		else if(int(begin) < 0)
		{
			begin = 0;
		}
		//std::cout << "begin: " << begin << " end: " << end << std::endl;
	 	stop = 0;
	 	tp.clear();
	 	unsigned int improvement = circ.num_gates();
	 	while(stop < neighborhood * 50)
	 	{
	 		std::cout << "++++++++++ ITERATION " << stop << " +++++++++++" << std::endl;
	 		std::cout << circ << std::endl;
	 		list_rules(circ, x, begin, end, verbose);
	 		sortrows(x, opt+3);
	 		std::cout << "++++++++++ BEGIN LIST OF RULES +++++++++++" << std::endl;
 			print_list( x );
	 		std::cout << "++++++++++   END LIST OF RULES +++++++++++" << std::endl;
	 		choosing_rule(circ, x, tp, verbose);
	 		update_tabu_list(tp, neighborhood);
 			std::cout << "++++++++++ BEGIN TABU LIST +++++++++++" << std::endl;
 			print_list( tp );
 			std::cout << "++++++++++   END TABU LIST +++++++++++" << std::endl;	
	 		if(!update_circuit(circ, min, opt))
	 		{
	 			++stop;
	 		}
	 		else
	 		{
	 			stop = 0;
	 			std::cout << "++++++++++ CIRCUIT UPDATED +++++++++++" << std::endl;
	 		}
	 		x.clear();
	 		if(step)
 				std::cin.get();

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
	 	begin = begin + (neighborhood - overlap) - improvement;
	 	end = begin + neighborhood;
	}
	clear_circuit(circ);
	copy_circuit(min, circ);
}


bool tabu_command::execute()
{
	circuit aux, circ;
	auto& circuits = env->store<circuit>();
	aux = circuits.current();
	copy_circuit(aux, circ);
 	unsigned int InitialGatesCost, InitialQuantumCost = 0, FinalQuantumCost = 0;
 	InitialGatesCost = circ.num_gates();
 	
 	InitialQuantumCost = circuit_ncv_cost(circ);

 	if ( is_set("verbose") && is_set("step") )
 		tabu_search( circ, overlap, neighborhood, statistics, opt, true, true );
 	else if ( is_set("verbose") && !is_set("step") )
 		tabu_search( circ, overlap, neighborhood, statistics, opt, true, false );
 	else
 		tabu_search( circ, overlap, neighborhood, statistics, opt );
	if ( is_set( "new" ) )
        circuits.extend();    
 	circuits.current() = circ;

 	FinalQuantumCost = circuit_ncv_cost(circ);

    std::cout << " Begin gates: " << InitialGatesCost <<std::endl; 
    std::cout << " Final gates: " << circ.num_gates() << std::endl;
    std::cout << " Begin quantum cost: " << InitialQuantumCost <<std::endl; 
    std::cout << " Final quantum cost: " << FinalQuantumCost << std::endl;
    print_runtime();
    
    if (!changed_circuit(circ))
    	std::cout << "BRRRRRR" << std::endl;

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
