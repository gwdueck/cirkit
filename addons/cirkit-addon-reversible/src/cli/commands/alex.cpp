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
    ( "penalty,p",   value_with_default( &numberPenalty ), "penalty number (times the number of combination)" )
    ( "iterations,i",   value_with_default( &numberIterations ), "number of iterations" )

	;

}

command::rules_t alex_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

void alcMatrix( matrix& m, unsigned size )
{
	// std::cout << "Allocating matrix..." << std::endl;
  	std::vector<unsigned> v;
	for (int i = 0; i < size; ++i)
		v.push_back(0);
	for (int i = 0; i < size; ++i)
		m.push_back(v);
}

void genMatrix( circuit& circ, matrix& m )
{
	// std::cout << "Generating matrix..." << std::endl;	
  	unsigned target, control;
	for ( const auto& gate : circ )
	{
		if( gate.controls().size() == 1 )
		{
		  target = gate.targets().front();
		  control = gate.controls().front().line();
		  if( is_toffoli( gate ) )
		  	++m[control][target];
		}
	}
}

void prtMatrix( matrix& m )
{
	// std::cout << "Printing matrix..." << std::endl;
	for (int i = 0; i < m.size(); ++i)
	{
		for (int j = 0; j < m[i].size(); ++j)
			std::cout << " " << m[i][j];
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

unsigned int costMatrix( matrix& m, matrix& c )
{
	// std::cout << "Printing matrix..." << std::endl;
	unsigned cost = 0;
	for (int i = 0; i < m.size(); ++i)
		for (int j = 0; j < m[i].size(); ++j)
			cost += m[i][j]*c[i][j];
	return cost;
}

matrix qx5 ={{0,4,10,13,19,29,39,51,61,55,45,41,29,20,10,4},
			{0,0,0,3,9,19,29,41,51,61,53,41,29,19,9,10},
			{10,4,0,0,3,13,23,35,45,55,45,35,23,13,3,4},
			{19,7,4,0,0,10,19,31,41,45,35,31,19,10,0,7},
			{25,13,7,4,0,4,7,19,29,33,23,19,7,4,10,13},
			{35,25,19,10,0,0,4,10,20,23,13,10,4,10,13,23},
			{45,33,23,13,3,0,0,0,10,13,3,0,10,13,25,33},
			{55,43,33,23,13,10,4,0,4,7,0,10,19,23,35,43},
			{67,55,45,35,25,22,10,0,0,4,3,13,23,35,45,55},
			{57,65,55,45,35,25,13,3,0,0,0,10,19,31,35,45},
			{45,53,43,33,23,19,7,4,7,4,0,4,7,19,23,33},
			{35,43,33,23,13,10,4,10,19,10,0,0,4,10,13,23},
			{25,33,23,13,3,0,10,13,23,13,3,0,0,0,3,13},
			{22,25,19,10,0,10,19,23,33,23,13,10,4,0,0,10},
			{10,13,7,4,10,19,29,33,43,33,23,19,7,4,0,4},
			{0,10,0,3,9,19,29,41,51,45,35,31,19,10,0,0}};

matrix qx4 ={{0,4,4,7,7},{0,0,4,7,7},{0,0,0,4,4},{3,3,0,0,0},{3,3,0,4,0}};

void permuteLines( matrix& cnots, unsigned x, unsigned y )
{
	for (int i = 0; i < cnots.size(); ++i)
		std::swap( cnots[i][x], cnots[i][y] );
	for (int i = 0; i < cnots.size(); ++i)
		std::swap( cnots[x][i], cnots[y][i] );
}

void genDelta( matrix& delta, unsigned s )
{
	delta.push_back({0,0,0,0});
	for (int i = 0; i < s; ++i)
		for (int j = i + 1; j < s; ++j)
			delta.push_back({i,j,0,0});
}

void calcDelta( matrix& cnots, matrix qxCost, matrix& delta )
{
	for (int i = 0; i < delta.size(); ++i)
	{
		// std::cout << "Permuting: " << delta[i][0] << " " << delta[i][1];		
		// std::cout << " Cost ";
		permuteLines( cnots, delta[i][0], delta[i][1] );
		delta[i][2] = costMatrix( cnots, qxCost );
		// std::cout << costMatrix( cnots, qxCost ) << std::endl;
		permuteLines( cnots, delta[i][0], delta[i][1] );
	}
}

std::pair<unsigned,unsigned> chooseDelta( matrix& delta, unsigned& bestCost )
{
	std::pair<unsigned,unsigned> d;
	// prtMatrix(m);
	std::sort(delta.begin(), delta.end(),
          [](const std::vector<unsigned>& a, const std::vector<unsigned>& b) {
  		return a[2] < b[2];
	});
	// prtMatrix(m);
	for (int i = 0; i < delta.size(); ++i)
	{
		// if ( delta[i][3] == 0 || delta[i][2] < bestCost )
		if ( delta[i][3] == 0 || delta[i][2] < bestCost )
		{
			d = std::make_pair( delta[i][0], delta[i][1] );
			++delta[i][3];
			break;
		}
	}	
	return d;
}

void updateDelta( matrix& delta, unsigned numberPenalty )
{
	for (int i = 0; i < delta.size(); ++i)
		if ( delta[i][3] > numberPenalty*delta.size() )
			delta[i][3] = 0;
		else if ( delta[i][3] > 0 )
			++delta[i][3];
}


bool alex_command::execute()
{
	auto& circuits = env->store<circuit>();
	circuit circ = circuits.current();
	
	matrix cnots, qxCost, delta;
	std::vector<unsigned> permutation, bestPermutation;
	std::pair<unsigned,unsigned> d;
	unsigned int actCost, bestCost = -1;

	if ( circ.lines() == 5 )
		qxCost = qx4;
	else
		qxCost = qx5;

	// allocate matrix
	alcMatrix( cnots, circ.lines() );
	// generate the matrix from the circuit
	genMatrix( circ, cnots );
	// generate the delta matrix
	genDelta( delta, circ.lines() );
	// std::cout << delta.size() << std::endl;
	// print cnots matrix
	// prtMatrix( cnots );
	// print delta matrix
	// prtMatrix( delta );
	// create the default permutation
	for (int i = 0; i < circ.lines(); ++i)
		permutation.push_back(i);

	// initialize the best permutation vector
	bestPermutation = permutation;
	// Start algorithm
	for (int i = 0; i < numberIterations; ++i)
	{
		// std::cout << "Iteration: " << i << std::endl;
		// prtMatrix( cnots );
		// calculate the delta
		calcDelta( cnots, qxCost, delta );
		// choose the best delta
		d = chooseDelta( delta, bestCost );
		// prtMatrix( delta );
		// update the penalty
		updateDelta( delta, numberPenalty );
		// prtMatrix( delta );
		// std::cout << " " << d.first << " " << d.second;
		// permute for the best delta
		permuteLines( cnots, d.first, d.second );
		// update the permutation
		std::swap( permutation[d.first], permutation[d.second] );
		// std::cout << "\npermutation: ";
		// for (int i = 0; i < permutation.size(); ++i)
			// std::cout << " " << permutation[i];
		// std::cout << std::endl;
		// update the cost
		actCost = costMatrix( cnots, qxCost );
		// if better cost is found, update the bestPermutation
		if ( actCost < bestCost )
		{
			bestCost = actCost;
			bestPermutation = permutation;
		}
		// std::cout << " Cost: " << actCost << std::endl;

		// getchar();
	}
		
	// std::cout << "BP: ";
	for (int i = 0; i < bestPermutation.size(); ++i)
		std::cout << " " << bestPermutation[i];
	std::cout << "\t" << bestCost << std::endl;
	return true; 
}

// bool alex_command::execute()
// {
// 	v.clear();
// 	auto& circuits = env->store<circuit>();
// 	circuit circ = circuits.current();
	

// 	unsigned target, control, vgate=0, toffoligate=0, cnotgate=0, others=0;
// 	for ( const auto& gate : circ )
// 	{
// 		if( gate.controls().size() == 1 )
// 		{
// 		  if( is_toffoli( gate ) )
// 		  	++cnotgate;
// 		  else if( is_v( gate ) )
// 		  	++vgate;
// 		}
// 		else if( gate.controls().size() == 2 )
// 		{
// 		  if( is_toffoli( gate ) )
// 		  	++toffoligate;
// 		}
// 		else
// 		{
// 			++others;
// 		}
// 	}
// 	std::cout << "cnot: " << cnotgate << " toffoli: " << toffoligate << " v: " << vgate << " others: " << others << std::endl;

// 	matrix cnots;
// 	alcMatrix( cnots, circ.lines() );
// 	genMatrix( circ, cnots );
// 	prtMatrix( cnots );
// 	pCnots( cnots );
// 	return true; 
// 	// std::cout << "	" << circ.num_gates() << std::endl;
// 	// std::sort(v.begin(), v.end(), [](const vector<int> & a, const vector<int> & b){ return a.size() < b.size(); });
// 	std::sort(v.begin(), v.end(), std::greater<>());
// 	// std::sort(v.begin(), v.end());
// 	// std::cout << std::endl;
// 	bool seg = false;
// 	unsigned segValue = 0;
// 	std::vector< std::pair< int, std::pair< int,int> >> f;
// 	std::vector< std::pair< int, std::pair< int,int> >> s;

// 	for (int i = 0; i < v.size(); ++i)
// 	{
// 		// std::cout << "[" << v[i].second.first << "," << v[i].second.second << "] = " << v[i].first << std::endl; 

// 		if( v[i].first == v[0].first)
// 		{
// 			// std::cout << "[" << v[i].second.first << "," << v[i].second.second << "] -> " << v[i].first << " : "; 
// 			// std::cout << "[" << v[i].second.first << "," << v[i].second.second << "] = " << v[i].first << std::endl; 
// 			unsigned c = 0, t = 0;
// 			for (int j = 0; j < cnots.size(); ++j)
// 			{
// 				c = c + cnots[j][v[i].second.first];
// 				t = t + cnots[v[i].second.second][j];
// 			}
// 			// std::cout << "c: " << c << " t: " << t << " T: " << t + c << std::endl;
// 			f.push_back(std::make_pair(t+c,std::make_pair(v[i].second.first,v[i].second.second)));
// 		}
// 		else if(seg == false)
// 		{
// 			seg = true;
// 			segValue = cnots[v[i].second.first][v[i].second.second];
// 		}
// 		if( v[i].first == segValue )
// 		{
// 			// std::cout << "{" << v[i].second.first << "," << v[i].second.second << "} -> " << v[i].first << " : "; 
// 			// std::cout << "[" << v[i].second.first << "," << v[i].second.second << "] = " << v[i].first << std::endl; 
// 			unsigned c = 0, t = 0;
// 			for (int j = 0; j < cnots.size(); ++j)
// 			{
// 				c = c + cnots[j][v[i].second.first];
// 				t = t + cnots[v[i].second.second][j];
// 			}
// 			// std::cout << "c: " << c << " t: " << t << " T: " << t + c << std::endl;
// 			s.push_back(std::make_pair(t+c,std::make_pair(v[i].second.first,v[i].second.second)));
// 		}
// 	}
// 	// std::cout << std::endl;

// 	std::sort(f.begin(), f.end(), std::greater<>());
// 	for (int i = 0; i < f.size(); ++i)
//  	{
// 		// if( f[i].first == f[0].first)
// 		// {
// 			unsigned c = 0, t = 0;
// 			for (int j = 0; j < cnots.size(); ++j)
// 			{
// 				c = c + cnots[j][f[i].second.first];
// 				t = t + cnots[f[i].second.second][j];
// 			}
// 			// std::cout << "[" << f[i].second.first << "," << f[i].second.second << "] ";
// 			std::cout << "[" << f[i].second.first << "," << f[i].second.second << "] -> " << c << " " << t << std::endl;
// 			// std::cout << "[" << f[i].second.first << "," << f[i].second.second << "] -> " << f[i].first << std::endl;
// 		// }
// 	}
// 	std::sort(s.begin(), s.end(), std::greater<>());
// 	for (int i = 0; i < s.size(); ++i)
// 	{
// 		// if( s[i].first == s[0].first)
// 		// {
// 			unsigned c = 0, t = 0;
// 			for (int j = 0; j < cnots.size(); ++j)
// 			{
// 				c = c + cnots[j][s[i].second.first];
// 				t = t + cnots[s[i].second.second][j];
// 			}
// 			// std::cout << "{" << s[i].second.first << "," << s[i].second.second << "} ";
// 			std::cout << "{" << s[i].second.first << "," << s[i].second.second << "} -> " << c << " " << t << std::endl;
// 			// std::cout << "{" << s[i].second.first << "," << s[i].second.second << "} -> " << s[i].first << std::endl; 
// 		// }
// 	}
// 	std::cout << std::endl;
// 	return true;
// }

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
