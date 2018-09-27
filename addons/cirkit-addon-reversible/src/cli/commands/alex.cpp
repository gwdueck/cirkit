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
#include <iostream>

#include <boost/program_options.hpp>

#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include <cli/reversible_stores.hpp>
#include <reversible/utils/matrix_utils.hpp>
#include <alice/rules.hpp>


namespace cirkit
{

using boost::program_options::value;

alex_command::alex_command( const environment::ptr& env )
    : cirkit_command( env, "Alex test" )
{
	opts.add_options()
	( "input,i",	value( &input ),	"print single input -- use it in binary" )
    ;

}

command::rules_t alex_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}


void print_all_matrix(xt::xarray<complex_t>& t, unsigned int& q)
{
  	for (int i = 0; i < t.size(); ++i)
  	{
  		if(i % q == 0)
 			std::cout << std::endl;
  		std::cout << " " << t[i];
  	}
 	std::cout << std::endl;
}

void print_line(xt::xarray<complex_t>& t, unsigned int& q, std::string input)
{
	unsigned int j = 0;
	bool first = true;
	std::string p;
	unsigned int k = log10(q)/log10(2);

	unsigned int line;
  	line = std::stoi(input, nullptr, 2);
  	if(line >= q)
  	{
  		std::cout << "Out of range. The circuit has " << k << " qubits." << std::endl;
  	}
  	else
  	{
  		for (int i = line*q; i < (line+1)*q; ++i, ++j)
	  	{
			if(t[i] != 0.0)
	  		{
	  			if(first)
	  			{
	  				p = std::bitset<8>(j).to_string();
	  				// std::cout << std::endl;
					std::string pp = std::bitset<8>(line).to_string();
					std::cout << "input: " << pp.substr(pp.length() - k) << " => " << t[i] << " " << p.substr(p.length() - k);
	  				first = false;
	  			}
	  			else
	  			{
	  				p = std::bitset<8>(j).to_string();
	  				std::cout << " + " << t[i] << " " << p.substr(p.length() - k);
	  			}
	  		}
	  	}
	  	std::cout << std::endl;
  	}
}

void print_solution(xt::xarray<complex_t>& t, const unsigned int& q)
{
	unsigned int j = 0;
	bool first = true;
	std::string p;
	unsigned int k = log10(q)/log10(2);
  	for (int i = 0; i < t.size(); ++i, ++j)
  	{
		if(i % q == 0)
  		{
  			p = std::bitset<8>(i/q).to_string();
 			std::cout << std::endl;
 			std::cout << "input: " <<  p.substr(p.length() - k) << " => ";
  			j = 0;
  			first = true;
  		}
  	  	if(t[i] != 0.0)
  		{
  			if(first)
  			{
  				p = std::bitset<8>(j).to_string();
  				std::cout << " " << t[i] << " " << p.substr(p.length() - k);
  				first = false;
  			}
  			else
  			{
  				p = std::bitset<8>(j).to_string();
  				std::cout << " + " << t[i] << " " << p.substr(p.length() - k);
  			}
  		}
  	}
 	std::cout << std::endl;

  // std::cout << t << std::endl;
}

bool alex_command::execute()
{
	const circuit circ = env->store<circuit>().current();

  unsigned target, control;
  std::vector<unsigned> v;
  std::vector<std::vector<unsigned>> output;
  int static const map_method_qx2[5][5] ={{0,0,0,10,10},
                                          {4,0,0,10,10},
                                          {4,4,0,4,4},
                                          {10,10,0,0,0},
                                          {10,10,0,4,0}};

  // Create a matrix with 0's
  for (int i = 0; i < circ.lines(); ++i){
    for (int j = 0; j < circ.lines(); ++j){
      v.push_back(0);
    }
    output.push_back(v);
  }
  
  // Create a matrix with the cnots 
  for ( const auto& gate : circ )
  {
    if( !gate.controls().empty() ) // if is not a NOT gate
    {
      target = gate.targets().front();
      control = gate.controls().front().line();
      ++output[control][target];
    }
  }

  // // Print the cost of the transformations
  // for (int i = 0; i < 5; ++i){
  //   for (int j = 0; j < 5; ++j){
  //     std::cout << "\t" << map_method_qx2[i][j];
  //   }
  //   std::cout << std::endl;
  // }
  
  // std::cout << std::endl;

  // Print the cnots in the circuit
  for (int i = 0; i < circ.lines(); ++i){
    for (int j = 0; j < circ.lines(); ++j){
      std::cout << "\t" << output[i][j];
    }
    std::cout << std::endl;
  }

  std::cout << std::endl;
  
  // // Print the total cost of the circuit
  // for (int i = 0; i < circ.lines(); ++i){
  //   for (int j = 0; j < circ.lines(); ++j){
  //     std::cout << "\t" << output[i][j]*map_method_qx2[i][j];
  //   }
  //   std::cout << std::endl;
  // }

  // std::cout << std::endl;
  
  // // Print a possible approach to hungarian algorithm
  // for (int i = 0; i < circ.lines(); ++i){
  //   for (int j = 0; j < circ.lines(); ++j){
  //     std::cout << "\t" << (output[i][j]*map_method_qx2[i][j])+map_method_qx2[i][j]+1;
  //   }
  //   std::cout << std::endl;
  // }  





	// xt::xarray<complex_t> table;

 //  	table = matrix_from_clifford_t_circuit( circ, is_set( "progress" ) );

 //  	// Number of qubits
 //  	unsigned int qubits = sqrt(table.size());
 // 	// std::cout << "qubits: " << qubits << std::endl;

	// // print_all_matrix(table, qubits);
	// if(is_set("input"))
	// 	print_line(table, qubits, input);
	// else
	// 	print_solution(table, qubits);


	return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
