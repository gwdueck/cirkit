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

void print_line(xt::xarray<complex_t>& t, unsigned int& q)
{
	unsigned int line;
	std::string in;
  	std::cout << "Input: ";
  	std::cin >> in;

  	line = std::stoi(in, nullptr, 2);
  	for (int i = line*q; i < (line+1)*q; ++i)
  		std::cout << " " << t[i];
  	
 	std::cout << std::endl;
}

void print_solution(xt::xarray<complex_t>& t, unsigned int& q)
{
	unsigned int j = 0;
  	for (int i = 0; i < t.size(); ++i, ++j)
  	{
		if(i % q == 0)
  		{
 			std::cout << std::endl;
 			std::cout << "input: " << std::bitset<4>(i/q).to_string() << " => ";
  			j = 0;
  		}
  	  	if(t[i] != 0.0)
  		{
  			std::cout << " " << t[i] << " " << std::bitset<4>(j).to_string();
  		}
  	}
 	std::cout << std::endl;
}

bool alex_command::execute()
{
	const auto& circuits = env->store<circuit>().current();

	xt::xarray<complex_t> teste;

  	teste = matrix_from_clifford_t_circuit( circuits, is_set( "progress" ) );

  	// Number of qubits
  	unsigned int qubits = sqrt(teste.size());
 	// std::cout << "qubits: " << qubits << std::endl;

	// print_all_matrix(teste, qubits);
	// print_line(teste, qubits);
	print_solution(teste, qubits);


	return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
