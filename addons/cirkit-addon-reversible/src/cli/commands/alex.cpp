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




std::string vector_to_string(std::vector<unsigned>& ppp)
{
	std::string p;
	for (int i = 0; i < ppp.size(); ++i)
		p.append(std::to_string(ppp[i]));
	return p;
}

bool alex_command::execute()
{
	const auto& circuits = env->store<circuit>().current();

  	// std::vector<std::pair<xt::xarray<complex_t>, unsigned>> matrices;

	xt::xarray<complex_t> teste;

  	teste = matrix_from_clifford_t_circuit( circuits, is_set( "progress" ) );

  	std::cout << "teste: " <<  teste.size() <<std::endl;
  	std::cout << teste << std::endl;
  	std::cout << "teste1: " << std::endl;
 	
 	unsigned qubits = sqrt(teste.size());
 	std::cout << "qubits: " << qubits << std::endl;

  	for (int i = 0; i < teste.size(); ++i)
  	{
  		if(i % qubits == 0)
 			std::cout << std::endl;
  		std::cout << " " << teste[i];
  	}
 	std::cout << std::endl;
	return true;
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
