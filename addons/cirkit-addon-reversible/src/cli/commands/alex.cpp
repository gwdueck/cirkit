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
	;

}

command::rules_t alex_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}


bool alex_command::execute()
{
	auto& circuits = env->store<circuit>();
	circuit circ = circuits.current();
	std::cout << "	" << circ.num_gates() << std::endl;
	
	return true;
}

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
