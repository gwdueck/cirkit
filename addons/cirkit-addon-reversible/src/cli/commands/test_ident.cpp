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

#include "test_ident.hpp"

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
 
test_ident_command::test_ident_command(const environment::ptr& env)
    : cirkit_command(env, "Tabu Search")
{
	opts.add_options()
    ( "verbose,v",	"be verbose")
    ;
    add_new_option();
}

command::rules_t test_ident_command::validity_rules() const
{
	return {has_store_element<circuit>( env )};
}



bool test_ident_command::execute()
{
	std::cout << "here we are\n";

	return true;
}

command::log_opt_t test_ident_command::log() const
{
  return log_opt_t({{"runtime", statistics->get<double>( "runtime" )}});
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
