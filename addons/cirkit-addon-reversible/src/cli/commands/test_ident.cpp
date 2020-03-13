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
#include <fstream>
#include <string>
#include <time.h>
#include <core/utils/timer.hpp>
#include <alice/rules.hpp>
#include <cli/reversible_stores.hpp>
#include <reversible/circuit.hpp>
#include <reversible/io/print_circuit.hpp>
#include <reversible/io/read_qc.hpp>
#include <reversible/io/write_qc.hpp>
#include <reversible/functions/add_circuit.hpp>
#include <reversible/functions/copy_circuit.hpp>
#include <reversible/functions/clear_circuit.hpp>
#include <reversible/functions/is_identity.hpp>
#include <reversible/functions/reverse_circuit.hpp>
#include <reversible/functions/remove_dup_gates.hpp>
#include <reversible/functions/match_templates.hpp>
#include <reversible/functions/clifford_templates.hpp>


#include <cli/commands/rules.hpp>
#include <core/utils/program_options.hpp>


namespace cirkit
{
 
test_ident_command::test_ident_command(const environment::ptr& env)
    : cirkit_command(env, "testing identities")
{
	opts.add_options()
    ( "verbose,v",	"be verbose")
    ( "filter,f", "filter identities")
    ( "read_templ,r", "read templates")
    ( "print_templ,p", "print templates")
    ( "temp,t", "test and develop new ideas")
    ( "circ_to_templ,c", "convert a circuit to a template")
    ( "experimental,e", "experiment with new functionality")
    ( "add_templ,a", "add new templates")
    ( "print_templates,g", "print all templates graphically as identity circuits")
    ;
 //   add_new_option();
}

/*
command::rules_t test_ident_command::validity_rules() const
{
	 return {has_store_element<circuit>( env )};
}
*/


bool test_ident_command::execute()
{
	if( is_set( "read_templ" ) )
	{
		std::string filename;
    	std::cout << "Enter the file name for the templates ";
      	std::cin >> filename;
      	std::ifstream templ_file ( filename );
      	int n_templs; // number of templates
      	templ_file >> n_templs;
      	Clifford_Template new_templ;
      	for( int i = 0; i < n_templs; i++ )
      	{
      		new_templ.read( templ_file );
      		cliff_templates.push_back( new_templ );
      		new_templ.clear();
      	}
	}

	if( is_set( "print_templ" ) )
		{
			for (std::vector<Clifford_Template>::iterator it = cliff_templates.begin() ; it != cliff_templates.end(); ++it)
			{
				it->print();
			}
		}


	if( is_set( "filter" ) )
	{
		std::ifstream fileList ("file_list.txt");
		std::ofstream filterfile ("filter/file_list.txt");
		std::string infile_qc;
		circuit circ_working, circ_reduced;
		int not_reduced = 0, n = 0;

	 	if ( !fileList.is_open() )
	 	{
	 		std::cout << "ERROR cannot open file file_list.txt\n";
	 		return false;
	 	}

	 	fileList >> infile_qc;
	 	while( !fileList.eof() )
	 	{
	 		n++;
	 		circ_working = read_qc( infile_qc );
	 		circ_reduced = remove_dup_gates( circ_working );
	 		if( circ_working.num_gates() == circ_reduced.num_gates() )
	 		{
	 			bool flag = match_any_template( circ_reduced, cliff_templates );
	 		}
	 		if( circ_working.num_gates() == circ_reduced.num_gates() )
	 		{
	 			filterfile << infile_qc << std::endl;
	 			write_qc( circ_working, "filter/" + infile_qc, false );
	 			not_reduced++;
	 		}
	 		fileList >> infile_qc;
	 		if ( n%10000 == 0 )
	 		{
	 			std::cout << not_reduced << "\n";
	 		}
	 	}
	 	filterfile.close();
	}

	// read templates and test if they are not reducible by previously defined templates
	if( is_set( "temp" ) )
	{
		// auto& circuits = env->store<circuit>();
		std::string filename;
    	std::cout << "Enter the file name for the templates ";
      	std::cin >> filename;
      	std::ifstream templ_file ( filename );
      	int n_templs; // number of templates
      	templ_file >> n_templs;
      	Clifford_Template new_templ;
      	for( int i = 0; i < n_templs; i++ )
      	{
      		new_templ.read( templ_file );
      		circuit circ_reduced;
      		circuit circ = new_templ.convert_to_circ( true );
      		//circuits.extend();
      		//circuits.current() = circ;
      		std::cout << "\n";
      		new_templ.print();
      		std::cout << circ;
      		circ = new_templ.convert_to_circ( false );
      		circ_reduced = remove_dup_gates( circ );
      		bool flag1 = circ.num_gates() > circ_reduced.num_gates();
      		bool flag2 =  match_any_template( circ_reduced, cliff_templates );
    		std::cout << "flag1 = " << flag1 <<  " flag2 = " << flag2 <<  std::endl;
    		if( !flag1 && !flag2 )
    		{
      			cliff_templates.push_back( new_templ );
      		}
      		else
      		{
      			std::cout << "reduced circuit\n";
      			std::cout << circ_reduced;
      		}
      		new_templ.clear();
      	}
	}
	if( is_set( "circ_to_templ" ) )
	{
		auto& circuits = env->store<circuit>();
    	circuit circ_working = circuits.current();
    	Clifford_Template my_temp;
    	my_temp.convert_circ( circ_working );
    	my_temp.print();
	}
	if( is_set( "experimental" ) )
	{
		circuit circ_working, circ_reduced;
		std::ifstream fileList ("file_list.txt");
		std::ofstream tex_table ("table.tex");
		std::string infile_qc;
		if ( !fileList.is_open() )
	 	{
	 		std::cout << "ERROR cannot open file file_list.txt\n";
	 		return false;
	 	}

	 	fileList >> infile_qc;
	 	while( !fileList.eof() )
	 	{
	 		std::cout << "read " << infile_qc << std::endl;
	 		circ_working = read_qc( infile_qc );
	 		if( circ_working.num_gates() < 100000)
	 		{
		 		circ_reduced = remove_dup_gates( circ_working );
		 		bool flag = match_any_template( circ_reduced, cliff_templates );
		 		if( flag )
		 		{
		 			circ_reduced = remove_dup_gates( circ_reduced );
		 		}
		 		if( circ_working.num_gates() > circ_reduced.num_gates() )
		 		{
		 			circ_reduced = remove_dup_gates( circ_reduced );
		 			std::cout << "Success " << infile_qc << " reduced from " << circ_working.num_gates() <<
		 			" gates to " << circ_reduced.num_gates() << std::endl;
		 		
		 			tex_table << infile_qc << " & " << circ_working.num_gates() << " & " << 
		 				circ_reduced.num_gates() << " & "  << circ_working.num_gates() - circ_reduced.num_gates() << " & " <<
		 				(int) ( circ_working.num_gates() - circ_reduced.num_gates() ) * 100 / circ_working.num_gates() 
		 				<<  "\\% \\\\ \\hline\n";
		 			tex_table.flush();
		 		}
		 		fileList >> infile_qc;
		 	}
	 	}
	 	tex_table.close();
	}
	if( is_set( "add_templ" ))
	{
		circuit circ_working, circ_reduced;
		Clifford_Template my_temp;
		std::ifstream fileList ("file_list.txt");
		std::string infile_qc;
		if ( !fileList.is_open() )
	 	{
	 		std::cout << "ERROR cannot open file file_list.txt\n";
	 		return false;
	 	}

	 	fileList >> infile_qc;
	 	while( !fileList.eof() )
	 	{
	 		std::cout << "read " << infile_qc << std::endl;
	 		circ_working = read_qc( infile_qc );
	 		circ_reduced = remove_dup_gates( circ_working );
	 		if( circ_working.num_gates() == circ_reduced.num_gates() )
	 		{
	 			bool flag = match_any_template( circ_reduced, cliff_templates );
	 		}
	 		if( circ_working.num_gates() == circ_reduced.num_gates() )
	 		{
	 			my_temp.convert_circ( circ_working );
	 			cliff_templates.push_back( my_temp );
	 			my_temp.print();
	 			std::cout << "add new template " << infile_qc << std::endl;
	 		}
	 		else
	 		{
	 			std::cout << "reduced circuit\n";
	 			std::cout << circ_reduced;
	 		}
	 		my_temp.clear();
	 		fileList >> infile_qc;
	 	}
	}

	if( is_set( "print_templates" ))
	{
		circuit circ;
		for( auto & templ : cliff_templates )
		{
			circ = templ.convert_to_circ( true );
			templ.print();
			std::cout << circ;
		}
	}

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
