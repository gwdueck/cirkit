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
	( "filename,f",    value( &filename ),    "name of the output file" )
	;

}

command::rules_t alex_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}


using matrix = std::vector<std::vector<unsigned>>;

// Return the higher value element in the matrix
int get_max_element(matrix& m, unsigned& l, unsigned& c)
{
    unsigned h = 0;
    for (int i = 0; i < m.size(); ++i)
    {
        for (int j = 0; j < m[i].size(); ++j)
        {
            if(i != j && m[i][j] > h)
            {
                h = m[i][j];
                l = i;
                c = j;
            }
        }
    }
    return h;
}

// return the number of different gates of the circuit
int getNumberDifGates( matrix& c )
{
	unsigned qtd = 0;
	for (int i = 0; i < c.size(); ++i)
		for (int j = 0; j < c[i].size(); ++j)
			if( c[i][j] > 0 )
				++qtd;
	return qtd;
}

// Function to print the objective function
void printObjectiveFunction( matrix& qx, matrix& cnot, unsigned difGates )
{
	unsigned aux = 0;
	std::cout << "/* Begin Objective Function */" << std::endl;
	std::cout << "min:\t";
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1)
							std::cout << qx[k][m] << "G" << i << j << "c" << k << m << ";";
						else
							std::cout << qx[k][m] << "G" << i << j << "c" << k << m << " + ";
						line = true;
					}
				}
			}
			if(line && aux < difGates-1)
			{
				std::cout << std::endl << "\t";
				++aux;
			}
			else if(line)
			{
				std::cout << std::endl;
			}
		}
	}
	std::cout << "/* End Objective Function */" << std::endl;
}

// Function to print the first restriction
void printFirstRestriction( matrix& qx, matrix& cnot )
{
	unsigned aux = 0;
	std::cout << "/* Begin First Restriction */" << std::endl;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 )
							std::cout << cnot[i][j] << "G" << i << j << "c" << k << m << " = " << cnot[i][j] << ";";
						else
							std::cout << cnot[i][j] << "G" << i << j << "c" << k << m << " + ";
						line = true;
					}
				}
			}
			if(line)
				std::cout << std::endl;
		}
	}
	std::cout << "/* End First Restriction */" << std::endl;
}

// Function to print the final restriction
void printEndRestriction( matrix& qx, matrix& cnot, unsigned difGates )
{
	unsigned aux = 0;
	std::cout << "/* Begin Final Restriction */" << std::endl;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1)
							std::cout << "G" << i << j << "c" << k << m << " = " << difGates << ";";
						else
							std::cout << "G" << i << j << "c" << k << m << " + ";
						line = true;
					}
				}
			}
			if(line)
			{
				std::cout << std::endl;
				++aux;
			}
		}
	}
	std::cout << "/* End Final Restriction */" << std::endl;
}

// Function to print the variables
void printIntegerVariables( matrix& qx, matrix& cnot, unsigned difGates )
{
	unsigned aux = 0;
	std::cout << "/* Begin Integer Variables */" << std::endl;
	std::cout << "int\t";
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnot[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1)
							std::cout << "G" << i << j << "c" << k << m << ";";
						else
							std::cout << "G" << i << j << "c" << k << m << "  ";
						line = true;
					}
				}
			}
			if(line && aux < difGates-1)
			{
				std::cout << std::endl << "\t";
				++aux;
			}
			else if(line)
			{
				std::cout << std::endl;
			}
		}
	}
	std::cout << "/* End Integer Variables */" << std::endl;
}

// Create a matrix with 0's
void createMatrix( matrix& m, unsigned size )
{
	std::cout << "Creating matrix..." << std::endl;
  	std::vector<unsigned> v;
	for (int i = 0; i < size; ++i)
		v.push_back(0);
	for (int i = 0; i < size; ++i)
		m.push_back(v);
}

// Create a matrix with the cnots 
void generateMatrixCnots( circuit& circ, matrix& m )
{
	std::cout << "Generating matrix..." << std::endl;	
  	unsigned target, control;
	for ( const auto& gate : circ )
	{
		if( !gate.controls().empty() )
		{
		  target = gate.targets().front();
		  control = gate.controls().front().line();
		  ++m[control][target];
		}
	}
}

// Print the cnots in the circuit
void printMatrixCnots( matrix& m )
{
	std::cout << "Printing matrix..." << std::endl;
	for (int i = 0; i < m.size(); ++i)
	{
		for (int j = 0; j < m[i].size(); ++j)
			std::cout << "\t" << m[i][j];
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// Only function to write almost everything
// type = 0 -> control - control
// type = 1 -> target - target
// type = 2 -> control - target
// type = 3 -> target - control
// type = 4 -> inverse (control -> target; target -> control)
void writeDep( unsigned l0, unsigned c0, unsigned l1, unsigned c1, unsigned size, unsigned type )
{
	switch ( type )
    {
    	case 0:
			std::cout << "/* Writing Control - Control dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 1:
			std::cout << "/* Writing Target - Target dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 2:
			std::cout << "/* Writing Control - Target dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 3:
			std::cout << "/* Writing Target - Control dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	case 4:
			std::cout << "/* Writing Inverse dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
    		break;
    	default:
			std::cout << "/* ERRORRRRRRRRRRRRRRRRR */" << std::endl;
    }
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			if(i != j && type != 4)
			{
				if(type == 0 || type == 2)
				{
					std::cout << "G" << l0 << c0;
					std::cout << "c" << i << j;
					std::cout << " <= ";
				}
				else
				{
					std::cout << "G" << l0 << c0;
					std::cout << "c" << j << i;
					std::cout << " <= ";
				}
				unsigned aux = 0;
				for (int k = 0; k < size; ++k)
				{
					if(k != j && k != i)
					{
						++aux;
						if(type == 0 || type == 3)
						{
							std::cout << "G" << l1 << c1;
							std::cout << "c" << i << k;	
						}
						else
						{
							std::cout << "G" << l1 << c1;
							std::cout << "c" << k << i;
						}
						if(aux == size-2)
							std::cout << ";";
						else	
							std::cout << " + ";	
					}
				}
				std::cout << std::endl;	
			}
			else if(i != j && type == 4)
			{
				std::cout << "G" << l0 << c0;
				std::cout << "c" << i << j;
				std::cout << " <= ";
				std::cout << "G" << l1 << c1;
				std::cout << "c" << j << i;
				std::cout << ";" << std::endl;
			}
		}
	}
}

void getAllCombinations(matrix& output)
{
	enum type { cc, tt, ct, tc, in };
	for (int i = 0; i < output.size(); ++i)
	{
		for (int j = 0; j < output[i].size(); ++j)
		{
			for (int m = 0; m < output.size(); ++m)
			{
				for (int n = 0; n < output[m].size(); ++n)
				{
					if(i != j && m != n && output[i][j] > 0 && output[m][n] > 0)
					{
						if(i == m && j != n)
							writeDep( i, j, m, n, output.size(),  cc);
						else if(i != m && j == n)
							writeDep( i, j, m, n, output.size(),  tt);
						else if(i == n && j != m)
							writeDep( i, j, m, n, output.size(),  ct);
						else if(i != n && j == m)
							writeDep( i, j, m, n, output.size(),  tc);
						else if(i == n && j == m)
							writeDep( i, j, m, n, output.size(),  in);
					}
				}
			}		
		}
	}

}

bool alex_command::execute()
{
	circuit circ = env->store<circuit>().current();
	matrix output;
	matrix qx4 ={{0,4,4,10,10},
			  {0,0,4,10,10},
			  {0,0,0,4,0},
			  {3,3,0,0,0},
			  {10,10,4,4,0}};

  	createMatrix( output, circ.lines() );
  	generateMatrixCnots( circ, output );
	printMatrixCnots( output );
	printObjectiveFunction( qx4, output, getNumberDifGates(output) );
	printFirstRestriction( qx4, output );

	getAllCombinations(output);

	printEndRestriction( qx4, output, getNumberDifGates( output ) );
	printIntegerVariables( qx4, output, getNumberDifGates( output ) );
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
