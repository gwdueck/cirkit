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

#include "lpqx.hpp"

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
#include <core/utils/program_options.hpp>


namespace cirkit
{

using boost::program_options::value;

lpqx_command::lpqx_command( const environment::ptr& env )
	: cirkit_command( env, "Linear Programming to find the best mapping for IBM QX architecture" )
{
	opts.add_options()
	( "filename,f",    value( &filename ),  "name of the output file" )
	( "cplex,c",  					    	"write in cplex lp format (lp_solve is default)" )
	( "architecture,a", value_with_default( &architecture ) ,"select architecture\n" 
															"4: qx4 (5  qubits) -> default)\n"
															"2: qx2 (5  qubits)\n"
															"5: qx5 (16 qubits)" )
	( "version,v",value_with_default( &version ), "write LP in different versions\n" 
												  "1: write the gate version (default)\n"
												  "2: write the qubit version")
	;

}

command::rules_t lpqx_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

std::ofstream outputFile;
unsigned architecture;
bool cplex = false;
using matrix = std::vector<std::vector<unsigned>>;

// Create a matrix with the cnots 
void generateMatrixCnots( circuit& circ, matrix& m )
{
	// std::cout << "Generating matrix..." << std::endl;	
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
	// std::cout << "Printing matrix..." << std::endl;
	for (int i = 0; i < m.size(); ++i)
	{
		for (int j = 0; j < m[i].size(); ++j)
			std::cout << "\t" << m[i][j];
		std::cout << std::endl;
	}
	std::cout << std::endl;
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
void printObjectiveFunction( matrix& qx, matrix& cnots, unsigned difGates )
{
	unsigned aux = 0;
	if(!cplex)
	{
		outputFile << "/* Begin Objective Function */" << std::endl;
		outputFile << "min:\t";
	}
	else
		outputFile << "Minimize" << std::endl;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1 && !cplex)
							outputFile << qx[k][m]*cnots[i][j] << "G" << i << "_" << j << "c" << k << "_" << m << ";";
						else if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1 && cplex)
							outputFile << qx[k][m]*cnots[i][j] << "G" << i << "_" << j << "c" << k << "_" << m;
						else
							outputFile << qx[k][m]*cnots[i][j] << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						line = true;
					}
				}
			}
			if(line && aux < difGates-1)
			{
				outputFile << std::endl;
				++aux;
			}
			else if(line)
			{
				outputFile << std::endl;
			}
		}
	}
	if(cplex)
		outputFile << "st" << std::endl;
	else
		outputFile << "/* End Objective Function */" << std::endl;
}

// Function to print the one gate restriction
void printOneGateRestriction( matrix& cnots )
{
	unsigned aux = 0;
	if(!cplex)
		outputFile << "/* Begin One Gate Restriction */" << std::endl;
	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < cnots.size(); ++k)
			{
				for (int m = 0; m < cnots.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						if(k == cnots.size()-1 && m == cnots.size()-2 && !cplex)
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " = " << "1" << ";";
						else if(k == cnots.size()-1 && m == cnots.size()-2 && cplex)
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " = " << "1";
						else
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						line = true;
					}
				}
			}
			if(line)
				outputFile << std::endl;
		}
	}
	if(!cplex)
		outputFile << "/* End One Gate Restriction */" << std::endl;
}

// Function to print the final restriction
void printEndRestriction( matrix& qx, matrix& cnots, unsigned difGates )
{
	unsigned aux = 0;
	if(!cplex)
		outputFile << "/* Begin Final Restriction */" << std::endl;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1 && !cplex)
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " = " << difGates << ";";
						else if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1 && cplex)
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " = " << difGates;
						else
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						line = true;
					}
				}
			}
			if(line)
			{
				outputFile << std::endl;
				++aux;
			}
		}
	}
	if(!cplex)
		outputFile << "/* End Final Restriction */" << std::endl;
}

// Function to print the variables limits
void printLimitVariables( matrix& qx, matrix& cnots, unsigned difGates )
{
	unsigned aux = 0;
	if(!cplex)
		outputFile << "/* Begin Limit Variables */" << std::endl;
	else
		outputFile << "Bounds" << std::endl;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						if(cplex)
							outputFile << "0 <= G" << i << "_" << j << "c" << k << "_" << m << " <= 1" << std::endl;
						else
							outputFile << "0 <= G" << i << "_" << j << "c" << k << "_" << m << " <= 1;" << std::endl;
					}
				}
			}
		}
	}
	if(!cplex)
		outputFile << "/* End Limit Variables */" << std::endl;
}

// Function to print the variables
void printIntegerVariables( matrix& qx, matrix& cnots, unsigned difGates )
{
	unsigned aux = 0;
	if(!cplex)
	{
		outputFile << "/* Begin Integer Variables */" << std::endl;
		outputFile << "int\t";
	}
	else
		outputFile << "General" << std::endl;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						if(k == qx.size()-1 && m == qx.size()-2 && aux == difGates-1 && !cplex)
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << ";";
						else
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << "  ";
						line = true;
					}
				}
			}
			if(line && aux < difGates-1)
			{
				outputFile << std::endl;
				++aux;
			}
			else if(line)
			{
				outputFile << std::endl;
			}
		}
	}
	if(!cplex)
		outputFile << "/* End Integer Variables */" << std::endl;
	else
		outputFile << "End" << std::endl;
}

// Create a matrix with 0's
void createMatrix( matrix& m, unsigned size )
{
	// std::cout << "Creating matrix..." << std::endl;
  	std::vector<unsigned> v;
	for (int i = 0; i < size; ++i)
		v.push_back(0);
	for (int i = 0; i < size; ++i)
		m.push_back(v);
}



// Only function to write almost everything
// type = 0 -> control - control
// type = 1 -> target - target
// type = 2 -> control - target
// type = 3 -> target - control
// type = 4 -> inverse (control -> target; target -> control)
void writeDep( unsigned l0, unsigned c0, unsigned l1, unsigned c1, unsigned size, unsigned type )
{
	if(!cplex)
	{
		switch ( type )
	    {
	    	case 0:
				outputFile << "/* Writing Control - Control dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
	    		break;
	    	case 1:
				outputFile << "/* Writing Target - Target dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
	    		break;
	    	case 2:
				outputFile << "/* Writing Control - Target dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
	    		break;
	    	case 3:
				outputFile << "/* Writing Target - Control dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
	    		break;
	    	case 4:
				outputFile << "/* Writing Inverse dependency (" << l0 << "," << c0 << ")(" << l1 << "," << c1 << ") */" << std::endl;
	    		break;
	    	default:
				outputFile << "/* ERRORRRRRRRRRRRRRRRRR */" << std::endl;
	    }
	}
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			if(i != j && type != 4)
			{
				if(type == 0 || type == 2)
				{
					outputFile << "G" << l0 << "_" << c0;
					outputFile << "c" << i << "_" << j;
					if(cplex)
						outputFile << " - ";
					else
						outputFile << " <= ";
				}
				else
				{
					outputFile << "G" << l0 << "_" << c0;
					outputFile << "c" << j << "_" << i;
					if(!cplex)
						outputFile << " <= ";
					else
						outputFile << " - ";
				}
				unsigned aux = 0;
				for (int k = 0; k < size; ++k)
				{
					if(k != j && k != i)
					{
						++aux;
						if(type == 0 || type == 3)
						{
							outputFile << "G" << l1 << "_" << c1;
							outputFile << "c" << i << "_" << k;	
						}
						else
						{
							outputFile << "G" << l1 << "_" << c1;
							outputFile << "c" << k << "_" << i;
						}
						if(aux == size-2 && cplex)
							outputFile << " <= 0";
						else if(aux == size-2 && !cplex)
							outputFile << ";";						
						else if(!cplex)
							outputFile << " + ";	
						else	
							outputFile << " - ";	
					}
				}
				outputFile << std::endl;	
			}
			else if(i != j && type == 4)
			{
				outputFile << "G" << l0 << "_" << c0 << "c" << i << "_" << j;
				if(cplex)
					outputFile << " -G" << l1 << "_" << c1 << "c" << j << "_" << i << " <= 0\n";
				else
					outputFile << " <= G" << l1 << "_" << c1 << "c" << j << "_" << i << ";\n";
			}
		}
	}
}

bool empty_line(matrix& cnots, unsigned line)
{
	for (int i = 0; i < cnots.size(); ++i)
	{
		if( cnots[i][line] > 0 )
			return false;
		if( cnots[line][i] > 0 )
			return false;
	}
	return true;
}

// Create the ghost gates
void artificialGates(matrix& cnots)
{
	unsigned single, difGates;
	std::vector<int> ghost;
	bool inverse;

	bool qi, qj;
	std::pair<int, int> qa, qb;

    //going trough the cnots matrix
    for (int i = 0; i < cnots.size()-1; ++i)
	{
		for (int j = i+1; j < cnots.size(); ++j)
		{
			//"clear" the variables
			qi = qj = false;
			qa.first = qa.second = qb.first = qb.second = -1;
			//the cnot(i,j) is fixed... 
			//now the 'k' for will go trough the rest of the matrix in the lines and columns of 'i' and 'j'
			for (int k = 0; k < cnots.size(); ++k)
			{
				//if exists a cnot in (i,j) or (j,i), no need for a artificial gate
				if( cnots[i][j] > 0 || cnots[j][i] > 0 )
				{
					// std::cout << "Dont need to create ghost gate " << i << "-" << j << std::endl;
					break;
				}
				//check if exists a gate in the line and column of 'i'
				if( cnots[i][k] > 0 )
				{
					qi = true;
					qa.first = k;
				}
				else if( cnots[k][i] > 0 )
				{
					qi = true;
					qa.second = k;
				}
				//check if exists a gate in the line and column of 'j'
				if( cnots[k][j] > 0 )
				{
					qj = true;
					qb.first = k;
				}
				else if( cnots[j][k] > 0 )
				{
					qj = true;
					qb.second = k;
				}
			}
			if( qi && qj )
			{
				// std::cout << "Need artificial gate: " << i << "-" << j << std::endl;
				cnots[i][j] = 1;
			}
		}
	}
}

// Naive solution
void getAllCombinations(matrix cnots)
{
	// std::cout << "Getting all the dependencies..." << std::endl;
	enum type { cc, tt, ct, tc, in, sc, st };
	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots[i].size(); ++j)
		{
			bool done = false;
			for (int m = 0; m < cnots.size(); ++m)
			{
				for (int n = 0; n < cnots[m].size(); ++n)
				{
					if(i != j && m != n && cnots[i][j] > 0 && cnots[m][n] > 0)
					{
						if(i == m && j != n)
							writeDep( i, j, m, n, cnots.size(),  cc);
						else if(i != m && j == n)
							writeDep( i, j, m, n, cnots.size(),  tt);
						else if(i == n && j != m)
							writeDep( i, j, m, n, cnots.size(),  ct);
						else if(i != n && j == m)
							writeDep( i, j, m, n, cnots.size(),  tc);
						else if(i == n && j == m)
							writeDep( i, j, m, n, cnots.size(),  in);
						done = true;
					}
				}
			}
			if(done)
				cnots[i][j] = 0;		
		}
	}
	// std::cout << "Done!" << std::endl;
}

// Better approach
void getCombinationAnotherApproach(matrix& cnots)
{
	std::vector< std::pair< int,int > > q;
	if(!cplex)
	{
		for (int i = 0; i < cnots.size(); ++i)
		{
			for (int j = 0; j < cnots.size(); ++j)
			{
				if( cnots[i][j] > 0 )
					q.push_back(std::make_pair(i,j));
				if( cnots[j][i] > 0 )
					q.push_back(std::make_pair(j,i));
			}
			// std::cout << "tamanho " << i << " -> " << q.size() << std::endl;
			for (int m = 0; m < cnots.size(); ++m)
			{
				bool first = true;
				for (int j = 0; j < q.size(); ++j)
				{
					for (int n = 0; n < cnots.size(); ++n)
					{
						if( m != n )
						{
							if(first)
								outputFile << q.size()-1;
							outputFile << "G" << q[j].first << "_" << q[j].second; 
							if(q[j].first == i)
								outputFile << "c" << m << "_" << n;
							else
								outputFile << "c" << n << "_" << m;
							if(first && n == cnots.size()-1 || first && n == cnots.size()-2  && m == cnots.size()-1)
								outputFile << " = ";
							else
								outputFile << " + ";
						}
					}
					first = false;
				}
				if(!first)
					outputFile << "0;" << std::endl;
			}
			outputFile << std::endl;
			q.clear(); 
		}
	}
	else
	{
		for (int i = 0; i < cnots.size(); ++i)
		{
			for (int j = 0; j < cnots.size(); ++j)
			{
				if( cnots[i][j] > 0 )
					q.push_back(std::make_pair(i,j));
				if( cnots[j][i] > 0 )
					q.push_back(std::make_pair(j,i));
			}
			// std::cout << "tamanho " << i << " -> " << q.size() << std::endl;
			for (int m = 0; m < cnots.size(); ++m)
			{
				bool first = true;
				for (int j = 0; j < q.size(); ++j)
				{
					for (int n = 0; n < cnots.size(); ++n)
					{
						if( m != n )
						{
							if(first)
								outputFile << q.size()-1;
							outputFile << "G" << q[j].first << "_" << q[j].second; 
							if(q[j].first == i)
								outputFile << "c" << m << "_" << n;
							else
								outputFile << "c" << n << "_" << m;
							
							if(first && n == cnots.size()-1 || first && n == cnots.size()-2  && m == cnots.size()-1)
								outputFile << " - ";
							else if(first)
								outputFile << " + ";
							else if(m!=cnots.size()-1 && n == cnots.size()-1 && j == q.size()-1 || m==cnots.size()-1 && n == cnots.size()-2 && j == q.size()-1)
								outputFile << " ";
							else
								outputFile << " - ";
						}
					}
					first = false;
				}
				if(!first)
					outputFile << "= 0" << std::endl;
			}
			outputFile << std::endl;
			q.clear(); 
		}
	}
	
}

bool lpqx_command::execute()
{
	//circuit from store
	circuit circ = env->store<circuit>().current();
	//matrix of the cnots in the circuit
	matrix cnots;
	//matrix that assumes which architecture will be used
	matrix arch;
	//matrix of cost for ibm qx5
	matrix qx5 = {{0,4,10,20,19,29,39,51,61,64,54,42,30,20,10,4},
				{0,0,0,3,9,19,29,41,51,61,53,41,29,19,9,10},
				{10,4,0,0,3,22,32,44,54,64,54,42,30,20,3,4},
				{20,14,4,0,0,10,20,32,42,52,44,32,20,10,0,10},
				{30,24,14,4,0,4,14,20,30,40,32,20,14,4,10,20},
				{42,30,20,10,0,0,4,10,20,30,22,10,4,10,22,30},
				{54,42,32,22,3,0,0,0,10,20,3,0,10,22,34,42},
				{64,52,42,32,22,10,4,0,4,10,0,10,20,32,44,52},
				{76,64,54,44,34,22,10,0,0,4,3,20,30,42,54,64},
				{66,74,64,54,44,32,20,3,0,0,0,10,20,32,44,54},
				{54,62,52,42,32,20,14,4,10,4,0,4,14,20,32,42},
				{44,52,42,32,22,10,4,10,20,10,0,0,4,10,22,32},
				{34,42,32,22,3,0,10,22,32,22,3,0,0,0,3,22},
				{22,30,20,10,0,10,20,32,42,32,22,10,4,0,0,10},
				{10,20,10,4,10,20,30,42,52,42,32,20,14,4,0,4},
				{0,10,0,3,9,19,29,41,51,54,44,32,20,10,0,0}};
	//matrix of cost for ibm qx2 -- That is not qx2, has to update it.
	matrix qx2 ={{0,4,4,10,10},{0,0,4,10,10},{0,0,0,4,0},{3,3,0,0,0},{10,10,4,4,0}};
	//matrix qx4 ={{0,4,4,10,10},{0,0,4,10,10},{0,0,0,4,0},{3,3,0,0,0},{10,10,4,4,0}};
	//matrix of cost for ibm qx4
	matrix qx4 ={{0,4,4,10,10},
				{0,0,4,10,10},
				{0,0,0,4,0},
				{12,12,0,0,0},
				{10,10,4,4,0}};

	//update the variable arch
	switch ( architecture )
    {
    	case 2:
			arch = qx2;
    		break;
    	case 4:
			arch = qx4;
    		break;
    	case 5:
			arch = qx5;
    		break;
    	default:
    		std::cout << "Wrong architecture" << std::endl;
    		return true;
    }
    //check the number of qubits and the architecture selected
    if( arch.size() < circ.lines() )
	{
		std::cout << "This circuit requires an architecture with more qubits." << std::endl;
		return true;
	}
	//check if the output filename is missing
	if(filename.empty())
	{
		std::cout << "Missing output file. Use -f" << std::endl;
		return true;
	}
	//check version
	if(version != 1 && version != 2)
	{
		std::cout << "Wrong version" << std::endl;
		return true;
	}

	//check the file format
	if( is_set("cplex") )
		cplex = true;
	else
		cplex = false;
	
	outputFile.open (filename);
	
	//create the matrix of cnots of the circuit
  	createMatrix( cnots, circ.lines() );
  	
  	//fullfill the matrix with the cnots of the circuit
  	generateMatrixCnots( circ, cnots );
	// printMatrixCnots( cnots );

	//print the objective function for the LP
	printObjectiveFunction( arch, cnots, getNumberDifGates(cnots) );
	
	//if the circuit has only one gate, no need to get combinations
	if( getNumberDifGates(cnots) > 1 )
	{
		artificialGates(cnots);
		if( version == 1 )
			getAllCombinations(cnots);
		else if( version == 2 )
			getCombinationAnotherApproach(cnots);
	}
	//print the restriction that limits to one gate
	printOneGateRestriction( cnots );

	//print the type of variables of the LP
	printIntegerVariables( arch, cnots, getNumberDifGates( cnots ) );
	
	//clear the variables
  	outputFile.close();
  	filename.clear();

	return true;
}

command::log_opt_t lpqx_command::log() const
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
