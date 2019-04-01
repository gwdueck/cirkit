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
#include <reversible/pauli_tags.hpp>
#include <reversible/target_tags.hpp>

namespace cirkit
{

using boost::program_options::value;

lpqx_command::lpqx_command( const environment::ptr& env )
	: cirkit_command( env, "Linear Programming to find the best mapping for IBM QX architecture" )
{
	opts.add_options()
	( "filename,f",    value( &filename ),  "name of the output file" )
	( "lp_solve,l",  					    	"write in lp_solve format (cplex is default)" )
	( "toffoli,t",  					    	"Toffoli circuit" )
	( "architecture,a", value_with_default( &architecture ) ,"select architecture\n" 
															"4: qx4 (5  qubits) -> default)\n"
															"2: qx2 (5  qubits)\n"
															"5: qx5 (16 qubits)" )
	;

}

command::rules_t lpqx_command::validity_rules() const
{
  return {has_store_element<circuit>( env )};
}

std::ofstream outputFile;
unsigned architecture;
bool cplex = true;;
using matrix = std::vector<std::vector<unsigned>>;

// Create a matrix with the cnots 
void generateMatrixCnots( circuit& circ, matrix& m, matrix& v )
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
		  else if ( is_v( gate ) )
		  	++v[control][target];
		}
	}
}

void generateMatrixCnots( circuit& circ, matrix& m, matrix& v, matrix &t )
{
	// std::cout << "Generating matrix..." << std::endl;	
  	unsigned target, controla, controlb, control;
	for ( const auto& gate : circ )
	{
		if( gate.controls().size() == 1 )
		{
		  target = gate.targets().front();
		  control = gate.controls().front().line();
		  if( is_toffoli( gate ) )
		  	++m[control][target];
		  else if ( is_v( gate ) )
		  	++v[control][target];
		}
		else if( gate.controls().size() == 2 )
		{
		  target = gate.targets().front();
		  controla = gate.controls().front().line();
		  controlb = gate.controls().back().line();
		  ++t[controla+target][controlb];
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
void printObjectiveFunction( matrix& qx, matrix& cnots, matrix& vgates )
{
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " Begin Objective Function" << std::endl;

	if(cplex)
		outputFile << "Minimize" << std::endl;
	else
		outputFile << "min:\t" << std::endl;

	unsigned vdifgates = getNumberDifGates(vgates);
	unsigned cnotdifgates = getNumberDifGates(cnots);
	unsigned tam = (qx.size()*qx.size())-qx.size();
	unsigned aux = 0;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && vgates[i][j] > 0 && k != m)
					{
						++end;
						if( qx[k][m] < qx[m][k])
							outputFile << qx[k][m]*vgates[i][j]*2 << "V" << i << "_" << j << "c" << k << "_" << m;
						else
							outputFile << qx[m][k]*vgates[i][j]*2 << "V" << i << "_" << j << "c" << k << "_" << m;

						if(end < tam)
							outputFile << " + ";
						else
							++aux;
						line = true;
					}
				}
			}
			if( (line && aux < vdifgates) || (line && cnotdifgates > 0) )
				outputFile << " + " << std::endl;
		}
	}
	aux = 0;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						++end;
						outputFile << qx[k][m]*cnots[i][j] << "G" << i << "_" << j << "c" << k << "_" << m;
						if(end < tam)
							outputFile << " + ";
						else
							++aux;
						line = true;
					}
				}
			}
			if(line && aux < cnotdifgates)
				outputFile << " + " << std::endl;
		}
	}
	if(!cplex)
		outputFile << ";" << std::endl;
	else if(cplex)
		outputFile << "" << std::endl;

	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";
	
	outputFile << " End Objective Function" << std::endl;

	if(cplex)
		outputFile << "st" << std::endl;
}

// Function to print the objective function
void printObjectiveFunction( matrix& qx, matrix& cnots, matrix& vgates, matrix& tgates )
{
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " Begin Objective Function" << std::endl;

	if(cplex)
		outputFile << "Minimize" << std::endl;
	else
		outputFile << "min:\t" << std::endl;

	unsigned vdifgates = getNumberDifGates(vgates);
	unsigned cnotdifgates = getNumberDifGates(cnots);
	unsigned tam = (qx.size()*qx.size())-qx.size();
	unsigned aux = 0;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && vgates[i][j] > 0 && k != m)
					{
						++end;
						if( qx[k][m] < qx[m][k])
							outputFile << qx[k][m]*vgates[i][j]*2 << "V" << i << "_" << j << "c" << k << "_" << m;
						else
							outputFile << qx[m][k]*vgates[i][j]*2 << "V" << i << "_" << j << "c" << k << "_" << m;

						if(end < tam)
							outputFile << " + ";
						else
							++aux;
						line = true;
					}
				}
			}
			if( (line && aux < vdifgates) || (line && cnotdifgates > 0) )
				outputFile << " + " << std::endl;
		}
	}
	aux = 0;
	for (int i = 0; i < qx.size(); ++i)
	{
		for (int j = 0; j < qx.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < qx.size(); ++k)
			{
				for (int m = 0; m < qx.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						++end;
						outputFile << qx[k][m]*cnots[i][j] << "G" << i << "_" << j << "c" << k << "_" << m;
						if(end < tam)
							outputFile << " + ";
						else
							++aux;
						line = true;
					}
				}
			}
			if(line && aux < cnotdifgates)
				outputFile << " + " << std::endl;
		}
	}
	if(!cplex)
		outputFile << ";" << std::endl;
	else if(cplex)
		outputFile << "" << std::endl;

	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";
	
	outputFile << " End Objective Function" << std::endl;

	if(cplex)
		outputFile << "st" << std::endl;
}

// Function to print the one gate restriction
void printOneGateRestriction( matrix& cnots, matrix& vgates )
{
	unsigned aux = 0;
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " Begin One Gate Restriction" << std::endl;

	for (int i = 0; i < vgates.size(); ++i)
	{
		for (int j = 0; j < vgates.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < vgates.size(); ++k)
			{
				for (int m = 0; m < vgates.size(); ++m)
				{
					if( i != j && vgates[i][j] > 0 && k != m)
					{
						++end;
						if( end < (vgates.size()*vgates.size())-vgates.size() )
							outputFile << "V" << i << "_" << j << "c" << k << "_" << m << " + ";
						else
							outputFile << "V" << i << "_" << j << "c" << k << "_" << m;
						line = true;
					}
				}
			}
			if(line && cplex)
				outputFile << " = 1" << std::endl;
			else if(line && !cplex)
				outputFile << " = 1;" << std::endl;
		}
	}

	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < cnots.size(); ++k)
			{
				for (int m = 0; m < cnots.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						++end;
						if( end < (cnots.size()*cnots.size())-cnots.size() )
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						else
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m;
						line = true;
					}
				}
			}
			if(line && cplex)
				outputFile << " = 1" << std::endl;
			else if(line && !cplex)
				outputFile << " = 1;" << std::endl;
		}
	}
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";
	outputFile << " End One Gate Restriction" << std::endl;
}

// Function to print the one gate restriction
void printOneGateRestriction( matrix& cnots, matrix& vgates, matrix& tgates )
{
	unsigned aux = 0;
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " Begin One Gate Restriction" << std::endl;

	for (int i = 0; i < vgates.size(); ++i)
	{
		for (int j = 0; j < vgates.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < vgates.size(); ++k)
			{
				for (int m = 0; m < vgates.size(); ++m)
				{
					if( i != j && vgates[i][j] > 0 && k != m)
					{
						++end;
						if( end < (vgates.size()*vgates.size())-vgates.size() )
							outputFile << "V" << i << "_" << j << "c" << k << "_" << m << " + ";
						else
							outputFile << "V" << i << "_" << j << "c" << k << "_" << m;
						line = true;
					}
				}
			}
			if(line && cplex)
				outputFile << " = 1" << std::endl;
			else if(line && !cplex)
				outputFile << " = 1;" << std::endl;
		}
	}

	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < cnots.size(); ++k)
			{
				for (int m = 0; m < cnots.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						++end;
						if( end < (cnots.size()*cnots.size())-cnots.size() )
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m << " + ";
						else
							outputFile << "G" << i << "_" << j << "c" << k << "_" << m;
						line = true;
					}
				}
			}
			if(line && cplex)
				outputFile << " = 1" << std::endl;
			else if(line && !cplex)
				outputFile << " = 1;" << std::endl;
		}
	}
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";
	outputFile << " End One Gate Restriction" << std::endl;
}

// Function to print the final restriction
void printEndRestriction( matrix& qx, matrix& cnots, unsigned difGates )
{
	unsigned aux = 0;
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";
	outputFile << " Begin Final Restriction" << std::endl;

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
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";
	if(!cplex)
		outputFile << " End Final Restriction" << std::endl;
}

// Function to print the variables
void printIntegerVariables( matrix& cnots, matrix& vgates )
{
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " Begin Integer Variables" << std::endl;
	if(cplex)
		outputFile << "General" << std::endl;
	else
		outputFile << "int" << std::endl;

	unsigned vdifgates = getNumberDifGates(vgates);
	unsigned cnotdifgates = getNumberDifGates(cnots);
	unsigned tam = (cnots.size()*cnots.size())-cnots.size();
	unsigned aux = 0;
	for (int i = 0; i < vgates.size(); ++i)
	{
		for (int j = 0; j < vgates.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < vgates.size(); ++k)
			{
				for (int m = 0; m < vgates.size(); ++m)
				{
					if( i != j && vgates[i][j] > 0 && k != m)
					{
						++end;
						outputFile << "V" << i << "_" << j << "c" << k << "_" << m;
						if(end < tam)
							outputFile << " ";
						else
							++aux;
						line = true;
					}
				}
			}
			if( (line && aux < vdifgates) || (line && cnotdifgates > 0) )
				outputFile << " " << std::endl;
		}
	}
	aux = 0;
	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < cnots.size(); ++k)
			{
				for (int m = 0; m < cnots.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						++end;
						outputFile << "G" << i << "_" << j << "c" << k << "_" << m;
						if(end < tam)
							outputFile << " ";
						else
							++aux;
						line = true;
					}
				}
			}
			if(line && aux < cnotdifgates)
				outputFile << " " << std::endl;
		}
	}
	if(!cplex)
		outputFile << ";" << std::endl;
	else if(cplex)
		outputFile << "" << std::endl;

	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " End Integer Variables" << std::endl;

	if(cplex)
		outputFile << "End" << std::endl;
}

void printIntegerVariables( matrix& cnots, matrix& vgates, matrix& tgates )
{
	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " Begin Integer Variables" << std::endl;
	if(cplex)
		outputFile << "General" << std::endl;
	else
		outputFile << "int" << std::endl;

	unsigned vdifgates = getNumberDifGates(vgates);
	unsigned cnotdifgates = getNumberDifGates(cnots);
	unsigned tam = (cnots.size()*cnots.size())-cnots.size();
	unsigned aux = 0;
	for (int i = 0; i < vgates.size(); ++i)
	{
		for (int j = 0; j < vgates.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < vgates.size(); ++k)
			{
				for (int m = 0; m < vgates.size(); ++m)
				{
					if( i != j && vgates[i][j] > 0 && k != m)
					{
						++end;
						outputFile << "V" << i << "_" << j << "c" << k << "_" << m;
						if(end < tam)
							outputFile << " ";
						else
							++aux;
						line = true;
					}
				}
			}
			if( (line && aux < vdifgates) || (line && cnotdifgates > 0) )
				outputFile << " " << std::endl;
		}
	}
	aux = 0;
	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots.size(); ++j)
		{
			bool line = false;
			unsigned end = 0;
			for (int k = 0; k < cnots.size(); ++k)
			{
				for (int m = 0; m < cnots.size(); ++m)
				{
					if( i != j && cnots[i][j] > 0 && k != m)
					{
						++end;
						outputFile << "G" << i << "_" << j << "c" << k << "_" << m;
						if(end < tam)
							outputFile << " ";
						else
							++aux;
						line = true;
					}
				}
			}
			if(line && aux < cnotdifgates)
				outputFile << " " << std::endl;
		}
	}
	if(!cplex)
		outputFile << ";" << std::endl;
	else if(cplex)
		outputFile << "" << std::endl;

	if(cplex)
		outputFile << "\\";
	else
		outputFile << "//";

	outputFile << " End Integer Variables" << std::endl;

	if(cplex)
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

// Better approach
void getCombinationAnotherApproach(matrix& cnots, matrix& vgates)
{
	std::vector< std::pair< int,int > > q;
	std::vector< std::pair< int,int > > v;

	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots.size(); ++j)
		{
			if( cnots[i][j] > 0 )
				q.push_back(std::make_pair(i,j));
			if( cnots[j][i] > 0 )
				q.push_back(std::make_pair(j,i));
			if( vgates[i][j] > 0)
				v.push_back(std::make_pair(i,j));
			if( vgates[j][i] > 0 )
				v.push_back(std::make_pair(j,i));
		}

		bool first = true;
		unsigned signal = 0;
		if(v.size() + q.size() > 1)
		{
			for (int m = 0; m < cnots.size(); ++m)
			{
				first = true;
				signal = 0;
				for (int j = 0; j < v.size(); ++j)
				{
					for (int n = 0; n < cnots.size(); ++n)
					{
						if(m != n)
						{
							++signal;
							if(first)
								outputFile << q.size()+v.size()-1;
							outputFile << "V" << v[j].first << "_" << v[j].second; 
							if(v[j].first == i)
								outputFile << "c" << m << "_" << n;
							else
								outputFile << "c" << n << "_" << m;
							
							if( first && signal < cnots.size()-1 )
								outputFile << " + ";
							else if( q.size() == 0 && signal == (cnots.size()-1)*v.size() )
								outputFile << " ";
							else
								outputFile << " - ";
						}
					}
					if(v.size() > 0 && first)
						first = false;
				}
				for (int j = 0; j < q.size(); ++j)
				{
					for (int n = 0; n < cnots.size(); ++n)
					{
						if(m != n)
						{
							++signal;
							if(first)
								outputFile << q.size()+v.size()-1;
							outputFile << "G" << q[j].first << "_" << q[j].second; 
							if(q[j].first == i)
								outputFile << "c" << m << "_" << n;
							else
								outputFile << "c" << n << "_" << m;
							
							if( first && signal < cnots.size()-1 )
								outputFile << " + ";
							else if( signal == (cnots.size()-1)*(v.size()+q.size()) )
								outputFile << " ";
							else
								outputFile << " - ";
						}
					}
					if(q.size() > 0 && first)
						first = false;
				}
				if(!first && cplex)
					outputFile << "= 0" << std::endl;
				else if(!first && !cplex)
					outputFile << "= 0;" << std::endl;
			}
		}		
		if(!first)
			outputFile << std::endl;	
		q.clear();
		v.clear();
	}
}

void getCombinationAnotherApproach(matrix& cnots, matrix& vgates, matrix& tgates)
{
	std::vector< std::pair< int,int > > q;
	std::vector< std::pair< int,int > > v;

	for (int i = 0; i < cnots.size(); ++i)
	{
		for (int j = 0; j < cnots.size(); ++j)
		{
			if( cnots[i][j] > 0 )
				q.push_back(std::make_pair(i,j));
			if( cnots[j][i] > 0 )
				q.push_back(std::make_pair(j,i));
			if( vgates[i][j] > 0)
				v.push_back(std::make_pair(i,j));
			if( vgates[j][i] > 0 )
				v.push_back(std::make_pair(j,i));
		}

		bool first = true;
		unsigned signal = 0;
		if(v.size() + q.size() > 1)
		{
			for (int m = 0; m < cnots.size(); ++m)
			{
				first = true;
				signal = 0;
				for (int j = 0; j < v.size(); ++j)
				{
					for (int n = 0; n < cnots.size(); ++n)
					{
						if(m != n)
						{
							++signal;
							if(first)
								outputFile << q.size()+v.size()-1;
							outputFile << "V" << v[j].first << "_" << v[j].second; 
							if(v[j].first == i)
								outputFile << "c" << m << "_" << n;
							else
								outputFile << "c" << n << "_" << m;
							
							if( first && signal < cnots.size()-1 )
								outputFile << " + ";
							else if( q.size() == 0 && signal == (cnots.size()-1)*v.size() )
								outputFile << " ";
							else
								outputFile << " - ";
						}
					}
					if(v.size() > 0 && first)
						first = false;
				}
				for (int j = 0; j < q.size(); ++j)
				{
					for (int n = 0; n < cnots.size(); ++n)
					{
						if(m != n)
						{
							++signal;
							if(first)
								outputFile << q.size()+v.size()-1;
							outputFile << "G" << q[j].first << "_" << q[j].second; 
							if(q[j].first == i)
								outputFile << "c" << m << "_" << n;
							else
								outputFile << "c" << n << "_" << m;
							
							if( first && signal < cnots.size()-1 )
								outputFile << " + ";
							else if( signal == (cnots.size()-1)*(v.size()+q.size()) )
								outputFile << " ";
							else
								outputFile << " - ";
						}
					}
					if(q.size() > 0 && first)
						first = false;
				}
				if(!first && cplex)
					outputFile << "= 0" << std::endl;
				else if(!first && !cplex)
					outputFile << "= 0;" << std::endl;
			}
		}		
		if(!first)
			outputFile << std::endl;	
		q.clear();
		v.clear();
	}
}

void writeBlockRestrictions(matrix& res, unsigned s)
{
	for (int i = 0; i < s; ++i)
	{
		unsigned end = 0;
		for (int r = 0; r < res.size(); ++r)
		{
			for (int j = 0; j < s; ++j)
			{
				if(i != j)
				{
					++end;
					if(res[r][2] == 0)
						outputFile << "G" << res[r][0] << "_" << res[r][1] << "c" << i << "_" << j;
					else if(res[r][2] == 1)
						outputFile << "G" << res[r][0] << "_" << res[r][1] << "c" << j << "_" << i;
					else if(res[r][2] == 10)
						outputFile << "V" << res[r][0] << "_" << res[r][1] << "c" << i << "_" << j;
					else if(res[r][2] == 11)
						outputFile << "V" << res[r][0] << "_" << res[r][1] << "c" << j << "_" << i;
					
					if( end < (s-1)*res.size() )
						outputFile << " + ";
				}
			}	
		}
		if(cplex)
			outputFile << " <= 1" << std::endl;
		else
			outputFile << " <= 1;" << std::endl;
	}
}

void getBlockLessEqualRestrictions(matrix& cnots, matrix& vgates )
{
	matrix res;
	std::vector<unsigned> c;
	bool insert;
	for (int i = 0; i < cnots.size(); ++i)
	{
		insert = false;
		c.clear();
		for (int j = 0; j < cnots.size(); ++j)
		{
			if( cnots[i][j] > 0 )
			{
				c.push_back(i);
				c.push_back(j);
				c.push_back(0);
				insert = true;
				break;
			}
			else if( cnots[j][i] > 0 )
			{
				c.push_back(j);
				c.push_back(i);
				c.push_back(1);
				insert = true;
				break;
			}
			if( vgates[i][j] > 0 )
			{
				c.push_back(i);
				c.push_back(j);
				c.push_back(10);
				insert = true;
				break;
			}
			else if( vgates[j][i] > 0 )
			{
				c.push_back(j);
				c.push_back(i);
				c.push_back(11);
				insert = true;
				break;
			}
		}
		for (int k = 0; k < res.size(); ++k)
			if(c == res[k])
				insert = false;
		if(insert)
			res.push_back(c);
	}

	writeBlockRestrictions(res, cnots.size());
}

void getBlockLessEqualRestrictions(matrix& cnots, matrix& vgates, matrix& tgates )
{
	matrix res;
	std::vector<unsigned> c;
	bool insert;
	for (int i = 0; i < cnots.size(); ++i)
	{
		insert = false;
		c.clear();
		for (int j = 0; j < cnots.size(); ++j)
		{
			if( cnots[i][j] > 0 )
			{
				c.push_back(i);
				c.push_back(j);
				c.push_back(0);
				insert = true;
				break;
			}
			else if( cnots[j][i] > 0 )
			{
				c.push_back(j);
				c.push_back(i);
				c.push_back(1);
				insert = true;
				break;
			}
			if( vgates[i][j] > 0 )
			{
				c.push_back(i);
				c.push_back(j);
				c.push_back(10);
				insert = true;
				break;
			}
			else if( vgates[j][i] > 0 )
			{
				c.push_back(j);
				c.push_back(i);
				c.push_back(11);
				insert = true;
				break;
			}
		}
		for (int k = 0; k < res.size(); ++k)
			if(c == res[k])
				insert = false;
		if(insert)
			res.push_back(c);
	}

	writeBlockRestrictions(res, cnots.size());
}

bool lpqx_command::execute()
{
	bool artificial = false;
	//circuit from store
	circuit circ = env->store<circuit>().current();
	//matrix of the cnots in the circuit
	matrix cnots;
	//matrix of the v gates in the circuit
	matrix vgates;
	//matrix of the toffoli gates in the circuit
	matrix tgates;
	//matrix that assumes which architecture will be used
	matrix arch;
	//matrix of cost for ibm qx5

	matrix qx5 ={{0,4,10,13,19,29,39,47,61,55,41,41,25,20,10,4},
				{0,0,0,3,9,19,29,41,47,61,53,41,29,19,9,10},
				{10,4,0,0,3,13,23,35,41,55,41,35,23,13,3,4},
				{19,7,4,0,0,10,19,31,41,45,35,31,19,10,0,7},
				{25,13,7,4,0,4,7,19,25,33,23,19,7,4,10,13},
				{31,25,19,10,0,0,4,10,20,23,13,10,4,10,13,23},
				{45,33,23,13,3,0,0,0,10,13,3,0,10,13,25,33},
				{55,43,33,23,13,10,4,0,4,7,0,10,19,23,35,43},
				{67,51,45,31,25,22,10,0,0,4,3,13,23,31,41,51},
				{53,61,55,41,35,25,13,3,0,0,0,10,19,31,35,45},
				{45,53,43,33,23,19,7,4,7,4,0,4,7,19,23,33},
				{31,43,33,23,13,10,4,10,19,10,0,0,4,10,13,23},
				{25,33,23,13,3,0,10,13,23,13,3,0,0,0,3,13},
				{22,25,19,10,0,10,19,23,33,23,13,10,4,0,0,10},
				{10,13,7,4,10,19,25,33,43,33,23,19,7,4,0,4},
				{0,10,0,3,9,19,29,41,47,41,35,31,19,10,0,0}};
			
	//matrix of cost for ibm qx2 -- That is not qx2, has to update it.
	matrix qx2 ={{0,4,4,10,10},{0,0,4,10,10},{0,0,0,4,0},{3,3,0,0,0},{10,10,4,4,0}};
	//matrix of cost for ibm qx4
	matrix qx4 ={{0,4,4,7,7},{0,0,4,7,7},{0,0,0,4,4},{3,3,0,0,0},{3,3,0,4,0}};

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

	//check the file format
	if( is_set("lp_solve") )
		cplex = false;
	else
		cplex = true;

	if( is_set("artificial") )
		artificial = true;
	else
		artificial = false;
	
	outputFile.open (filename);
	
	//create the matrix of cnots of the circuit
  	createMatrix( cnots, circ.lines() );
  	createMatrix( vgates, circ.lines() );
	if( is_set("toffoli") )
  		createMatrix( tgates, circ.lines() * circ.lines() );
  	
  	//fullfill the matrix with the cnots and vgates of the circuit
	if( is_set("toffoli") )
  	  	generateMatrixCnots( circ, cnots, vgates, tgates );
  	else
  	  	generateMatrixCnots( circ, cnots, vgates );

	// printMatrixCnots( cnots );
	// std::cout << "v gates" << std::endl;
	// printMatrixCnots( vgates );

	//print the objective function for the LP
	if( is_set("toffoli") )
		printObjectiveFunction( arch, cnots, vgates, tgates );
	else
		printObjectiveFunction( arch, cnots, vgates );
	
	unsigned ndfg = 0;
	//if the circuit has only one gate, no need to get combinations
	if( is_set("toffoli") )
		ndfg = getNumberDifGates(cnots) + getNumberDifGates(vgates) + getNumberDifGates(tgates);
	else
		ndfg = getNumberDifGates(cnots) + getNumberDifGates(vgates);

	if( ndfg > 1 )
	{
		if( is_set("toffoli") )
		{
			getCombinationAnotherApproach(cnots, vgates, tgates);
			getBlockLessEqualRestrictions(cnots, vgates, tgates);
		}
		else
		{
			getCombinationAnotherApproach(cnots, vgates);
			getBlockLessEqualRestrictions(cnots, vgates);
		}
	}
	//print the restriction that limits to one gate
	if( is_set("toffoli") )
		printOneGateRestriction( cnots, vgates, tgates );
	else
		printOneGateRestriction( cnots, vgates );

	//print bounds
	// printBounds( arch, cnots );

	//print the type of variables of the LP
	if( is_set("toffoli") )
		printIntegerVariables( cnots, vgates, tgates );
	else
		printIntegerVariables( cnots, vgates );
	
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
