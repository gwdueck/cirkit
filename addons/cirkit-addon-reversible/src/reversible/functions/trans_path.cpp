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

#include <iostream>
#include <cmath>
#include "trans_path.hpp"
#include "move_qubit.hpp"

namespace cirkit
{
   
    void TransPath::add( MoveQubit q ){
        tpath.push_back( q );
    }
    void TransPath::print(){
        for ( auto &p : tpath ) {
            p.print();
        }
        std::cout << "cost = " << cost() << std::endl;
    }

    // Print to a file
    void TransPath::print( std::ofstream& graphfile ){
        for ( auto &p : tpath ) {
            p.print( graphfile );
        }
        graphfile << "cost = " << cost() << std::endl;
    }

    // Check if the movement has a cnot3 and update it
    void TransPath::movCnot3(){
        for (int i = tpath.size()-1; i > 0; --i)
        {
            if( tpath[i].getType() == nop && tpath[i-1].getType() == cab ) // two consecutive arrows
            {
                unsigned a = tpath[i-1].getA();  
                unsigned b = tpath[i-1].getB();  
                unsigned c = tpath[i].getB();  
                tpath.erase(tpath.begin()+i-1);                              
                tpath.erase(tpath.begin()+i-1);                                                 
                tpath.insert(tpath.begin()+i-1, MoveQubit( cnot3, a, b, c ) ); 
            }
            else if( tpath[i].getType() == flip && tpath[i-1].getType() == cba ) // two consecutive arrows backwards
            {
                unsigned a = tpath[i-1].getA();  
                unsigned b = tpath[i-1].getB();  
                unsigned c = tpath[i].getB();  
                tpath.erase(tpath.begin()+i-1);                                                 
                tpath.erase(tpath.begin()+i-1);                                                 
                tpath.insert(tpath.begin()+i-1, MoveQubit( cnot3i, a, b, c ) ); 
            }
            else if( tpath[i].getType() == cnot3 && tpath[i-1].getType() == cab ) //cnot3 with more arrows
            {
                unsigned a = tpath[i-1].getA();  
                unsigned b = tpath[i-1].getB(); 
                unsigned c = tpath[i].getB();  
                tpath.erase(tpath.begin()+i-1);                              
                tpath.insert(tpath.begin()+i-1, MoveQubit( cnot3, a, b, c ) );
            }
            else if( tpath[i].getType() == cnot3i && tpath[i-1].getType() == cba ) //cnot3 with more arrows
            {
                unsigned a = tpath[i-1].getA();  
                unsigned b = tpath[i-1].getB(); 
                unsigned c = tpath[i].getB();  
                tpath.erase(tpath.begin()+i-1);                              
                tpath.insert(tpath.begin()+i-1, MoveQubit( cnot3i, a, b, c ) );
            }
        }
    }

    // Possible optimization in the transformations - two adjacents hadamards
    int TransPath::opt(){
        int total = 0;
        for (int i = 0; i < tpath.size()/2; i=i+2){
            if( tpath[i].getType() == cba && tpath[i+1].getType() == cab) 
                total += 4;
            else if( tpath[i].getType()== cba && tpath[i+1].getType() == tba) 
                total += 4;
            else if( tpath[i].getType() == cba && tpath[i+1].getType() == flip)
                total += 4;
            else if( tpath[i].getType() == tab && tpath[i+1].getType() == cab)
                total += 4;
            else if( tpath[i].getType() == tab && tpath[i+1].getType() == tba)
                total += 4;
            else if( tpath[i].getType() == tab && tpath[i+1].getType() == flip)
                total += 4;
        }
        return total;
    }

    // Calculate the cost of the cnot3 movement
    int TransPath::cnot3Cost(){
        int res = 0, resi = 0;
        for ( auto &p : tpath ) {
            if(p.getType() == cnot3)
                ++res;
            if(p.getType() == cnot3i)
                ++resi;
        }
        if(res > 0 && resi == 0) // Formula: 2^n + 2^(n+1) - 2
            return (pow(2,res)+pow(2,++res)-2);
        else if(res > 0 && resi > 0)
            return ((pow(2,res)+pow(2,++res)-2) + (pow(2,resi)+pow(2,++resi)+2));
        else if(res == 0 && resi > 0)
            return (pow(2,resi)+pow(2,++resi)+2);
        else
            return 1;
    }

    int TransPath::cost(){
        int res = 0;
        for ( auto &p : tpath ) {
            res += p.cost();
        }
        return res + ( cnot3Cost() - 1 ); 
    }
    
    // add the cost of the inverse path
    int TransPath::costPlus(){
        int res = 0;
        for ( auto &p : tpath ) {
            res += p.cost();
        }
        return 2*res - tpath.back().cost() + ( cnot3Cost() - 1 );
    }
    
    void TransPath::addInverse(){
        MoveQubit q;
        for(int i = tpath.size() - 2; i >= 0; i-- )
        {
            if(tpath[i].getType() != cnot3 && tpath[i].getType() != cnot3i) // cnot3 has no inverse
            {
                q = tpath[i];
                q.invert();
                tpath.push_back( q );    
            }
        }
    }
}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
