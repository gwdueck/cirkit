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

/**
 * @file move_qubit.hpp
 *
 * @brief Move a qubit from a to b (this is a partial swap)
 *
 * @author Gerhard Dueck
 * @since  2.3
 */

#ifndef MOVE_QUBIT_HPP
#define MOVE_QUBIT_HPP

#include <string>



namespace cirkit
{
    /* cab = move control from a to b given cnot(a,b)
     cba = move control from a to b given cnot(b,a)
     tab = move target from a to b given cnot(a,b)
     tba = move target from a to b given cnot(b,a)
     flip =
     (appended by i = inverse)
     */
     enum move_qubit_type { cab, cba, tab, tba, cabi, cbai, tabi, tbai, nop, flip, cnot3 };
    static move_qubit_type inverse_type[10] = {cabi, cbai, tabi, tbai, cab, cba, tab, tba, nop, flip };
        
    
    inline move_qubit_type invert_type( move_qubit_type a);

    class MoveQubit{
    private:
        
        static const std::string type_name[11];
        static const int move_cost[11];
        move_qubit_type mv_type;
        int v,w;
    public:
        MoveQubit() {};
        MoveQubit( move_qubit_type, int, int );
        void set( move_qubit_type, int, int );
        void print();
        int cost();
        void invert();
        move_qubit_type getType() { return mv_type; };
        unsigned int getA() { return v; };
        unsigned int getB() { return w; };
    };
    


}

#endif /* MOVE_QUBIT_HPP */

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
