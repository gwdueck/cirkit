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

#include "rotation_tags.hpp"

#include <reversible/target_tags.hpp>

namespace cirkit
{

bool is_rotation( const gate& g )
{
  return is_type<rotation_tag>( g.type() );
}

gate& create_rotation( gate& g, unsigned target, rotation_axis axis, double rotation )
{
  g.add_target( target );
  g.set_type( rotation_tag( axis, rotation ) );
  return g;
}

gate& append_rotation( circuit& circ, unsigned target, rotation_axis axis, double rotation )
{
  return create_rotation( circ.append_gate(), target, axis, rotation );
}

gate& prepend_rotation( circuit& circ, unsigned target, rotation_axis axis, double rotation )
{
  return create_rotation( circ.prepend_gate(), target, axis, rotation );
}

gate& insert_rotation( circuit& circ, unsigned n, unsigned target, rotation_axis axis, double rotation )
{
  return create_rotation( circ.insert_gate( n ), target, axis, rotation );
}

}

// Local Variables:
// c-basic-offset: 2
// eval: (c-set-offset 'substatement-open 0)
// eval: (c-set-offset 'innamespace 0)
// End:
