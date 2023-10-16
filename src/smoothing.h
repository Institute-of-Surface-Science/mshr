// Copyright (C) 2016 Benjamin Kehlet
//
// This file is part of mshr.
//
// mshr is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mshr is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mshr. If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SMOOTHING_H
#define __SMOOTHING_H

namespace mshr
{

class LaplacianSmoothing
{
 public:
  template<typename Polyhedron>
  static void smooth(Polyhedron& p, double c=1.0)
  {
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    typedef typename Polyhedron::Halfedge_around_vertex_const_circulator HV_const_circulator;

    assert(p.is_valid());

    std::vector<std::pair<Vertex_handle, Point_3>> smoothed;

    for (typename Polyhedron::Vertex_iterator vit = p.vertices_begin();
         vit != p.vertices_end();
         vit++)
    {
      const Point_3 current = vit->point();

      Vector_3 delta;

      const HV_const_circulator h_start = vit->vertex_begin();
      HV_const_circulator h_current = h_start;
      do
      {
        delta = delta + (h_current->opposite()->vertex()->point()-current);
        h_current++;
      } while (h_current != h_start);

      Point_3 p = current + c*delta/vit->vertex_degree();

      // Evaluate exact value to reduce memory usage
      p.exact();

      smoothed.push_back(std::make_pair(Vertex_handle(vit), p));
    }

    // Apply the smoothing
    for (const std::pair<Vertex_handle, Point_3>& s : smoothed)
    {
      s.first->point() = s.second;
    }
  }
};

}
#endif
