// Copyright (C) 2013 Benjamin Kehlet
//
// This file is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2013-10-03
// Last changed: 2013-10-03

#include <intersection_segments.h>
#include <PolyhedronFactory.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Triangle_3 Triangle;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef typename Polyhedron::Facet_handle Facet_handle;
typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
typedef typename Kernel::Segment_3 Segment;

void two_tetrahedrons()
{
  Polyhedron a;

  make_tetrahedron(a, 
                   Point(1.0, 0.0, 0.0),
                   Point(2.0, 0.0, 0.0),
                   Point(1.5, 1.0, 0.0),
                   Point(1.5, .5, 10.0));

  Polyhedron b;
  make_tetrahedron(b,
                   Point(0.0, 0., .5),
                   Point(0.0, 0.0, 1.5),
                   Point(0.0, 1.0, 1.0),
                   Point(10.0, .5, 1.0));

  if (a.is_pure_triangle())
    std::cout << "a is pure triangle" << std::endl;

  if (b.is_pure_triangle())
    std::cout << "b is pure triangle" << std::endl;

  Polyhedron &biggest = a.size_of_facets() > b.size_of_facets() ? a : b;
  Polyhedron &smallest = a.size_of_facets() > b.size_of_facets() ? b : a;

  std::list<std::list<boost::tuple<Facet_handle, Facet_handle, Segment> > > polylines;
  {
    std::list<boost::tuple<Facet_handle, Facet_handle, Segment> > intersections;
    compute_intersections(biggest, smallest, std::back_inserter(intersections));

    for (std::list<boost::tuple<Facet_handle, Facet_handle, Segment> >::iterator it = intersections.begin();
         it != intersections.end(); it++)
    {
      {
        Halfedge_handle h = it->get<0>()->halfedge();
        Triangle t(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        assert(t.has_on(it->get<2>().source()));
        assert(t.has_on(it->get<2>().target()));
      }
      {
        Halfedge_handle h = it->get<1>()->halfedge();
        Triangle t(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        assert(t.has_on(it->get<2>().source()));
        assert(t.has_on(it->get<2>().target()));
      }
    }
    sort_polylines<Polyhedron>(biggest, smallest, intersections, polylines);
  }

  std::list<std::vector<typename Polyhedron::Halfedge_handle> > intersection_list;

  split_facets<Polyhedron, 0>(biggest,  polylines, intersection_list);
  //split_facets<Polyhedron, 1>(smallest, polylines);

}

void two_boxes()
{
  Polyhedron a;
  make_box(0,0,0, 4, 5, 2, a);

  Polyhedron b;
  make_box(1, 1, -1, 2, 2, 1, b);

  if (a.is_pure_triangle())
    std::cout << "a is pure triangle" << std::endl;

  if (b.is_pure_triangle())
    std::cout << "b is pure triangle" << std::endl;

  Polyhedron &biggest = a.size_of_facets() > b.size_of_facets() ? a : b;
  Polyhedron &smallest = a.size_of_facets() > b.size_of_facets() ? b : a;

  std::list<std::list<boost::tuple<Facet_handle, Facet_handle, Segment> > > polylines;
  {
    std::list<boost::tuple<Facet_handle, Facet_handle, Segment> > intersections;
    compute_intersections(biggest, smallest, std::back_inserter(intersections));

    for (std::list<boost::tuple<Facet_handle, Facet_handle, Segment> >::iterator it = intersections.begin();
         it != intersections.end(); it++)
    {
      {
        Halfedge_handle h = it->get<0>()->halfedge();
        Triangle t(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        assert(t.has_on(it->get<2>().source()));
        assert(t.has_on(it->get<2>().target()));
      }
      {
        Halfedge_handle h = it->get<1>()->halfedge();
        Triangle t(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        assert(t.has_on(it->get<2>().source()));
        assert(t.has_on(it->get<2>().target()));
      }
    }
    sort_polylines<Polyhedron>(biggest, smallest, intersections, polylines);
  }

  std::list<std::vector<Halfedge_handle> > a_edges;
  split_facets<Polyhedron, 0>(biggest, polylines, a_edges);
  check_splitting<Polyhedron, 0>(biggest, polylines, a_edges);
  //split_facets<Polyhedron, 1>(smallest, /* smallest, */ polylines);
}

int main(int argc, char** argv)
{
  two_boxes();
  return 0;
}
