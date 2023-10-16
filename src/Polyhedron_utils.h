// Copyright (C) 2014-2015 Benjamin Kehlet
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
// along with mshr.  If not, see <http://www.gnu.org/licenses/>.


#ifndef POLYHEDRON_UTILS_H__
#define POLYHEDRON_UTILS_H__

#include <dolfin/math/basic.h>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
//#include <CGAL/Self_intersection_polyhedron_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/corefinement_operations.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <cmath>
#include <deque>
#include <fstream>

namespace mshr
{

class PolyhedronUtils
{
 public:

  //-----------------------------------------------------------------------------
  template<typename CDT>
  static void dump_2D_triangulation(const CDT& cdt, std::string filename)
  {
    std::cout << "Dumping 2D triangulation" << std::endl;

    // Count valid cells and connected vertices
    std::size_t num_cells = 0;
    std::map<typename CDT::Vertex_handle, std::size_t> vertex_map;

    for (typename CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
         cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
    {
      if (cgal_cell->is_in_domain())
      {
        num_cells++;

        for (std::size_t i = 0; i < 3; i++)
        {
          const typename CDT::Vertex_handle v = cgal_cell->vertex(i);
          if (vertex_map.count(v) == 0)
          {
            const std::size_t s = vertex_map.size();
            vertex_map[v] = s;
          }
        }
      }
    }

    std::cout << "Adding " << vertex_map.size() << " vertices and " << num_cells << " cells" << std::endl;

    std::ofstream outfile(filename);
    outfile << std::setprecision(16);
    outfile << "OFF" << std::endl;
    outfile << vertex_map.size() << " " << num_cells << " 0" << std::endl;
    outfile << std::endl;

    std::vector<typename CDT::Vertex_handle> vertices(vertex_map.size());
    for (const std::pair<typename CDT::Vertex_handle, std::size_t>& vertex : vertex_map)
      vertices[vertex.second] = vertex.first;

    for (const typename CDT::Vertex_handle v : vertices)
    {
      outfile << v->point()[0] << " "
              << v->point()[1] << " 0 " << std::endl;
    }

    for (typename CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
         cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
    {
      // Add cell if it is in the domain
      if (cgal_cell->is_in_domain())
      {
        outfile << "3 "
                << vertex_map[cgal_cell->vertex(0)] << " "
                << vertex_map[cgal_cell->vertex(1)] << " "
                << vertex_map[cgal_cell->vertex(2)] << std::endl;
      }
    }
  }
  //-----------------------------------------------------------------------------
  // Scans the vertices of the polyhedron the polyhedron and returns a
  // Polyhedron::Vertex_const_handle for each disconnected component.
  template <typename Polyhedron, typename OutputIterator>
  static void get_disconnected_components(const Polyhedron& p, OutputIterator it)
  {
    //typedef Polyhedron Polyhedron_t;
    typedef typename Polyhedron::Halfedge_around_vertex_const_circulator HV_const_circulator;
    typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;

    // store all vertices in a set
    std::set<Vertex_const_handle> v;
    for (typename Polyhedron::Vertex_const_iterator vit = p.vertices_begin();
         vit != p.vertices_end(); vit++)
      v.insert(vit);

    while (!v.empty())
    {
      // Add the component to the output
      typename std::set<Vertex_const_handle>::iterator start_it = v.begin();
      Vertex_const_handle start = *start_it;

      *it = start;
      it++;

      // Remove all vertices belonging to component from v
      std::deque<Vertex_const_handle> queue;
      queue.push_back(start);
      while (!queue.empty())
      {
        const Vertex_const_handle current = queue.front();
        queue.pop_front();

        if (v.count(current) > 0)
        {
          v.erase(current);

          const HV_const_circulator h_start = current->vertex_begin();
          HV_const_circulator h_current = h_start;
          do
          {
            queue.push_back(h_current->opposite()->vertex());
            h_current++;
          } while (h_current != h_start);
        }
      }
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::vector<typename Polyhedron::Halfedge_handle>
  get_holes(Polyhedron& P)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    P.normalize_border();

    std::set<Halfedge_handle> border_edges;
    for (typename Polyhedron::Halfedge_iterator hit = P.border_halfedges_begin(); hit != P.halfedges_end(); hit++)
    {
      Halfedge_handle b = hit;
      dolfin_assert(b->is_border_edge());
      if (b->is_border())
        border_edges.insert(b);
      else if (b->opposite()->is_border())
        border_edges.insert(b->opposite());
      else
        dolfin_assert(false);
    }

    std::vector<Halfedge_handle> border_begins;
    while (!border_edges.empty())
    {
      Halfedge_handle current = *(border_edges.begin());
      border_begins.push_back(current);

      std::size_t counter = 0;
      const Halfedge_handle start = current;
      do
      {
        dolfin_assert(border_edges.count(current) == 1);
        border_edges.erase(current);
        counter++;
        current = current->next();
      } while (current != start);
    }

    return std::move(border_begins);
  }
  //-----------------------------------------------------------------------------
  // Count the number of edges in a facet or along a hole
  template<typename Halfedge_handle>
  static std::size_t edge_count(Halfedge_handle h)
  {
    Halfedge_handle current;
    std::size_t counter = 0;
    do
    {
      counter++;
      current = current->next();
    } while (current != h);

    return counter;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Traits::Triangle_3 get_facet_triangle(typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    dolfin_assert(h->facet()->is_triangle());

    return Triangle_3(h->vertex()->point(),
                      h->next()->vertex()->point(),
                      h->next()->next()->vertex()->point());
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static double cos_min_angle(typename Polyhedron::Traits::Vector_3 v,
                                typename Polyhedron::Halfedge_handle h1,
                                typename Polyhedron::Halfedge_handle h2)
  {
    v /= v.squared_length();

    double max = 0;
    typename Polyhedron::Halfedge_handle current = h1;
    do
    {
      typename Polyhedron::Traits::Vector_3 w(current->vertex()->point(),
                                              current->next()->vertex()->point());
      w /= std::sqrt(CGAL::to_double(w.squared_length()));
      max = std::max(std::abs(CGAL::to_double(w*v)));
    } while (current != h2);

    return max;
  }
  //-----------------------------------------------------------------------------
  // Given a facet and a vertex (assumed to be incident to the facet), find the
  // corresponding halfedge
  template<typename Polyhedron>
  static typename Polyhedron::Halfedge_handle
  find_edge(typename Polyhedron::Vertex_handle v,
            typename Polyhedron::Face_handle f)
  {
    dolfin_assert(v != typename Polyhedron::Vertex_handle());
    dolfin_assert(f != typename Polyhedron::Face_handle());

    typename Polyhedron::Halfedge_around_vertex_circulator start = v->vertex_begin();
    typename Polyhedron::Halfedge_around_vertex_circulator current = start;

    do
    {
      if (current->facet() == f)
      {
        dolfin_assert(current->vertex() == v && current->facet()  == f);
        return current;
      }

      current++;
    } while(current != start);

    dolfin_assert(false);
    return typename Polyhedron::Halfedge_handle();
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Facet_handle
  find_common_facet(typename Polyhedron::Vertex_handle v0,
                    typename Polyhedron::Vertex_handle v1,
                    typename Polyhedron::Vertex_handle v2)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_circulator He_circulator;
    std::set<typename Polyhedron::Facet_handle> facets0;
    {
      He_circulator start = v0->vertex_begin();
      He_circulator current = start;
      do
      {
        if (!current->is_border())
          facets0.insert(current->facet());

        current++;
      } while (current != start);
    }

    std::set<typename Polyhedron::Facet_handle> facets1;
    {
      He_circulator start = v1->vertex_begin();
      He_circulator current = start;
      do
      {
        if (!current->is_border() && facets0.count(current->facet()) > 0)
          facets1.insert(current->facet());

        current++;
      } while (current != start);
    }

    std::set<typename Polyhedron::Facet_handle> facets2;
    {
      He_circulator start = v2->vertex_begin();
      He_circulator current = start;
      do
      {
        if (!current->is_border() && facets1.count(current->facet()) > 0)
          facets2.insert(current->facet());

        current++;
      } while (current != start);
    }

    dolfin_assert(facets2.size() < 2);
    dolfin_assert(facets2.size() > 0);

    return *(facets2.begin());
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void insert_edge(Polyhedron& P,
                          typename Polyhedron::Vertex_handle h,
                          typename Polyhedron::Vertex_handle g,
                          typename Polyhedron::Facet_handle f)
  {
    //std::pair<typename Polyhedron::Halfedge_handle, typename Polyhedron::Halfedge_handle>
    //edges = find_edges(h_edge, g_edge);
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    Halfedge_handle h_edge, g_edge;
    const Halfedge_handle start = g->halfedge();
    Halfedge_handle current = start;
    do
    {
      if (current->vertex() == h)
        h_edge = current;
      else if(current->vertex() == g)
        g_edge == current;
    } while (current != start);

    dolfin_assert(h_edge != Halfedge_handle() && g_edge != Halfedge_handle());
    P.split_facet(h_edge, g_edge);
    dolfin_assert(P.is_valid());
  }
  //-----------------------------------------------------------------------------
  template<typename Triangle_3>
  static double get_triangle_cos_angle(Triangle_3 t1,
                                       Triangle_3 t2)
  {
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Vector_3 Vector_3;

    const Vector_3 v1 = t1.supporting_plane().orthogonal_vector();
    const Vector_3 v2 = t2.supporting_plane().orthogonal_vector();

    return CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())));
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static double get_edge_cos_angle(typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;

    // std::cout << "get edge cos theta" << std::endl;

    Vector_3 h1_vec(h->vertex()->point(), h->prev()->vertex()->point());
    h1_vec = h1_vec/std::sqrt(CGAL::to_double(h1_vec.squared_length()));

    // std::cout << "h1_vec_normalized: " << h1_vec << std::endl;
    Vector_3 h2_vec(h->vertex()->point(), h->next()->vertex()->point());
    h2_vec = h2_vec/std::sqrt(CGAL::to_double(h2_vec.squared_length()));
    // std::cout << "h2_vec_normalized: " << h2_vec << std::endl;
    const double cos_theta = CGAL::to_double(h1_vec*h2_vec);
    // std::cout << "Cos theta: " << cos_theta << std::endl;
    return cos_theta;
  }
  //-----------------------------------------------------------------------------
  // Compute the plane fit quality of the vertices from h1 to h2 both included
  // Return fitting quality, max projection distance, max cos angle
  template<typename Polyhedron>
  static std::array<double, 3>
  get_plane_fit(const typename Polyhedron::Halfedge_handle h1,
                const typename Polyhedron::Halfedge_handle h2)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Plane_3 Plane_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    typedef typename Polyhedron::Traits::FT FT;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    //typedef typename InexactKernel::Segment_3 InexactSegment_3;
    //typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    // std::cout << "Get plane fit" << std::endl;

    // std::cout << "Polygon ";
    std::vector<InexactPoint_3> points;
    //std::vector<InexactSegment_3> segments;
    Halfedge_handle current = h1;
    do
    {
      const Point_3& p = current->vertex()->point();
      // std::cout << ", " << p;
      points.push_back(InexactPoint_3(CGAL::to_double(p[0]),
                                      CGAL::to_double(p[1]),
                                      CGAL::to_double(p[2])));
      current = current->next();
    } while (current != h2);

    {
      const Point_3& p = h2->vertex()->point();
      // std::cout << ", " << p;
      points.push_back(InexactPoint_3(CGAL::to_double(p[0]),
                                      CGAL::to_double(p[1]),
                                      CGAL::to_double(p[2])));
    }

    dolfin_assert(points.size() > 2);
    // std::cout << "Size: " << points.size() << std::endl;
    //std::cout << std::endl;
    InexactPlane_3 fitting_plane_inexact;
    const double fit_quality = CGAL::linear_least_squares_fitting_3(points.begin(),
                                                                    points.end(),
                                                                    fitting_plane_inexact,
                                                                    CGAL::Dimension_tag<0>());
    Plane_3 fitting_plane(fitting_plane_inexact.a(),
                          fitting_plane_inexact.b(),
                          fitting_plane_inexact.c(),
                          fitting_plane_inexact.d());
    // std::cout << "Plane: " << fitting_plane << std::endl;
    // std::cout << "Length of normal: " << fitting_plane.orthogonal_vector().squared_length() << std::endl;
    const Vector_3 normal = fitting_plane.orthogonal_vector()/std::sqrt(CGAL::to_double(fitting_plane.orthogonal_vector().squared_length()));

    FT max_distance = (h1->vertex()->point()-fitting_plane.projection(h1->vertex()->point())).squared_length();
    FT max_angle = 0;

    Halfedge_handle prev = h1;
    current = h1->next();
    do
    {
      const Vector_3 v = current->vertex()->point()-prev->vertex()->point();
      const FT cos_angle = v/std::sqrt(CGAL::to_double(v.squared_length())) * normal;
      const Point_3 projection = fitting_plane.projection(current->vertex()->point());
      max_angle = std::max(max_angle, cos_angle);
      max_distance = std::max(max_distance, (current->vertex()->point()-projection).squared_length());

      prev = current;
      current = current->next();
    } while (prev != h2);

    // std::cout << "Fit quality: " << fit_quality << ", max distance: " << max_distance << ", cos_angle: " << max_angle << std::endl;
    return std::array<double, 3>{fit_quality, CGAL::to_double(max_distance), CGAL::to_double(max_angle)};
    // return -max_distance;
    //return CGAL::to_double(fit_quality - max_angle);
  }
  //-----------------------------------------------------------------------------
  // Compute the fit quality heuristic of the vertices from h1 to h2 both
  // included.

  // Some experiences:
  // * The plane fit quality as returned from
  // CGAL::linear_least_squares_fitting_3() is not suitable here. It measures
  // how distinct the best fitting plane is, rather than the actual quality of
  // the fit.
  // * The max angle between the segments and the normal work good is the hole
  // has no more than one "kink".
  // * The max distance from the points to their projection is also a pretty
  // good measure, but needs to be normalized in some clever way if not used
  // solely.
  template<typename Polyhedron>
  static double evaluate_hole_subdivision(const typename Polyhedron::Halfedge_handle h1,
                                          const typename Polyhedron::Halfedge_handle h2)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef CGAL::Simple_cartesian<double> InexactKernel;
    //typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    typedef typename InexactKernel::Segment_3 InexactSegment_3;
    typedef typename InexactKernel::Vector_3 InexactVector_3;

    // double plane1fit;
    double max_angle1 = 0;
    InexactPlane_3 fitting_plane1;
    double avg_distance_squared1 = 0;
    {
      std::vector<InexactSegment_3> segments;
      Halfedge_handle current = h1;
      do
      {
        const Point_3& p = current->vertex()->point();
        const Point_3& next = current->next()->vertex()->point();
        segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));
        current = current->next();
      } while (current != h2);

      {
        const Point_3& p = h2->vertex()->point();
        const Point_3& next = h1->vertex()->point();
        segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));
      }

      dolfin_assert(segments.size() > 2);

      /* plane1fit = */ CGAL::linear_least_squares_fitting_3(segments.begin(),
                                                             segments.end(),
                                                             fitting_plane1,
                                                             CGAL::Dimension_tag<1>());

      // std::cout << "  Plane 1: " << fitting_plane1 << ", " << plane1fit << std::endl;

      // Compute max angle between segment and fitting plane
      const InexactVector_3 normal = fitting_plane1.orthogonal_vector();
      // std::cout << "Length: " << normal.squared_length() << std::endl;
      dolfin_assert(dolfin::near(normal.squared_length(), 1, DOLFIN_EPS_LARGE));
      for (auto sit = segments.begin(); sit != segments.end(); sit++)
      {
        const InexactVector_3 v = InexactVector_3(*sit) / std::sqrt(sit->squared_length());
        // std::cout << "Length: " << (v.squared_length()-1) << std::endl;
        dolfin_assert(dolfin::near(v.squared_length(), 1, DOLFIN_EPS_LARGE));
        max_angle1 = std::max(max_angle1, std::abs(v*normal));
        // max_distance_squared1 = std::max(max_distance_squared1,
        avg_distance_squared1 += InexactVector_3(sit->source(), fitting_plane1.projection(sit->source())).squared_length();
      }
      avg_distance_squared1 /= segments.size();
    }

    /* double plane2fit; */
    InexactPlane_3 fitting_plane2;
    double max_angle2 = 0;
    //double max_distance_squared2 = 0;
    double avg_distance_squared2 = 0;
    {
      std::vector<InexactSegment_3> segments;
      Halfedge_handle current = h2;
      do
      {
        const Point_3& p = current->vertex()->point();
        const Point_3& next = current->next()->vertex()->point();
        segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));
        current = current->next();
      } while (current != h1);

      const Point_3& p = h1->vertex()->point();
      const Point_3& next = h2->vertex()->point();

      segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                          InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));

      dolfin_assert(segments.size() > 2);
      // std::cout << "  Size: " << segments.size() << std::endl;

      /* plane2fit = */ CGAL::linear_least_squares_fitting_3(segments.begin(),
                                                             segments.end(),
                                                             fitting_plane2,
                                                             CGAL::Dimension_tag<1>());

      // std::cout << "  Plane 1: " << fitting_plane1 << ", " << plane1fit << std::endl;

      // Compute max angle between plane and segments
      const InexactVector_3 normal = fitting_plane2.orthogonal_vector();
      dolfin_assert(dolfin::near(normal.squared_length(), 1, DOLFIN_EPS_LARGE));
      for (auto sit = segments.begin(); sit != segments.end(); sit++)
      {
        const InexactVector_3 v = InexactVector_3(*sit) / std::sqrt(sit->squared_length());
        dolfin_assert(dolfin::near(v.squared_length(), 1, DOLFIN_EPS_LARGE));
        max_angle2 = std::max(max_angle2, std::abs(v*normal));
        /* max_distance_squared2 = std::max(max_distance_squared2, */
        /*                                  InexactVector_3(sit->source(), fitting_plane2.projection(sit->source())).squared_length()); */
        avg_distance_squared2 += InexactVector_3(sit->source(), fitting_plane2.projection(sit->source())).squared_length();
      }
      avg_distance_squared2 /= segments.size();
    }

    // const double cos_angle = fitting_plane1.orthogonal_vector()*fitting_plane2.orthogonal_vector();
    // std::cout << "  Angle: " << cos_angle << "(" << acos(cos_angle)/(2*DOLFIN_PI)*360 << ")" << std::endl;
    // std::cout << "Max angles: " << max_angle1 << "(" << acos(max_angle1)/(2*DOLFIN_PI)*360 << "), " << max_angle2 << " (" << acos(max_angle2)/(2*DOLFIN_PI)*360 << ")" << std::endl;
    // return std::min(plane1fit, plane2fit); // - .04*cos_angle - .05*max_angle1 - .05*max_angle2;// + triangulation_extra;

    //return (-avg_distance_squared1 - avg_distance_squared2)/std::max(max_angle1, max_angle2); ///std::min(plane1fit, plane2fit);
    return -std::max(max_angle1, max_angle2);
  }

  //-----------------------------------------------------------------------------
  // Compute the transformation that rotates a given vector (assumed to be of
  // unit length) into (0,0,1)
  template<typename Vector_3>
  static CGAL::Aff_transformation_3<typename CGAL::Kernel_traits<Vector_3>::Kernel>
    rotate_to_xy(Vector_3 a)
  {
    typedef typename CGAL::Kernel_traits<Vector_3>::Kernel::RT RT;
    typedef typename CGAL::Aff_transformation_3<typename CGAL::Kernel_traits<Vector_3>::Kernel> Aff_transformation_3;

    // Inner product of a and target vector (0,0,1)
    const RT cos_theta = a[2];

    // Cross product of a and target vector (0,0,1)
    // const RT sine_theta = CGAL::sqrt(a[1]*a[1] + a[0]*a[0]);
    const RT sine_theta = sqrt(CGAL::to_double(a[1]*a[1] + a[0]*a[0]));

    const RT ux = a[1]/sine_theta;;
    const RT uy = -a[0]/sine_theta;
    const RT uz = 0;
    dolfin_assert(CGAL::abs(ux*ux + uy*uy + uz*uz) - 1 < DOLFIN_EPS);

    return Aff_transformation_3(
      cos_theta+ux*ux*(1-cos_theta),     ux*uy*(1-cos_theta)-uz*sine_theta, ux*uz*(1-cos_theta)+uy*sine_theta,
      uy*ux*(1-cos_theta)+uz*sine_theta, cos_theta+uy*uy*(1-cos_theta),     uy*uz*(1-cos_theta)-ux*sine_theta,
      uz*ux*(1-cos_theta)-uy*sine_theta, uz*uy*(1-cos_theta)+ux*sine_theta, cos_theta+uz*uz*(1-cos_theta));
  }
  //-----------------------------------------------------------------------------
  template <typename HDS, typename CDT>
  class Add2DTriangulation : public CGAL::Modifier_base<HDS>
  {
  public:
 Add2DTriangulation(const CDT& cdt,
                    const typename HDS::Traits::Aff_transformation_3& from_xy,
                    const typename CDT::Vertex_handle v0,
                    const typename CDT::Vertex_handle v1,
                    const typename HDS::Traits::Vector_3 displacement=typename HDS::Traits::Vector_3(CGAL::Null_vector()))
   : cdt(cdt),
      rotate_from_xy(from_xy),
      displacement(displacement),
      v0_2d(v0),
      v1_2d(v1),
      v0_initialized(false),
      v1_initialized(false)
    { }

    void operator()(HDS& hds)
    {
      CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);

      for (typename CDT::Finite_vertices_iterator cgal_vertex = cdt.finite_vertices_begin();
           cgal_vertex != cdt.finite_vertices_end(); ++cgal_vertex)
      {
        cgal_vertex->info() = std::make_pair(typename HDS::Vertex_handle(),
                                             std::numeric_limits<std::size_t>::max());
      }

      // Find the face incident to both v0 and v1
      // We will use this to check the orientation of the faces when they are
      // inserted into the 3D polyhedron
      typename CDT::Face_handle f;

      // Count valid cells and connected vertices
      std::size_t num_cells = 0;
      std::size_t num_vertices = 0;
      for (typename CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
           cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
      {
        if (cgal_cell->is_in_domain())
        {
          num_cells++;

          for (std::size_t i = 0; i < 3; i++)
          {
            typename CDT::Vertex_handle v = cgal_cell->vertex(i);
            if (v->info().second == std::numeric_limits<std::size_t>::max())
            {
              v->info().second = num_vertices;
              num_vertices++;
            }
          }

          if (cgal_cell->has_vertex(v0_2d) && cgal_cell->has_vertex(v1_2d))
          {
            std::cout << "!!!! Found the face!!!" << std::endl;
            std::cout << cgal_cell->index(v0_2d) << " " << cgal_cell->index(v1_2d) << std::endl;
            f = cgal_cell;
          }
        }
      }

      // Check the orientation of the facets
      const bool flip = ((f->index(v1_2d)+1)%3 == f->index(v0_2d));

      const typename HDS::Traits::Vector_3 displacement_flipped = flip ? -displacement : displacement;

      builder.begin_surface(num_vertices, num_cells);

      std::cout << "Adding " << num_vertices << " vertices and " << num_cells << " cells" << std::endl;

      // Add vertices
      std::size_t vertex_index = 0;
      for (typename CDT::Finite_vertices_iterator cgal_vertex = cdt.finite_vertices_begin();
           cgal_vertex != cdt.finite_vertices_end(); ++cgal_vertex)
      {
        const typename CDT::Vertex_handle current = cgal_vertex;

        // Transform point from xy plane to where it belongs in the polyhendron (2D point is EPICK)
        typename HDS::Traits::Point_3 p = rotate_from_xy(typename HDS::Traits::Point_3(cgal_vertex->point()[0],
                                                                                       cgal_vertex->point()[1],
                                                                                       0))+displacement;

        /* std::cout << "  2D point: " << cgal_vertex->point() << " (" << z << ")" << std::endl; */
        /* std::cout << "  Rotated:  " << p << std::endl; */
        /* std::cout << "  Index: " << vertex_index << std::endl; */
// FRom remove-null-facets
          /* // std::cout << "Now creating face" << std::endl; */

          /* // This typedef for some reason gives a "unused local typedef" warning */
          /* //from clang (...?) --> write out the typename in the statements */
          /* //below. */
          /* // typedef typename Halfedge::Base HBase; */
          /* edges[0]->Halfedge::Base::set_next(edges[1]); */
          /* decorator.set_prev(edges[1], edges[0]); */
          /* edges[1]->Halfedge::Base::set_next(edges[2]); */
          /* decorator.set_prev(edges[2], edges[1]); */
          /* edges[2]->Halfedge::Base::set_next(edges[0]); */
          /* decorator.set_prev(edges[0], edges[2]); */
// end remove-null-facets

        // Add vertex (convert point to EPECK)
        typename HDS::Vertex_handle h = builder.add_vertex(p);
        if (current == v0_2d)
        {
          v0 = h;
          v0_initialized = true;
        }
        else if (current == v1_2d)
        {
          v1 = h;
          v1_initialized = true;
        }

        // Attach index to vertex and increment
        cgal_vertex->info().second = vertex_index++;
      }

      // Add cells to mesh and build domain marker mesh function
      // std::size_t cell_index = 0;
      for (typename CDT::Finite_faces_iterator cgal_cell = cdt.finite_faces_begin();
           cgal_cell != cdt.finite_faces_end(); ++cgal_cell)
      {
        // Add cell if it is in the domain
        if (cgal_cell->is_in_domain())
        {
          const typename CDT::Geom_traits::Triangle_2 t_2d(cgal_cell->vertex(0)->point(),
                                                           cgal_cell->vertex(1)->point(),
                                                           cgal_cell->vertex(2)->point());

          // Transform point from xy plane to where it belongs in the polyhendron (2D point is EPICK)
          typename HDS::Traits::Point_3 p0 = rotate_from_xy(typename HDS::Traits::Point_3(cgal_cell->vertex(0)->point()[0],
                                                                                          cgal_cell->vertex(0)->point()[1],
                                                                                          0));

          typename HDS::Traits::Point_3 p1 = rotate_from_xy(typename HDS::Traits::Point_3(cgal_cell->vertex(1)->point()[0],
                                                                                          cgal_cell->vertex(1)->point()[1],
                                                                                          0));

          typename HDS::Traits::Point_3 p2 = rotate_from_xy(typename HDS::Traits::Point_3(cgal_cell->vertex(2)->point()[0],
                                                                                          cgal_cell->vertex(2)->point()[1],
                                                                                          0));

          typename HDS::Traits::Triangle_3 t_3d(p0, p1, p2);

          const double diff = sqrt(CGAL::to_double(t_3d.squared_area()))-CGAL::to_double(CGAL::abs(t_2d.area()));
          if (std::abs(diff) > 1e-10)
            std::cout << "Rotated triangle differs: " << diff << std::endl;



          builder.begin_facet();
          builder.add_vertex_to_facet(cgal_cell->vertex(0)->info().second);
          if (flip)
          {
            builder.add_vertex_to_facet(cgal_cell->vertex(2)->info().second);
            builder.add_vertex_to_facet(cgal_cell->vertex(1)->info().second);
            /* std::cout << "  Adding vertex: (" << cgal_cell->vertex(0)->info().second */
            /*           << ", " << cgal_cell->vertex(2)->info().second */
            /*           << ", " << cgal_cell->vertex(1)->info().second << std::endl; */
          }
          else
          {
            builder.add_vertex_to_facet(cgal_cell->vertex(1)->info().second);
            builder.add_vertex_to_facet(cgal_cell->vertex(2)->info().second);
            /* std::cout << "  Adding vertex: (" << cgal_cell->vertex(0)->info().second */
            /*           << ", " << cgal_cell->vertex(1)->info().second */
            /*           << ", " << cgal_cell->vertex(2)->info().second << std::endl; */

          }
          new_facets.insert(builder.end_facet());
        }
      }

      builder.end_surface();

      dolfin_assert(v0_initialized);
      dolfin_assert(v1_initialized);
    }
    const CDT& cdt;
    const typename HDS::Traits::Aff_transformation_3 rotate_from_xy;
    const typename HDS::Traits::Vector_3 displacement;
    const typename CDT::Vertex_handle v0_2d, v1_2d;
    typename HDS::Vertex_handle v0, v1;
    bool v0_initialized, v1_initialized;
    std::set<typename HDS::Halfedge_handle> new_facets;
  };

  /// Attempts to triangulate a polygon in 3d by projecting vertices into the
  /// best fitting plane and triangulating in 2d.
  /// Return the new added edges.
  /// If the triangulation is not possible (the boundary self intersects in this 2d plane) then
  /// the return vector is empty
  template <typename Polyhedron>
  static bool triangulate_polygon_3d(Polyhedron& P,
                                     const typename Polyhedron::Halfedge_handle h,
                                     bool check_for_intersections = true,
                                     bool refine = true)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    // typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Traits::FT FT;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    typedef CGAL::Aff_transformation_3<typename Polyhedron::Traits> Aff_transformation_3;


    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename InexactKernel::Point_2 InexactPoint_2;
    // typedef typename InexactKernel::Segment_2 InexactSegment_2;
    typedef typename InexactKernel::Segment_3 InexactSegment_3;
    typedef typename InexactKernel::Triangle_2 InexactTriangle_2;

    typedef CGAL::Triangulation_vertex_base_with_info_2<std::pair<Vertex_handle, std::size_t>, InexactKernel> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<InexactKernel> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
    typedef CGAL::No_intersection_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<InexactKernel, TDS, Itag> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Mesh_criteria_2;
    typedef CGAL::Delaunay_mesher_no_edge_refinement_2<CDT, Mesh_criteria_2>  CGAL_Mesher_2;

    std::cout << "Triangulating hole as 2d polygon" << std::endl;
    dolfin_assert(P.is_valid());

    {
      std::cout << "Polygon ";
      typename Polyhedron::Halfedge_handle current = h;
      do
      {
        std::cout << current->vertex()->point() << ", ";
        current = current->next();
      } while (current != h);

       std::cout << std::endl;
    }

    Aff_transformation_3 to_xy;
    Vector_3 plane_normal;
    {
      // Compute the best fitting plane of the points of the hole
      InexactPlane_3 fitting_plane;


      std::vector<InexactSegment_3> boundary;
      Halfedge_handle current = h;
      do
      {
        const Point_3& p = current->vertex()->point();
        const Point_3& next = current->next()->vertex()->point();
        boundary.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2])),
                                            InexactPoint_3(CGAL::to_double(next[0]), CGAL::to_double(next[1]), CGAL::to_double(next[2]))));

        current = current->next();
      } while (current != h);

      /* const double fit_quality = */
      CGAL::linear_least_squares_fitting_3(boundary.begin(),
                                           boundary.end(),
                                           fitting_plane,
                                           CGAL::Dimension_tag<1>());

      // Compute rotation that will rotate the fitting plane to the xy plane
      const InexactVector_3 orthogonal = fitting_plane.orthogonal_vector();
      const Vector_3 orthogonal_exact(orthogonal[0], orthogonal[1], orthogonal[2]);
      // std::cout << "Orthogonal vector: " << orthogonal << ", " << orthogonal.squared_length() << std::endl;
      const Aff_transformation_3 rotation = rotate_to_xy(orthogonal_exact);
      // std::cout << "Rotation: " << rotation << std::endl;
      std::cout << "Rotated: " << rotation.transform(orthogonal_exact) << std::endl;
      const InexactPoint_3 point_on_plane = fitting_plane.point();
      const FT z = rotation.transform(Point_3(point_on_plane[0],
                                              point_on_plane[1],
                                              point_on_plane[2]))[2];

      const Aff_transformation_3 zero_z(CGAL::Translation(), Vector_3(0,0,-z));

      double max_cos_normal_angle = 0;
      const InexactVector_3 normal = fitting_plane.orthogonal_vector();
      std::cout << "Normal: " << plane_normal << std::endl;
      plane_normal = Vector_3(normal.x(), normal.y(), normal.z());

      double max_squared_distance = 0;
      dolfin_assert(dolfin::near(normal.squared_length(), 1, DOLFIN_EPS_LARGE));

      // Compute 2D bounding box
      std::vector<InexactPoint_2> points_2D;
      for (const InexactSegment_3& s : boundary)
      {
        const Point_3 current_exact(s.source()[0], s.source()[1], s.source()[2]);
        const Point_3 rotated = rotation.transform(current_exact);
        points_2D.push_back(InexactPoint_2(CGAL::to_double(rotated[0]),
                                           CGAL::to_double(rotated[1])));
      }
      CGAL::Bbox_2 bbox = CGAL::bbox_2(points_2D.begin(), points_2D.end());
      std::cout << "Bounding box: " << bbox << std::endl;

      // InexactVector_3 prev = InexactVector_3(boundary[boundary.size()-1]) / std::sqrt(boundary[boundary.size()-1].squared_length());
      for (const InexactSegment_3& s : boundary)
      {
        InexactVector_3 current = InexactVector_3(s) / std::sqrt(s.squared_length());
        max_squared_distance = std::max(max_squared_distance, (s.source()-fitting_plane.projection(s.source())).squared_length());

        // std::cout << "Length: " << sit->squared_length() << ", " << current.squared_length() << ", " << (current*normal) << std::endl;
        dolfin_assert(dolfin::near(current.squared_length(), 1, DOLFIN_EPS_LARGE));
        max_cos_normal_angle = std::max(max_cos_normal_angle, CGAL::abs(current*normal));
        // prev = current;

        const Point_3 current_exact(s.source()[0], s.source()[1], s.source()[2]);
        const Point_3 rotated = rotation.transform(current_exact);
        // bbox += InexactPoint_2(CGAL::to_double(rotated[0]), CGAL::to_double(rotated[1])).bbox();
        // std::cout << "  " << bbox << std::endl;
      }

      std::cout << "Max abs cos normal angle: " << max_cos_normal_angle << std::endl;
      // std::cout << "Plane quality: " << fit_quality << std::endl;
      std::cout << "Max distance: " << std::sqrt(max_squared_distance) << std::endl;

      if (max_cos_normal_angle > .2)
      {
        std::cout << "ERROR: Rejecting 2d triangulating, max_cos_normal_angle: " << max_cos_normal_angle << std::endl;
        return false;
      }
      if (max_squared_distance > 1e-8)
      {
        std::cout << "ERROR: Rejecting 2d triangulating, max squared_distance: " << max_squared_distance << std::endl;
        return false;
      }

      const Aff_transformation_3 normalization(CGAL::Translation(), Vector_3(-(bbox.xmax()+bbox.xmin())/2,
                                                                             -(bbox.ymax()+bbox.ymin())/2,
                                                                             0));

      to_xy = normalization*zero_z*rotation;
    }


    // std::cout << "Rotate normal: " << rotation.transform(fitting_plane.orthogonal_vector()) << std::endl;
    // dolfin_assert(dolfin::near(fitting_plane.orthogonal_vector().squared_length(), 1, DOLFIN_EPS_LARGE));
    double max_z = 0.;

    CDT cdt;

    std::cout << "Projected polygon" << std::endl;
    std::cout << "Polygon";

    // Insert vertices into 2D triangulation
    std::vector<typename CDT::Vertex_handle> vertices;
    double max_squared_edge_length = 0;
    double min_squared_edge_length = std::numeric_limits<double>::max();

    {
      Halfedge_handle current = h;
      Point_3 prev = current->prev()->vertex()->point();
      std::stringstream ss;
      ss << "Polygon ";
      do
      {

        const Point_3& p = current->vertex()->point();
        // InexactPoint_3 p_inexact(CGAL::to_double(p[0]),
//                                  CGAL::to_double(p[1]),
        //                               CGAL::to_double(p[2]));
        // InexactPoint_3 p_projected = fitting_plane.projection(p_inexact);
        // std::cout << " " << p_projected << ", ";

        const double length_current = CGAL::to_double(Segment_3(prev, p).squared_length());
        max_squared_edge_length = std::max(max_squared_edge_length, length_current);
        min_squared_edge_length = std::min(min_squared_edge_length, length_current);
        const Point_3 rotated = to_xy.transform(p);

        max_z = std::max(max_z, CGAL::to_double(CGAL::abs(rotated.z())));
        ss << " " << rotated << ", ";

        const InexactPoint_2 p_2d(CGAL::to_double(rotated[0]), CGAL::to_double(rotated[1]));

        // std::cout << " " << p_2d << ", ";

        vertices.push_back(cdt.insert(p_2d));
        vertices.back()->info().first = current->vertex();

        prev = p;
        current = current->next();
      } while (current != h);

      // std::cout << std::endl;
      std::cout << std::endl;
      std::cout << ss.str() << std::endl;
    }

    std::cout << "Size of points: " << vertices.size() << std::endl;
    // std::cout << "z = " << z << std::endl;
    std::cout << "Max z : " << max_z << std::endl;
    std::cout << "Longest edge:  " << max_squared_edge_length << std::endl;
    std::cout << "Shortest edge: " << min_squared_edge_length << std::endl;


    // Check if any of the edges intersect (before actually adding the
    // constrained edges to the triangulation)
#if 0
    if (check_for_intersections)
    {
      for (std::size_t i = 0; i < vertices.size()-1; i++)
      {
        const Point_3& a = vertices[i]->info().first->point(), b = vertices[+1]->info().first->point();
        const InexactSegment_2 s(vertices[i]->point(), vertices[i+1]->point());
        const Segment_3 original(a, b);

        const InexactSegment_3 s2(fitting_plane.projection(InexactPoint_3(CGAL::to_double(a[0]),
                                                                          CGAL::to_double(a[1]),
                                                                          CGAL::to_double(a[2]))),
                                  fitting_plane.projection(InexactPoint_3(CGAL::to_double(b[0]),
                                                                          CGAL::to_double(b[1]),
                                                                          CGAL::to_double(b[2]))));

        for (std::size_t j = i+1; j < vertices.size(); j++)
        {
          InexactSegment_2 s2(vertices[j]->point(), vertices[(j+1)%vertices.size()]->point());

          const auto intersection = CGAL::intersection(s, s2);

          if (intersection)
          {
            if (boost::get<InexactPoint_2>(&*intersection))
            {
              if (j != i+1 && i != (j+1)%vertices.size())
              {
                std::cout << "Non-neighbors (" << i << ", " << j << ")/" << vertices.size()
                          << " intersect in single point" << std::endl;

                return false;
              }
            }
            else if (boost::get<InexactSegment_2>(&*intersection))
            {
              std::cout << "Intersects in segment" << std::endl;
              return false;
            }
            else
            {
              dolfin_assert(false);
              return false;
            }
          } // end if intersection
        } // end inner loop
      } // end outer loop

      // No edges intersect, so we can safely insert then as constraints to the
      // triangulation
    }
#endif

    // Insert the edges around the facet as constraints to the triangulation
    for (std::size_t i = 0; i < vertices.size(); i++)
    {
      // std::cout << "Insert constraint: (" << i << ") " << vertices[i]->point() << ", (" << (i+1)%vertices.size() << ") " << vertices[(i+1)%vertices.size()]->point() << std::endl;
      cdt.insert_constraint(vertices[i], vertices[(i+1)%vertices.size()]);
    }

    // std::cout << "Done triangulating" << std::endl;
    // std::cout << "Num vertices: " << cdt.number_of_vertices() << std::endl;

    if (refine)
    {
      // Create mesher
      CGAL_Mesher_2 mesher(cdt);

      // Set shape and size criteria
      mesher.set_criteria(Mesh_criteria_2(.125, 2*std::sqrt(max_squared_edge_length)));

      // No size criteria, only shape
      //mesher.set_criteria(Mesh_criteria_2());

      std::cout << "Max edge length: " << 2*std::sqrt(max_squared_edge_length) << std::endl;

      // Refine CGAL mesh/triangulation
      std::cout << "Refining 2D mesh" << std::endl;
      mesher.refine_mesh();
    }
    else
    {
      for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
      {
        f->set_in_domain(true);
      }
    }

    std::cout << "Done meshing. Num vertices: " << cdt.number_of_vertices() << std::endl;

    {
      double shortest_edge = std::numeric_limits<double>::max();
      double smallest_triangle = std::numeric_limits<double>::max();
      for (auto f = cdt.finite_faces_begin(); f != cdt.finite_faces_end(); f++)
      {
        if (f->is_in_domain())
        {
          for (int i = 0; i < 3; i++)
            shortest_edge = std::min(shortest_edge,
                                     CGAL::to_double((f->vertex(i)->point()-f->vertex((i+1)%3)->point()).squared_length()));
          InexactTriangle_2 t(f->vertex(0)->point(),
                              f->vertex(1)->point(),
                              f->vertex(2)->point());
          smallest_triangle = std::min(smallest_triangle, CGAL::to_double(CGAL::abs(t.area())));
        }
      }

      std::cout << "Shortest edge in 2D triangulation: " << shortest_edge << std::endl;
      std::cout << "Smallest triangle in 2D triangulation: " << smallest_triangle << std::endl;
    }

    //P.normalize_border();
    // dolfin_assert(P.is_valid(false, 1));

    // add the triangulation to the polyhedron

    dump_2D_triangulation(cdt, "triangulation2D.off");

    Add2DTriangulation<typename Polyhedron::HalfedgeDS, CDT> builder(cdt,
                                                                     to_xy.inverse(),
                                                                     vertices[0],
                                                                     vertices[1],
                                                                     plane_normal*.001);
    P.delegate(builder);

    std::cout << "Inserted: " << builder.new_facets.size() << " new facets" << std::endl;

    // Find the border edge incident to "inserted h"
    const typename Polyhedron::Halfedge_around_vertex_circulator start = builder.v0->vertex_begin();
    typename Polyhedron::Halfedge_around_vertex_circulator current = start;

    do
    {
      if (current->opposite()->vertex() == builder.v1)
      {
        std::cout << "Found the edge" << std::endl;
        break;
      }
      std::cout << "Advancing" << std::endl;
      current++;
    } while (current != start);

    const typename Polyhedron::Halfedge_handle h_new = current->is_border() ? current : current->opposite();

    if (!h_new->is_border())
    {
      std::cout << "Is border edge: " << h_new->is_border_edge() << std::endl;
      dolfin::dolfin_error("Polyhedron_utils.h",
                           "Locating border edge",
                           "locating border edge");
    }

    std::cout << "2D: " << vertices[0]->point() << std::endl;
    std::cout << "Original: " << h->vertex()->point() << std::endl;
    std::cout << "Original Next: " << h->next()->vertex()->point() << std::endl;
    std::cout << "Original prev: " << h->prev()->vertex()->point() << std::endl;
    std::cout << "Inserted: " << h_new->vertex()->point() << std::endl;
    std::cout << "Inserted next: " << h_new->next()->vertex()->point() << std::endl;
    std::cout << "Inserted prev: " << h_new->prev()->vertex()->point() << std::endl;
    std::cout << "Vertex degree: " << h_new->vertex()->vertex_degree() << std::endl;


    std::cout << "Joining loop" << std::endl;

    {
      Halfedge_handle a = h_new;
      Halfedge_handle b = h;
      double max_distance = 0.;
      do
      {
        max_distance = std::max(max_distance, CGAL::to_double((a->vertex()->point()-b->vertex()->point()).squared_length()));

        a = a->next();
        b = b->prev();
      } while (a != h_new);
      std::cout << "Max distance between merged vertices: " << max_distance << std::endl;
    }

    {
      double shortest_edge = std::numeric_limits<double>::max();
      for (auto e = P.halfedges_begin(); e != P.halfedges_end(); e++)
      {
        shortest_edge = std::min(shortest_edge, CGAL::to_double((e->vertex()->point()-e->opposite()->vertex()->point()).squared_length()));
      }

      std::cout << "Shortest edge: " << shortest_edge << std::endl;
    }

    {
      double smallest_triangle = std::numeric_limits<double>::max();;
      for (auto f = P.facets_begin(); f != P.facets_end(); f++)
      {
        auto h = f->halfedge();
        if (f->is_triangle())
        {
          typename Polyhedron::Traits::Triangle_3 t(h->vertex()->point(),
                                                    h->next()->vertex()->point(),
                                                    h->next()->next()->vertex()->point());
          smallest_triangle = std::min(smallest_triangle,
                                      CGAL::to_double(t.squared_area()));
        }
      }
      std::cout << "Smallest triangle before join: " << smallest_triangle << std::endl;
    }


    //
    P.join_loop(h, h_new);

    {
      double shortest_edge = std::numeric_limits<double>::max();
      for (auto e = P.halfedges_begin(); e != P.halfedges_end(); e++)
      {
        shortest_edge = std::min(shortest_edge, CGAL::to_double((e->vertex()->point()-e->opposite()->vertex()->point()).squared_length()));
      }

      std::cout << "Shortest edge after merge: " << shortest_edge << std::endl;
    }

    {
      double smallest_triangle = std::numeric_limits<double>::max();;
      for (auto f = P.facets_begin(); f != P.facets_end(); f++)
      {
        auto h = f->halfedge();
        if (f->is_triangle())
        {
          typename Polyhedron::Traits::Triangle_3 t(h->vertex()->point(),
                                                    h->next()->vertex()->point(),
                                                    h->next()->next()->vertex()->point());
          smallest_triangle = std::min(smallest_triangle,
                                      CGAL::to_double(t.squared_area()));
        }
      }
      std::cout << "Smallest triangle: " << smallest_triangle << std::endl;
    }



    std::cout << "H degree: " << h->facet()->facet_degree() << std::endl;
    std::cout << "Hole: " << (h->is_border_edge() ? "Yes" : "No") << std::endl;

    return true;
  }
  //-----------------------------------------------------------------------------
  template <typename Polyhedron>
  static std::pair<double, double> evaluate_planarity(const typename Polyhedron::Halfedge_handle from,
                                                      const typename Polyhedron::Halfedge_handle to)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Segment_3 InexactSegment_3;
    typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    std::vector<InexactSegment_3> segments;
    Halfedge_handle current = from;
    while (current != to)
    {
      const Point_3& s1 = current->vertex()->point();
      const Point_3& s2 = current->next()->vertex()->point();
      segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(s1[0]),
                                                         CGAL::to_double(s1[1]),
                                                         CGAL::to_double(s1[2])),
                                          InexactPoint_3(CGAL::to_double(s2[0]),
                                                         CGAL::to_double(s2[1]),
                                                         CGAL::to_double(s2[2]))));
      current = current->next();
    }

    InexactPlane_3 fitting_plane;
    CGAL::linear_least_squares_fitting_3(segments.begin(),
                                         segments.end(),
                                         fitting_plane,
                                         CGAL::Dimension_tag<1>());

    const InexactVector_3 orthogonal_vector = fitting_plane.orthogonal_vector()/sqrt(CGAL::to_double(fitting_plane.orthogonal_vector().squared_length()));

    dolfin_assert(dolfin::near(orthogonal_vector.squared_length(), 1., DOLFIN_EPS_LARGE));

    double max_distance = 0.;
    double min_cos_angle = 1.;

    current = from;;
    while (current != to)
    {
      const Point_3& p1 = current->vertex()->point();
      const InexactPoint_3 p1_inexact(CGAL::to_double(p1[0]),
                                      CGAL::to_double(p1[1]),
                                      CGAL::to_double(p1[2]));

      const Point_3& next = current->next()->vertex()->point();
      const InexactPoint_3 next_inexact(CGAL::to_double(next[0]),
                                        CGAL::to_double(next[1]),
                                        CGAL::to_double(next[2]));

      const InexactVector_3 current_vector = (p1_inexact-next_inexact)/sqrt(CGAL::to_double((p1_inexact-next_inexact).squared_length()));
      dolfin_assert(dolfin::near(current_vector.squared_length(), 1., DOLFIN_EPS_LARGE));

      // Check if the angle between this vector and the fitting plane is acceptable
      min_cos_angle = std::min(min_cos_angle,
                               CGAL::abs(current_vector*orthogonal_vector));

      max_distance = std::max(max_distance,
                              sqrt((p1_inexact-fitting_plane.projection(p1_inexact)).squared_length()));

      current = current->next();
    }

    return std::make_pair(max_distance, min_cos_angle);
  }
  //-----------------------------------------------------------------------------
  template <typename Polyhedron>
  static std::pair<typename Polyhedron::Halfedge_handle, typename Polyhedron::Halfedge_handle>
  find_best_cut (Polyhedron& P,
                 const typename Polyhedron::Halfedge_handle h,
                 double cos_angle_tolerance)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel;
    typedef typename InexactKernel::Plane_3 InexactPlane_3;
    typedef typename InexactKernel::Point_3 InexactPoint_3;
    typedef typename InexactKernel::Vector_3 InexactVector_3;
    typedef typename InexactKernel::Segment_3 InexactSegment_3;

    Halfedge_handle h1_distance, h2_distance;
    double min_distance = std::numeric_limits<double>::max();

    Halfedge_handle h1_angle, h2_angle;
    double min_cos_normal = std::numeric_limits<double>::max();

    Halfedge_handle current1 = h->next()->next();
    const Halfedge_handle current1_end = h->prev()->prev();
    do
    {
      // std::cout << "Current1 " << current1->vertex()->point() << std::endl;
      Halfedge_handle current2 = current1->next()->next();
      const Halfedge_handle current2_end = current1->prev();
      do
      {
        /* std::cout << "  Current 2" << std::endl; */
        /* std::cout << "Segment " << current1->vertex()->point() << ", " << current2->vertex()->point() << std::endl; */

        double max_distance_local = 0.;
        double max_cos_normal_local = 0;

        // Compute fitting plane of segments from current1 to current2
        // std::cout << "Compute fit plane 1" << std::endl;
        {
          InexactPlane_3 fitting_plane;
          std::vector<InexactSegment_3> segments;
          Halfedge_handle current3 = current1;
          while (current3 != current2)
          {
            const Point_3& s1 = current3->vertex()->point();
            const Point_3& s2 = current3->next()->vertex()->point();
            segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(s1[0]),
                                                               CGAL::to_double(s1[1]),
                                                               CGAL::to_double(s1[2])),
                                                InexactPoint_3(CGAL::to_double(s2[0]),
                                                               CGAL::to_double(s2[1]),
                                                               CGAL::to_double(s2[2]))));
            current3 = current3->next();
          }

          /* const double fit_quality = */
          CGAL::linear_least_squares_fitting_3(segments.begin(),
                                               segments.end(),
                                               fitting_plane,
                                               CGAL::Dimension_tag<1>());
          const InexactVector_3 orthogonal_vector = fitting_plane.orthogonal_vector()/sqrt(CGAL::to_double(fitting_plane.orthogonal_vector().squared_length()));

          // std::cout << "Plane vector length: " << (CGAL::abs(orthogonal_vector.squared_length())-1) << std::endl;
          dolfin_assert(dolfin::near(orthogonal_vector.squared_length(), 1., DOLFIN_EPS_LARGE));

          current3 = current1;
          while (current3 != current2)
          {
            const Point_3& p1 = current3->vertex()->point();
            const InexactPoint_3 p1_inexact(CGAL::to_double(p1[0]),
                                            CGAL::to_double(p1[1]),
                                            CGAL::to_double(p1[2]));

            const Point_3& next = current3->next()->vertex()->point();
            const InexactPoint_3 next_inexact(CGAL::to_double(next[0]),
                                              CGAL::to_double(next[1]),
                                              CGAL::to_double(next[2]));

            const InexactVector_3 current_vector = (p1_inexact-next_inexact)/sqrt(CGAL::to_double((p1_inexact-next_inexact).squared_length()));
            dolfin_assert(dolfin::near(current_vector.squared_length(), 1., DOLFIN_EPS_LARGE));

            // Check if angle between this vector and the fitting plane is acceptable
            max_cos_normal_local = std::max(max_cos_normal_local, CGAL::abs(current_vector*orthogonal_vector));

            max_distance_local = std::max(max_distance_local,
                                          CGAL::to_double((p1_inexact-fitting_plane.projection(p1_inexact)).squared_length()));

            current3 = current3->next();
          }
        }

        // std::cout << "Done evaluating candidate: " << std::endl;

        // Compute fitting plane of segments from current2 to current1
        {
          InexactPlane_3 fitting_plane;
          std::vector<InexactSegment_3> segments;
          Halfedge_handle current3 = current2;
          while (current3 != current1)
          {
            const Point_3& s1 = current3->vertex()->point();
            const Point_3& s2 = current3->next()->vertex()->point();
            segments.push_back(InexactSegment_3(InexactPoint_3(CGAL::to_double(s1[0]),
                                                               CGAL::to_double(s1[1]),
                                                               CGAL::to_double(s1[2])),
                                                InexactPoint_3(CGAL::to_double(s2[0]),
                                                               CGAL::to_double(s2[1]),
                                                               CGAL::to_double(s2[2]))));
            current3 = current3->next();
          }

          /* const double fit_quality = */ CGAL::linear_least_squares_fitting_3(segments.begin(),
                                                                          segments.end(),
                                                                          fitting_plane,
                                                                          CGAL::Dimension_tag<1>());

          const InexactVector_3 orthogonal_vector = fitting_plane.orthogonal_vector()/sqrt(CGAL::to_double(fitting_plane.orthogonal_vector().squared_length()));
          // const double angle_tolerance = .5;
          // std::cout << "Plane vector length: " << (CGAL::abs(orthogonal_vector.squared_length())-1) << std::endl;
          dolfin_assert(dolfin::near(CGAL::to_double(orthogonal_vector.squared_length()), 1., DOLFIN_EPS_LARGE));

          current3 = current2;
          while (current3 != current1)
          {
            const Point_3& p1 = current3->vertex()->point();
            const InexactPoint_3 p1_inexact(CGAL::to_double(p1[0]),
                                            CGAL::to_double(p1[1]),
                                            CGAL::to_double(p1[2]));

            const Point_3& next = current3->next()->vertex()->point();
            const InexactPoint_3 next_inexact(CGAL::to_double(next[0]),
                                              CGAL::to_double(next[1]),
                                              CGAL::to_double(next[2]));

            const InexactVector_3 current_vector = (p1_inexact-next_inexact)/sqrt(CGAL::to_double((p1_inexact-next_inexact).squared_length()));
            dolfin_assert(dolfin::near(CGAL::to_double(current_vector.squared_length()), 1., DOLFIN_EPS_LARGE));

            // Check if the angle between this vector and the fitting plane is acceptable
            max_cos_normal_local = std::max(max_cos_normal_local, CGAL::abs(current_vector*orthogonal_vector));

            max_distance_local = std::max(max_distance_local,
                                          CGAL::to_double((p1_inexact-fitting_plane.projection(p1_inexact)).squared_length()));

            current3 = current3->next();
          }
        }

        // std::cout << "Done evaluating other half: " << std::endl;

        if (max_distance_local < min_distance)
        {
          // std::cout << "New min split" << std::endl;
          min_distance = max_distance_local;
          h1_distance = current1;
          h2_distance = current2;
        }

        if (max_cos_normal_local < min_cos_normal)
        {
          min_cos_normal = max_cos_normal_local;
          h1_angle = current1;
          h2_angle = current2;
        }

        current2 = current2->next();
      } while (current2 != current2_end);

      current1 = current1->next();
    } while (current1 != current1_end);

    std::cout << "Best cut, distance=" << min_distance << " : Segment " << h1_distance->vertex()->point() << ", " << h2_distance->vertex()->point() << std::endl;
    std::cout << "Best cut, angle=" << min_cos_normal << " : Segment " << h1_angle->vertex()->point() << ", " << h2_angle->vertex()->point() << std::endl;


    return std::make_pair(h1_angle, h2_angle);
  }
  //-----------------------------------------------------------------------------
  template <typename Polyhedron>
  static bool split_hole_planar(Polyhedron& P,
                                typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;

    std::cout << "Split hole" << std::endl;

    const std::pair<Halfedge_handle, Halfedge_handle> best_cut = find_best_cut(P, h, .3);

    std::cout << "Found best cut" << std::endl;
    if (best_cut.first == Halfedge_handle() || best_cut.second == Halfedge_handle())
      return false;


    const std::pair<double, double> planarity1 = evaluate_planarity<Polyhedron>(best_cut.first, best_cut.second);
    const std::pair<double, double> planarity2 = evaluate_planarity<Polyhedron>(best_cut.second, best_cut.first);

    std::cout << "Planarity: " << planarity1.first << " vs " << planarity2.first << std::endl;

    const Halfedge_handle start = planarity1.first < planarity2.first ? best_cut.first : best_cut.second;
    const Halfedge_handle end   = planarity1.first < planarity2.first ? best_cut.second : best_cut.first;

    // Compute average edge length
    double length = 0.;
    int num_edges = 0;
    Halfedge_handle current = best_cut.first;
    do
    {
      length += sqrt(CGAL::to_double((current->vertex()->point()-current->next()->vertex()->point()).squared_length()));
      num_edges++;
    } while (current != best_cut.first);

    length /= num_edges;

    const Vector_3 cut = start->vertex()->point() - end->vertex()->point();
    const double cut_length = sqrt(CGAL::to_double(cut.squared_length()));
    const int num_new_edges = static_cast<int>(cut_length/length);

    P.fill_hole(h);

    Halfedge_handle new_diagonal = P.split_facet(end, start);
    // Now new_diagonal is on the facet to be

    std::cout << "Inserting " << num_new_edges << " vertices from " << start->vertex()->point() << " to " << end->vertex()->point() << std::endl;
    for (int i = 1; i < num_new_edges; i++)
    {
      Halfedge_handle hnew = P.split_edge(new_diagonal);
      hnew->vertex()->point() = end->vertex()->point() + cut*i/num_new_edges;
      // std::cout << "  " << i << ": " << hnew->vertex()->point() << (hnew->vertex()->point()-hnew->opposite()->vertex()->point()).squared_length() << std::endl;
    }

    /* Halfedge_handle c = new_diagonal; */
    /* do */
    /* { */
    /*   std::cout << "Edge: " << c->vertex()->point() << " <--> " << c->opposite()->vertex()->point() << " : " << (c->vertex()->point()-c->opposite()->vertex()->point()).squared_length() << std::endl; */
    /*   c = c->next(); */
    /* } while (c != new_diagonal); */


    const Halfedge_handle opposite = new_diagonal->opposite();
    P.make_hole(new_diagonal);
    const bool success = triangulate_polygon_3d(P, new_diagonal);
    //dolfin_assert(success);

    P.make_hole(opposite);
    if (!success)
      return false;

    return true;
  }
  //-----------------------------------------------------------------------------

  /* template <typename Polyhedron> */
  /* void min_vertex_degree(const Polyhedron& p) */
  /* { */
  /*   std::size_t min_degree = std::numeric_limits<std::size_t>::max(); */
  /*   std::size_t min_degree_non_border = min_degree; */

  /*   for (typename Polyhedron::Vertex_const_iterator vit = p.vertices_begin(); */
  /*        vit != p.vertices_end(); vit++) */
  /*   { */
  /*     min_degree = std::min(min_degree, vit->vertex_degree); */

  /*   } */

  /*   std::cout << "Min vertex_degree: " << min_degree << std::endl; */
  /* } */
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void list_hole(const typename Polyhedron::Halfedge_handle h)
  {
    std::size_t counter = 0;
    std::cout << "Polygon";

    {
      typename Polyhedron::Halfedge_handle current = h;
      do
      {
        counter++;
        current = current->next();
      } while(current != h);
    }

    // if (counter < 250)
    // {
      typename Polyhedron::Halfedge_handle current = h;
      do
      {
        std::cout << " " << current->vertex()->point() << ",";

        current = current->next();
      } while(current != h);
      // }
      std::cout << std::endl;

      // std::cout << " size: " << counter << std::endl;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::string print_triangle(typename Polyhedron::Halfedge_handle h)
  {
    std::stringstream ss;
    ss << "Triangle "
              << h->prev()->vertex()->point() << ", "
              << h->vertex()->point() << ", "
              << h->next()->vertex()->point() << std::endl;
    ss << "Area: " << triangle_area<Polyhedron>(h) << std::endl;
    return ss.str();
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static double triangle_area(typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    Triangle_3 t(h->prev()->vertex()->point(),
                 h->vertex()->point(),
                 h->next()->vertex()->point());

    return CGAL::to_double(t.squared_area());
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Vertex_handle get_common_vertex(typename Polyhedron::Facet_handle f1,
                                                              typename Polyhedron::Facet_handle f2)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;

    // Find common vertex
    Halfedge_handle h1 = f1->halfedge();
    Halfedge_handle current1 = h1;
    do
    {
      Halfedge_handle h2 = f2->halfedge();
      Halfedge_handle current2 = h2;
      do
      {
        if (current2->vertex() == current1->vertex())
          return current2->vertex();

        current2 = current2->next();
      } while (h2 != current2);

      current1 = current1->next();
    } while (h1 != current1);

    return Vertex_handle();
  }
  //-----------------------------------------------------------------------------
  template <typename Polyhedron>
  static double facet_angle(typename Polyhedron::Halfedge_handle h)
  {
    dolfin_assert(h->is_border());
    dolfin_assert(h->next()->is_border());

    typedef typename Polyhedron::Traits::Line_3 Line_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    const Line_3 l(h->prev()->vertex()->point(),
                   h->vertex()->point());

    const Point_3& p11 = h->next()->vertex()->point();
    const Point_3  p12 = l.projection(p11);
    const Vector_3 v1(p12, p11);

    dolfin_assert(h->opposite()->facet()->is_triangle());

    const Point_3& p21 = h->opposite()->next()->vertex()->point();
    const Point_3  p22 = l.projection(p21);
    const Vector_3 v2(p22, p21);

    return CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())));
  }
  //-----------------------------------------------------------------------------
  template<typename Segment_3>
  static bool segment_intersects_triangle(const Segment_3& s,
                                          const typename CGAL::Kernel_traits<Segment_3>::Kernel::Triangle_3& t)
  {
    typedef typename CGAL::Kernel_traits<Segment_3>::Kernel::Point_3 Point_3;

    auto result = CGAL::intersection(s, t);
    if (!result)
      return false;

    if (const Point_3* p = boost::get<Point_3>(&*result))
    {
      if (*p == t[0] || *p == t[1] || *p == t[2])
        return false;
      else
        return true;
    }
    else if (const Segment_3* s_ = boost::get<Segment_3>(&*result))
    {
      if ( (s.source() == t[0] || s.source() == t[1] || s.source() == t[2]) &&
           (s.target() == t[0] || s.target() == t[1] || s.target() == t[2]) )
        return false;
      else
        return true;
    }

    dolfin_assert(false);
    return false;
  }
  //-----------------------------------------------------------------------------
  // Check if two triangles intersect.
  // Neighbor triangles (share a vertice or an edge) do not intersect
  // if t1 == t2 (geometrically), the triangles do not intersect
  template<typename Triangle_3>
  static bool triangles_intersect(const Triangle_3& t1, const Triangle_3& t2)
  {
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Point_3 Point_3;
    typedef typename CGAL::Kernel_traits<Triangle_3>::Kernel::Segment_3 Segment_3;

    if ( (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2]) &&
         (t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2]) &&
         (t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2]))
      return false;

    auto result = CGAL::intersection(t1, t2);

    if (!result)
      return false;

    if (const Point_3* p = boost::get<Point_3>(&*result))
    {
      if (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2] ||
          t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2] ||
          t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2])
        return false;
      else
        return true;
    }
    else if (const Segment_3* s = boost::get<Segment_3>(&*result))
    {
      std::size_t common_vertices = 0;
      if (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2])
        common_vertices++;

      if (t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2])
        common_vertices++;

      if (t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2])
        common_vertices++;

      if (common_vertices > 1)
        return false;
      else
        return true;
    }
    else if (const Triangle_3* t = boost::get<Triangle_3>(&*result))
      return true;
    else if (const std::vector<Point_3>* v = boost::get<std::vector<Point_3> >(&*result))
      return true;

    dolfin_assert(false);
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Triangle_3>
  static bool triangle_set_intersects(const std::vector<Triangle_3>& t)
  {
    for (std::size_t i = 0; i < t.size(); i++)
    {
      for (std::size_t j = i+1; j < t.size(); j++)
      {
        if (triangles_intersect<Triangle_3>(t[i], t[j]))
          return true;
      }
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Halfedge_handle>
  static std::size_t vertex_count_halfedges(Halfedge_handle h)
  {
    std::size_t count = 0;
    Halfedge_handle current = h;
    do
    {
      ++count;
      current = current->next()->opposite();
    } while (current != h);

    dolfin_assert(count > 0);
    return count;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static std::size_t total_vertex_count_halfedges(const Polyhedron& P, std::map<typename Polyhedron::Vertex_const_handle, std::size_t>& m)
  {
    std::size_t total_count = 0;
    for (typename Polyhedron::Vertex_const_iterator it = P.vertices_begin(); it != P.vertices_end(); it++)
    {
      const std::size_t c = vertex_count_halfedges(it->halfedge());
      m[it] = c;
      // std::cout << " Halfedge count: " << c << std::endl;
      if (c == 0)
      {int tmp; std::cin >> tmp; }
      total_count += c;
    }
    return total_count;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool halfedge_is_in_polyhedron(const Polyhedron& P,
                                        typename Polyhedron::Halfedge_const_handle h)
  {
    for (auto hit = P.halfedges_begin(); hit != P.halfedges_end(); hit++)
    {
      if (hit == h)
        return true;
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool vertex_is_in_polyhedron(const Polyhedron& P,
                                      typename Polyhedron::Vertex_const_handle v)
  {
    for (auto vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    {
      if (vit == v)
        return true;
    }

    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
    static bool check_vertex_consistency(const Polyhedron& P)
  {
    // std::cout << "Checking vertex consistency" << std::endl;
    std::size_t counter = 0;

    // Build a set of the list of vertices for faster lookup
    std::set<typename Polyhedron::Vertex_const_handle> vertex_set;
    for (auto vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    {
      dolfin_assert(vertex_set.count(vit) == 0);
      vertex_set.insert(vit);
    }

    std::deque<typename Polyhedron::Halfedge_const_handle> queue;
    std::set<typename Polyhedron::Halfedge_const_handle> visited;

    queue.push_back(P.halfedges_begin());
    while (!queue.empty())
    {
      typename Polyhedron::Halfedge_const_handle current = queue.back();
      queue.pop_back();
      if (visited.count(current) == 0)
      {
        counter++;
        typename Polyhedron::Halfedge_const_handle start = current;
        visited.insert(current);

        // Walk around the facet (or hole). Check halfedges and queue opposites
        do
        {
          if (vertex_set.count(current->vertex()) == 0)
          {
            // std::cout << "Vertex not in vertex list: " << current->vertex()->point() << std::endl;
            return false;
          }

          // TODO: Add the opposite check: All vertices should be reachable via halfedges
          queue.push_back(current->opposite());
          current = current->next();
        } while(current != start);
      }
    }

    // std::cout << "  Checked " << counter << " halfedges" << std::endl;
    return true;
  }
  //-----------------------------------------------------------------------------
  template<typename Triangle_3>
  static bool triangle_set_intersect_triangle(Triangle_3 t,
                                              const std::vector<Triangle_3>& triangle_set)
  {
    for (auto tit = triangle_set.begin(); tit != triangle_set.end(); tit++)
    {
      if (triangles_intersect(*tit, t))
        return true;
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Segment_3>
    static bool segment_intersects_triangle_set(const Segment_3& s,
                                                const std::vector<typename CGAL::Kernel_traits<Segment_3>::Kernel::Triangle_3>& triangle_set)
  {
    for (auto tit = triangle_set.begin(); tit != triangle_set.end(); tit++)
    {
      if (segment_intersects_triangle(s, *tit))
        return true;
    }
    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool center_vertex_triangulation(Polyhedron& P,
                                          const typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    // check if the triangulation with a center vertex intersects any of the
    // neighbor triangles
    // std::cout << "Attempting center vertex triangulation" << std::endl;
    dolfin_assert(P.is_valid(false, 0));
    dolfin_assert(halfedge_is_in_polyhedron(P, h));
    dolfin_assert(vertex_is_in_polyhedron(P, h->vertex()));

    const std::array<double, 3> plane_fit = get_plane_fit<Polyhedron>(h, h->prev());
    if (plane_fit[0] < .85)
    {
      // std::cout << "  Rejected. Not sufficiently planar: " << plane_fit[0] << std::endl;
      return false;
    }

    // Collect set of neighbor triangles and compute centroid
    std::vector<Triangle_3> triangles;
    Point_3 centroid = CGAL::ORIGIN;
    std::size_t counter = 0;
    {
      Halfedge_handle current = h;
      do
      {
        if (!current->opposite()->is_border() && current->opposite()->facet()->facet_degree() == 3)
        {
          triangles.push_back(get_facet_triangle<Polyhedron>(current->opposite()));
        }
        centroid = centroid + (current->vertex()->point()-CGAL::ORIGIN);
        counter++;

        current = current->next();
      } while (current != h);
    }
    // std::cout << "Number of triangles: " << triangles.size() << std::endl;
    centroid = CGAL::ORIGIN + (centroid-CGAL::ORIGIN)/counter;
    // std::cout << "Centroid: " << centroid << std::endl;

    Halfedge_handle current = h;
    do
    {
      Triangle_3 current_triangle(current->vertex()->point(), current->next()->vertex()->point(), centroid);
      for (auto tit = triangles.begin(); tit != triangles.end(); tit++)
      {
        if (triangles_intersect(current_triangle, *tit))
        {
          // std::cout << "No: Triangle " << current_triangle[0] << " " << current_triangle[1] << " " << current_triangle[2] << std::endl;
          // std::cout << "Triangle " << (*tit)[0] << " " << (*tit)[1] << " " << (*tit)[2] << std::endl;
          return false;
        }
      }

      triangles.push_back(current_triangle);

      current = current->next();
    } while (current != h);

    // std::cout << "Facet degree before center vertex: " << h->facet()->facet_degree() << std::endl;

    P.normalize_border();
    dolfin_assert(P.is_valid(false, 1));
    dolfin_assert(!h->is_border_edge());
    dolfin_assert(h->facet()->facet_degree() > 3);

    Halfedge_handle center = P.create_center_vertex(h);
    /* Halfedge_handle g = h->next()->next(); */
    /* dolfin_assert(check_vertex_consistency(P)); */
    /* dolfin_assert(halfedge_is_in_polyhedron(P, h)); */
    /* dolfin_assert(vertex_is_in_polyhedron(P, h->vertex())); */
    /* dolfin_assert(vertex_is_in_polyhedron(P, g->vertex())); */
    /* std::cout << "Splitting: Segment " << h->vertex()->point() << ", " << g->vertex()->point() << std::endl; */
    /* std::cout << "Vertex degrees: " << h->vertex()->vertex_degree() << ", " << g->vertex()->vertex_degree() << std::endl; */
    /* std::cout << "My count: " << vertex_count_halfedges(h) << ", " << vertex_count_halfedges(g) << std::endl; */

    /* std::cout << "Facet degree: " << h->facet()->facet_degree() << std::endl; */
    /* std::cout << "Opposite facet degree: " << h->opposite()->facet()->facet_degree() << std::endl; */
    /* std::cout << "Num halfedges: " << P.size_of_halfedges() << std::endl; */
    /* std::map<typename Polyhedron::Vertex_const_handle, std::size_t> vertex_degrees; */
    /* std::cout << "My total count: " << total_vertex_count_halfedges(P, vertex_degrees); */
    /* std::cout << "Mapped: " << vertex_degrees.at(h->vertex()) << ", " << vertex_degrees.at(g->vertex()) << std::endl; */
    /* std::cout << "Size of map: " << vertex_degrees.size() << std::endl; */
    /* std::cout << "-- Splitting facet --" << std::endl; */

    /* Halfedge_handle diagonal = P.split_facet(h, g); */

    /* std::cout << "Vertex degrees: " << h->vertex()->vertex_degree() << ", " << g->vertex()->vertex_degree() << std::endl; */
    /* std::cout << "My count: " << vertex_count_halfedges(h) << ", " << vertex_count_halfedges(g) << std::endl; */
    /* std::map<typename Polyhedron::Vertex_const_handle, std::size_t> vertex_degrees_after; */
    /* std::cout << "My total count: " << total_vertex_count_halfedges(P, vertex_degrees_after) << std::endl; */
    /* std::cout << "Size of map: " << vertex_degrees_after.size() << std::endl; */
    /* std::cout << "Num halfedges: " << P.size_of_halfedges() << std::endl; */
    /* std::cout << "Mapped: " << vertex_degrees_after.at(h->vertex()) << ", " << vertex_degrees_after.at(g->vertex()) << std::endl; */
    /* for (auto it = vertex_degrees_after.begin(); it != vertex_degrees_after.end(); ++it) */
    /* { */
    /*   if (vertex_degrees.at(it->first) != it->second) */
    /*   { */
    /*     std::cout << "DIFF!!!" << vertex_degrees[it->first] << " " << it->second << std::endl; */
    /*   } */

    /*   if (vertex_degrees.at(it->first) == 0 || it->second == 0) */
    /*   { */
    /*     std::cout << "zero degree " << it->first->point() << std::endl; */
    /*   } */
    /* } */
    /* dolfin_assert(P.is_valid(false)); */
    /* Halfedge_handle c = diagonal->vertex()->halfedge(); */
    /* Halfedge_handle start = c; */
    /* do */
    /* { */
    /*   if (c->opposite()->vertex() == diagonal->opposite()->vertex()) */
    /*     std::cout << "Yes!!!" << std::endl; */
    /*   c = c->opposite()->next(); */
    /* } while (c != start); */
    /* dolfin_assert(P.is_valid(false, 0)); */
    /* std::cout << "Splitting edge: Segment " << diagonal->opposite()->vertex()->point() << ", " << diagonal->vertex()->point() << std::endl; */
    /* Halfedge_handle center = P.split_edge(diagonal); */
    center->vertex()->point() = centroid;
    dolfin_assert(P.is_valid(false, 0));

    /* std::cout << "Splitting facet: " << diagonal->opposite()->vertex()->point() << ", " << diagonal->opposite()->prev()->prev()->vertex()->point() << std::endl; */
    /* P.split_facet(diagonal->opposite(), diagonal->opposite()->prev()->prev()); */

    /* do */
    /* { */
    /*   dolfin_assert(P.is_valid(false, 0)); */
    /*   std::cout << "adding diagonal: " << diagonal->next()->vertex()->point() << ", " << center->vertex()->point() << std::endl; */
    /*   diagonal = P.split_facet(diagonal->next(), center); */
    /*   diagonal = diagonal->opposite(); */

    /* } while (diagonal->next()->next() != center); */
    /* std::cout << "Center vertex degree: " << center->vertex()->vertex_degree() << std::endl; */
    /* //P.normalize_border(); */
    /* { */
    /*   std::ofstream ofile("center-vertex.off"); */
    /*   ofile << P; */
    /* } */

    /* dolfin_assert(P.is_valid()); */

    return true;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void close_hole(Polyhedron& P,
                         typename Polyhedron::Halfedge_handle h)
  {
    //typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    dolfin_assert(P.is_valid(false, 0));
    dolfin_assert(P.is_pure_triangle());
    dolfin_assert(h->is_border());

    P.fill_hole(h);
    P.normalize_border();

    dolfin_assert(h->facet()->facet_degree() > 2);

    // Since the facet may be split, we push the facets to a fifo queue.
    // Neighbor facets are now not guaranteed to be triangles.
    std::deque<Halfedge_handle> queue;
    queue.push_back(h);

    while (!queue.empty())
    {
      // std::cout << "--- Popping facet from queue (" << queue.size() << ")" << std::endl;
      const Halfedge_handle current = queue.front();
      queue.pop_front();

      // list_hole<Polyhedron>(current);

      dolfin_assert(P.is_valid(false, 0));
      dolfin_assert(halfedge_is_in_polyhedron(P, current));

      if (current->facet()->facet_degree() == 3)
      {
        //P.fill_hole(h);
        dolfin_assert(current->opposite()->facet()->facet_degree() != 3  ||
                      !triangles_intersect(get_facet_triangle<Polyhedron>(current),
                                           get_facet_triangle<Polyhedron>(current->opposite())));
        dolfin_assert(current->next()->opposite()->facet()->facet_degree() != 3  ||
                      !triangles_intersect(get_facet_triangle<Polyhedron>(current),
                                           get_facet_triangle<Polyhedron>(current->next()->opposite())));
        dolfin_assert(current->prev()->opposite()->facet()->facet_degree() != 3  ||
                      !triangles_intersect(get_facet_triangle<Polyhedron>(current),
                                           get_facet_triangle<Polyhedron>(current->prev()->opposite())));
      }
      else
      {
        // std::cout << "Attempting to triangulate in 2D" << std::endl;
        dolfin_assert(halfedge_is_in_polyhedron(P, current));
        if (!triangulate_polygon_3d(P, current, false))
        {
          dolfin_assert(halfedge_is_in_polyhedron(P, current));

          Halfedge_handle facet = subdivide_facet(P, current);

          queue.push_back(facet->opposite());
          queue.push_back(facet);
        }
      }
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static typename Polyhedron::Halfedge_handle
  subdivide_facet(Polyhedron& P, typename Polyhedron::Halfedge_handle h)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    typedef typename Polyhedron::Traits::Vector_3 Vector_3;

    // Search for segments that divide the hole, such that the dividing segment
    // does not intersect triangles next to the hole and the facets on each side
    // of the split is as planar as possible.

    // Store all triangles around the hole and compute max edge length
    std::vector<Triangle_3> border_triangles;
    double max_squared_edge_length = 0;
    {
      Halfedge_handle current = h;
      do
      {
        max_squared_edge_length = std::max(max_squared_edge_length,
                                           CGAL::to_double(Segment_3(current->prev()->vertex()->point(),
                                                                     current->vertex()->point()).squared_length()));
        if (current->opposite()->facet()->facet_degree() == 3)
        {
          border_triangles.push_back(get_facet_triangle<Polyhedron>(current->opposite()));
        }

        current = current->next();
      } while (current != h);
    }

    dolfin_assert(border_triangles.size() > 4);

    // Search for the best dividing segment
    double best_quality = -1000;
    Halfedge_handle best_outer;
    Halfedge_handle best_inner;

    Halfedge_handle current_outer = h;

    // FIXME: This loop should run to h->prev()->prev(), but need
    // handle the specially
    const Halfedge_handle outer_end = h->prev()->prev();
    do
    {
      Halfedge_handle current_inner = current_outer->next()->next();
      const Halfedge_handle inner_end = h;
      do
      {
        if (current_inner->next() != current_outer &&
            current_inner->prev() != current_outer)
        {

          Segment_3 current_segment(current_outer->vertex()->point(),
                                    current_inner->vertex()->point());

          // Check that this does not introduce an intersection
          if (!segment_intersects_triangle_set(current_segment, border_triangles))
          {
            const double candidate_quality = evaluate_hole_subdivision<Polyhedron>(current_inner, current_outer);

            if (candidate_quality > best_quality)
            {
              best_outer = current_outer;
              best_inner = current_inner;
              best_quality = candidate_quality;
            }
          }
        }

        current_inner = current_inner->next();
      } while (current_inner != inner_end);

      current_outer = current_outer->next();
    } while (current_outer != outer_end);

    dolfin_assert(best_outer != Halfedge_handle());
    dolfin_assert(best_inner != Halfedge_handle());

    // std::cout << "Found best subdivision: " << std::endl;

    // list_hole<Polyhedron>(best_outer);
    // std::cout << "Segment " << best_outer->vertex()->point()
    //           << ", " << best_inner->vertex()->point() << std::endl;
    // std::cout << "Quality: " << best_quality << std::endl;

    dolfin_assert(P.is_valid(false, 0));
    const Halfedge_handle new_diagonal = P.split_facet(best_inner, best_outer);
    const Point_3& p = new_diagonal->opposite()->vertex()->point();

    const Vector_3 new_edge(new_diagonal->opposite()->vertex()->point(),
                            new_diagonal->vertex()->point());

    const int num_segments = static_cast<int>(sqrt(CGAL::to_double(new_edge.squared_length())/max_squared_edge_length)+.5);

    // Note: Don't use std::size_t as 0-1 becomes very large...
    for (int i = 1; i < num_segments; i++)
    {
      //std::cout << "Splitting segment" << std::endl;
      Halfedge_handle new_segment = P.split_edge(new_diagonal);
      new_segment->vertex()->point() = p + static_cast<double>(i)/num_segments * new_edge;
    }

    P.normalize_border();
    dolfin_assert(P.is_valid(false, 0));

    return new_diagonal;
  }
  //-----------------------------------------------------------------------------
  /* template<typename Polyhedron> */
  /* static double evaluate_heuristic(const Polyhedron& P, */
  /*                                  typename Polyhedron::Halfedge_handle h, */
  /*                                  double plane_fit) */
  /* { */
  /*   typedef typename Polyhedron::Traits::Triangle_3 Triangle_3; */
  /*   typedef typename Polyhedron::Traits::Vector_3 Vector_3; */
  /*   typedef CGAL::Exact_predicates_inexact_constructions_kernel InexactKernel; */
  /*   typedef typename InexactKernel::Plane_3 InexactPlane_3; */


  /*   // const double distance_to_plane_weight = 1.0; */
  /*   const double planarity_weight = 1.0; */
  /*   const double dihedral_weight  = 1.0; */
  /*   const double ear_angle_weight = 1.0; */

  /*   // Compute the planarity of the points excluding the current point */
  /*   InexactPlane_3 p; */
  /*   const double plane_fit_quality = get_plane_fit<Polyhedron>(h->next(), */
  /*                                                              h->prev(), */
  /*                                                              &p); */
  /*   // Compute the maximum of the dihedral angle to the neighbors */
  /*   const Triangle_3 candidate_triangle(h->prev()->vertex()->point(), */
  /*                                       h->vertex()->point(), */
  /*                                       h->next()->vertex()->point()); */
  /*   const double cos_dihedral = (std::min(get_triangle_cos_angle(candidate_triangle, */
  /*                                                                get_facet_triangle<Polyhedron>(h->opposite())), */
  /*                                         get_triangle_cos_angle(candidate_triangle, */
  /*                                                                get_facet_triangle<Polyhedron>(h->next()->opposite())))+1)/2; */

  /*   // Compute the angle of the cutted ear */
  /*   const Vector_3 v1(h->vertex()->point(), */
  /*                     h->prev()->vertex()->point()); */
  /*   const Vector_3 v2(h->vertex()->point(), */
  /*                     h->next()->vertex()->point()); */

  /*   const double cos_ear_angle = (CGAL::to_double((v1*v2)/std::sqrt(CGAL::to_double(v1.squared_length()*v2.squared_length())))+1)/2.0; */
  /*   const double ear_angle_quality = cos_ear_angle; */

  /*   std::cout << "Triangle " << candidate_triangle[0] */
  /*             << ", " << candidate_triangle[1]  */
  /*             << "," << candidate_triangle[2] << std::endl; */

  /*   std::cout << "Evaluate: planarity: " << (plane_fit_quality/plane_fit) */
  /*             << ", dihedral: " << cos_dihedral */
  /*             << ", ear angle: " << ear_angle_quality << std::endl; */

  /*   return planarity_weight*plane_fit_quality + dihedral_weight*cos_dihedral + ear_angle_quality*ear_angle_weight; */

  /*   /\* return planarity_weight*plane_fit_quality/plane_fit + *\/ */
  /*   /\*   dihedral_weight*cos_dihedral + *\/ */
  /*   /\*   (1-std::exp(-ear_angle_weight*cos_ear_angle))*cos_dihedral; *\/ */
  /* } */
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool vertex_is_border(typename Polyhedron::Vertex_const_handle v)
  {
    //typename Polyhedron::Vertex::Halfedge_around_vertex_circulator h_start = v->vertex_begin();
    auto h_start = v->vertex_begin();
    auto h_current = h_start;
    do
    {
      if (h_current->is_border_edge())
        return true;
      h_current++;
    } while(h_current != h_start);

    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool facets_are_neighbors(typename Polyhedron::Facet_handle f1,
                                   typename Polyhedron::Facet_handle f2)
  {
    typename Polyhedron::Halfedge_handle h1 = f1->halfedge();
    typename Polyhedron::Halfedge_handle start1 = h1;
    do
    {
      typename Polyhedron::Halfedge_handle h2 = f2->halfedge();
      typename Polyhedron::Halfedge_handle start2 = h2;
      do
      {
        if (h2->vertex() == h1->vertex())
          return true;

        h2 = h2->next();
      } while (h2 != start2);

      h1 = h1->next();
    } while (h1 != start1);

    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static bool segment_intersects_facets(typename Polyhedron::Vertex_handle v1,
                                        typename Polyhedron::Vertex_handle v2)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_circulator Vertex_circulator;
    typedef typename Polyhedron::Traits::Segment_3 Segment_3;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;
        typedef typename Polyhedron::Traits::Point_3 Point_3;

    Segment_3 s(v1->point(), v2->point());

    Vertex_circulator start = v1->vertex_begin();
    Vertex_circulator current = start;
    do
    {
      if (!current->is_border() && current->facet()->is_triangle())
      {
        Triangle_3 t = get_facet_triangle<Polyhedron>(current);
        auto result = CGAL::intersection(t, s);
        dolfin_assert(result);

        if (const Point_3* p = boost::get<Point_3>(&*result))
        {
          dolfin_assert(*p == s.source() || *p == s.target());
        }
        else if (const Segment_3* s = boost::get<Segment_3>(&*result))
        {
          return true;
        }
        else
        {
          dolfin_assert(false);
        }
      }

      current++;
    } while (start != current);

    return false;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::string list_self_intersections(Polyhedron& p)
  {
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits Polyhedron_traits;
    typedef typename Polyhedron_traits::Triangle_3 Triangle_3;
    typedef typename Polyhedron_traits::Point_3 Point_3;
    typedef typename Polyhedron_traits::Segment_3 Segment_3;

    std::stringstream ss;

    std::vector<std::pair<Facet_handle, Facet_handle> > intersections;
    CGAL::Polygon_mesh_processing::self_intersections(p, std::back_inserter(intersections));

    for (const std::pair<Facet_handle, Facet_handle>& iit : intersections)
    {
      ss << "Intersection (neighbors: " << (facets_are_neighbors<Polyhedron>(iit.first, iit.second) ? "Yes" : "No") << ")" << std::endl;
      const Halfedge_handle h1 = iit.first->halfedge();
      const Halfedge_handle h2 = iit.second->halfedge();
      print_triangle<Polyhedron>(h1);
      print_triangle<Polyhedron>(h2);

      // Compute intersection

      const Triangle_3 t1(h1->vertex()->point(),
                          h1->next()->vertex()->point(),
                          h1->next()->next()->vertex()->point());

      const Triangle_3 t2(h2->vertex()->point(),
                          h2->next()->vertex()->point(),
                          h2->next()->next()->vertex()->point());

      const auto result = CGAL::intersection(t1, t2);

      dolfin_assert(result);
      if (const Segment_3* s = boost::get<Segment_3>(&*result))
      {
        ss << "Segment: " << *s << std::endl;
      }
      else if (const Point_3* p = boost::get<Point_3>(&*result))
      {
        ss << "Point: " << *p << std::endl;
      }
      else if (const Triangle_3* t = boost::get<Triangle_3>(&*result))
      {
        ss << "Triangle: " << *t << std::endl;
      }
      else
      {
        ss << "Polygon" << std::endl;
      }
    }

    return ss.str();
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static double closest_vertices(const Polyhedron& p)
  {
    // TODO: Use a better optimal algorithm for closest pair problem
    std::cout << "Computing closest vertices" << std::endl;
    double min_distance = std::numeric_limits<double>::max();

    std::size_t counter = 0;
    std::cout << "Vertices: " << p.size_of_vertices() << std::endl;
    for (typename Polyhedron::Vertex_const_iterator v1 = p.vertices_begin();
         v1 != p.vertices_end(); v1++)
    {
      if (counter % 1000 == 0)
        std::cout << counter << std::endl;

      typename Polyhedron::Vertex_const_handle v2 = v1;
      v2++;
      std::size_t counter2 = 0;
      for (;v2 != p.vertices_end(); v2++)
      {
        min_distance = std::min(min_distance,
                                CGAL::to_double(CGAL::squared_distance(v1->point(), v2->point())));

        if (min_distance == 0)
          return 0.0;

        counter2++;
      }

      counter++;
    }

    std::cout << "  Done computing closest" << std::endl;
    return min_distance;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::size_t min_vertex_degree(const Polyhedron& p)
  {
    std::size_t min_degree = std::numeric_limits<std::size_t>::max();
    for (typename Polyhedron::Vertex_const_iterator it = p.vertices_begin();
         it != p.vertices_end(); it++)
    {
      if (!vertex_is_border<Polyhedron>(it))
      {
        min_degree = std::min(min_degree, it->vertex_degree());
      }
    }
    return min_degree;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void remove_vertex(Polyhedron& P, typename Polyhedron::Vertex_handle v)
  {
    typedef typename Polyhedron::Halfedge_around_vertex_circulator Vertex_circulator;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

    // std::cout << "Removing vertex" << std::endl;

    Vertex_circulator h = v->vertex_begin();
    Vertex_circulator start = h;

    std::vector<Halfedge_handle> to_be_removed;

    do
    {
      if (!h->is_border())
        to_be_removed.push_back(h);

      h++;
    } while (h != start);


    // std::cout << "Removing " << to_be_removed.size() << " halfedges" << std::endl;
    for (auto it = to_be_removed.begin(); it != to_be_removed.end(); it++)
    {
      P.erase_facet(*it);
    }

    // std::cout << "  done removing vertex" << std::endl;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static std::size_t remove_self_intersections(Polyhedron& P)
  {
    std::size_t removed = 0;

    // typedef typename Polyhedron::Traits Polyhedron_traits;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    std::vector<std::pair<Facet_handle, Facet_handle> > intersections;
    CGAL::Polygon_mesh_processing::self_intersections(P, std::back_inserter(intersections));

    while (intersections.size() > 0)
    {
      std::cout << "Removing self intersection (" << intersections.size() << ")" << std::endl;
      const typename Polyhedron::Facet_handle f1 = intersections.front().first;
      const typename Polyhedron::Facet_handle f2 = intersections.front().second;

      std::deque<Facet_handle> queue1;
      queue1.push_back(f1);

      std::deque<Facet_handle> queue2;
      queue2.push_back(f2);

      std::set<Facet_handle> to_be_removed1;
      to_be_removed1.insert(f1);

      std::set<Facet_handle> to_be_removed2;
      to_be_removed2.insert(f2);

      bool f1_done = false;

      while (!f1_done)
      {
        // Pop from queue 1
        if (queue1.size() > 0)
        {
          to_be_removed1.insert(queue1.front());
          Halfedge_handle start = queue1.front()->halfedge();

          queue1.pop_front();
          Halfedge_handle current = start;

          do
          {
            std::cout << "Spreading out 1" << std::endl;
            if (!current->is_border_edge())
            {
              //if (to_be_removed2.count(current->opposite()->facet()) > 0)
              if (current->opposite()->facet() == f2)
              {
                f1_done = true;
                break;
              }
              else
              {
                queue1.push_back(current->opposite()->facet());
              }
            }
            current = current->next();
          } while (current != start);
        }
        else
        {
          f1_done = true;
        }
      }

      bool f2_done = false;

      while (!f2_done)
      {
        // Pop from queue 2
        if (queue2.size() > 0)
        {
          to_be_removed2.insert(queue2.front());
          Halfedge_handle start = queue2.front()->halfedge();
          queue2.pop_front();
          Halfedge_handle current = start;

          do
          {
            std::cout << "Spreading out 2" << std::endl;
            if (!current->is_border_edge())
            {
              // if (to_be_removed1.count(current->opposite()->facet()) > 0)
              if (current->opposite()->facet() == f1)
              {
                f2_done = true;
                break;
              }
              else
              {
                queue2.push_back(current->opposite()->facet());
              }
            }
            current = current->next();
          } while (current != start);
        }
        else
        {
          f2_done = true;
        }
      }

      std::cout << "To be removed 1: " << to_be_removed1.size() << std::endl;
      std::cout << "To be removed 2: " << to_be_removed2.size() << std::endl;

      to_be_removed1.insert(to_be_removed2.begin(),
                            to_be_removed2.end());

      for (typename Polyhedron::Face_handle f : to_be_removed1)
      {
        P.erase_facet(f->halfedge());
        dolfin_assert(P.is_valid());
        removed++;
      }


      /* for (typename Polyhedron::Face_handle f : to_be_removed2) */
      /* { */
      /*   P.erase_facet(f->halfedge()); */
      /*   dolfin_assert(P.is_valid()); */
      /*   removed++; */
      /* } */

      dolfin_assert(P.is_valid());

      intersections.clear();
      CGAL::Polygon_mesh_processing::self_intersections(P, std::back_inserter(intersections));

      break;
    }

    return removed;
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void filter_sharp_features(Polyhedron& P, int start_facet, double tolerance)
  {
    typedef typename Polyhedron::Facet_iterator Facet_iterator;
    typedef typename Polyhedron::Facet_handle Facet_handle;
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Triangle_3 Triangle_3;

    /* { */
    /*   std::ofstream outfile("before_filtering.off"); */
    /*   outfile << P; */
    /* } */

    /* std::cout << "Filter sharp features" << std::endl; */

    const double cos_tolerance = std::cos(tolerance);
    // std::cout << "tolerance: " << cos_tolerance << std::endl;
    Facet_iterator fit = P.facets_begin();
    for (int i = 0; i < start_facet; i++)
      fit++;

    std::deque<Facet_handle> queue;
    {
      Halfedge_handle h = fit->halfedge();
      // std::cout << "Starting facet: Triangle " << h->vertex()->point() << ", ";
      h = h->next();
      // std::cout << h->vertex()->point() << ", ";
      h = h->next();
      // std::cout << h->vertex()->point() << std::endl;
    }

    std::set<Facet_handle> visited;
    std::set<Facet_handle> to_be_removed;
    for (Facet_iterator fit = P.facets_begin(); fit != P.facets_end(); fit++)
      to_be_removed.insert(fit);

    // std::cout << "Number of facets: " << P.size_of_facets() << std::endl;

    queue.push_back(fit);
    while (!queue.empty())
    {
      // std::cout << "In queue" << std::endl;
      Facet_handle f = queue.front();
      queue.pop_front();

      if (visited.count(f) > 0)
      {
        // std::cout << "Already handled" << std::endl;
        continue;
      }

      visited.insert(f);
      to_be_removed.erase(f);

      const Halfedge_handle start = f->halfedge();
      Halfedge_handle current = start;
      do
      {
        // std::cout << "Exploring neighbor" << std::endl;
        if (!current->opposite()->is_border())
        {

          Triangle_3 t1 = get_facet_triangle<Polyhedron>(current);
          Triangle_3 t2 = get_facet_triangle<Polyhedron>(current->opposite());
          if (get_triangle_cos_angle(t1, t2) > cos_tolerance)
          {
            queue.push_back(current->opposite()->facet());
          }
        }

        current = current->next();
      } while (current != start);
    }

    // std::cout << "Remove " << to_be_removed.size() << " facets" << std::endl;

    for (auto fit = to_be_removed.begin(); fit != to_be_removed.end(); fit++)
    {
      P.erase_facet( (*fit)->halfedge() );
    }

    /* { */
    /*   std::ofstream outfile("after_filtering.off"); */
    /*   outfile << P; */
    /* } */
  }

//-----------------------------------------------------------------------------
  template <typename Polyhedron>
  static std::pair<typename Polyhedron::Traits::Point_3, typename Polyhedron::Traits::Point_3>
  get_aabb(const Polyhedron& P)
  {
    typedef typename Polyhedron::Vertex_const_iterator Vertex_const_iterator;
    typedef typename Polyhedron::Traits::Point_3 Point_3;
    typedef typename Polyhedron::Traits::FT FT;
    Vertex_const_iterator it = P.vertices_begin();

    FT x_min = it->point().x();
    FT y_min = it->point().y();;
    FT z_min = it->point().z();
    FT x_max = x_min;
    FT y_max = y_min;
    FT z_max = z_min;
    it++;

    for (;it != P.vertices_end(); it++)
    {
      const Point_3& current = it->point();
      x_min = std::min(x_min, current.x());
      y_min = std::min(y_min, current.y());
      z_min = std::min(z_min, current.z());

      x_max = std::max(x_max, current.x());
      y_max = std::max(y_max, current.y());
      z_max = std::max(z_max, current.z());
    }

    return std::make_pair(Point_3(x_min, y_min, z_min),
                          Point_3(x_max, y_max, z_max));
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static double sharpest_edge(const Polyhedron& P)
  {
    typedef typename Polyhedron::Traits Kernel;
    typedef typename Polyhedron::Edge_const_iterator Edge_const_iterator;
    typedef typename Kernel::Vector_3 Vector_3;

    double min_cos = 1.0;
    for (Edge_const_iterator it = P.edges_begin(); it != P.edges_end(); it++)
    {
      //const Point_3 a = it->vertex()->point();
      const Vector_3 n1 = CGAL::normal<Kernel>(it->vertex()->point(),
                                               it->next()->vertex()->point(),
                                               it->next()->next()->vertex()->point());
      const Vector_3 n2 = CGAL::normal<Kernel>(it->opposite()->vertex()->point(),
                                               it->opposite()->next()->vertex()->point(),
                                               it->opposite()->next()->next()->vertex()->point());
      min_cos = std::min(min_cos,
                         CGAL::to_double(n1*n2)/sqrt(CGAL::to_double(n1.squared_length()*n2.squared_length())));
    }

    return acos(min_cos);
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void center_vertex_refinement(Polyhedron& P)
  {
    // Precollect facet handles since we are adding facets within the loop
    std::vector<typename Polyhedron::Facet_handle> facets;
    facets.reserve(P.size_of_facets());
    for (typename Polyhedron::Facet_iterator fit = P.facets_begin();
         fit != P.facets_end();
         fit++)
    {
      facets.push_back(fit);
    }

    for (typename Polyhedron::Facet_handle f : facets)
    {
      typename Polyhedron::Halfedge_handle h = f->halfedge();
      typename Polyhedron::Traits::Point_3 c = CGAL::centroid(h->vertex()->point(),
                                                              h->next()->vertex()->point(),
                                                              h->next()->next()->vertex()->point());
      c.exact();
      typename Polyhedron::Halfedge_handle h_new = P.create_center_vertex(h);
      h_new->vertex()->point() = c;
    }
  }
  //-----------------------------------------------------------------------------
  template<typename Polyhedron>
  static void split_edge_refinement(Polyhedron& P)
  {
    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Traits::Point_3 Point_3;

    // Collect vertices from the original
    std::set<Halfedge_handle> facets;
    // facets.reserve(P.size_of_facets());
    for (typename Polyhedron::Facet_iterator fit = P.facets_begin();
         fit != P.facets_end();
         fit++)
    {
      const Halfedge_handle h = fit->halfedge();
      facets.insert(h);
    }

    // Split all edges
    {
      std::vector<Halfedge_handle> edges;
      for (typename Polyhedron::Edge_iterator eit = P.edges_begin();
           eit != P.edges_end(); eit++)
      {
        edges.push_back(eit);
      }

      for (Halfedge_handle h : edges)
      {
        const Point_3 p = CGAL::midpoint(h->vertex()->point(),
                                         h->opposite()->vertex()->point());
        Halfedge_handle opposite = h->opposite();
        Halfedge_handle h_new = P.split_edge(h);
        h_new->vertex()->point() = p;

        // This is suboptimal, but split_edge(h) changes h->opposite() to point
        // to the newly inserted vertex, so we have to check and update if
        // h->opposite() is stored in facets
        if (facets.erase(opposite))
        {
          facets.insert(h_new->opposite());
        }
      }
    }

    for (Halfedge_handle h : facets)
    {
      Halfedge_handle h2 = h->next()->next();
      Halfedge_handle h3 = h2->next()->next();
      P.split_facet(h->next(), h->prev());
      P.split_facet(h2->next(), h2->prev());
      P.split_facet(h3->next(), h3->prev());
    }
  }
};
  //-----------------------------------------------------------------------------
// Taken from demo/Polyhedron/Scene_nef_polyhedron_item.cpp in the
// CGAL source tree.
// Quick hacks to convert polyhedra from exact to inexact and
// vice-versa
template <class Polyhedron_input, class Polyhedron_output>
struct Copy_polyhedron_to
  : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
  Copy_polyhedron_to(const Polyhedron_input& in_poly)
    : _in_poly(in_poly) {}

  void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
  {
    typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
    //typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

    CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

    typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
    typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
    typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

    builder.begin_surface(_in_poly.size_of_vertices(),
      _in_poly.size_of_facets(),
      _in_poly.size_of_halfedges());

    for(Vertex_const_iterator
      vi = _in_poly.vertices_begin(), end = _in_poly.vertices_end();
      vi != end ; ++vi)
    {
      typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
	::CGAL::to_double( vi->point().y()),
	::CGAL::to_double( vi->point().z()));
      builder.add_vertex(p);
    }

    typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
    Index index(_in_poly.vertices_begin(), _in_poly.vertices_end());

    for(Facet_const_iterator
      fi = _in_poly.facets_begin(), end = _in_poly.facets_end();
      fi != end; ++fi)
    {
      HFCC hc = fi->facet_begin();
      HFCC hc_end = hc;
      //     std::size_t n = circulator_size( hc);
      //     CGAL_assertion( n >= 3);
      builder.begin_facet ();
      do
      {
        builder.add_vertex_to_facet(index[hc->vertex()]);
        ++hc;
      } while( hc != hc_end);
      builder.end_facet();
    }
    builder.end_surface();
  } // end operator()(..)
private:
  const Polyhedron_input& _in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_A, class Poly_B>
void copy_to(const Poly_A& poly_a, Poly_B& poly_b)
{
  Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
  poly_b.delegate(modifier);
  // CGAL_assertion(poly_b.is_valid());
}
}

#endif
