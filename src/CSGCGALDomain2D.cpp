// Copyright (C) 2013-2017 Benjamin Kehlet
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

// This file must be included to get the compiler flags
// They should ideally have been added via the command line,
// but since CGAL configure time is at mshr compile time, we don't
// have access to CGALConfig.cmake at configure time...
#ifndef CGAL_HEADER_ONLY
#include <CGAL/compiler_config.h>
#endif
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <mshr/CSGCGALDomain2D.h>
#include <mshr/CSGPrimitives2D.h>
#include <mshr/CSGOperators.h>

#include <dolfin/common/constants.h>
#include <dolfin/log/LogStream.h>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_set_2.h>

#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>

// Polygon typedefs
//typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef CGAL::Quotient<CGAL::MP_Float>                    FT;
typedef CGAL::Cartesian<FT>                               Exact_Kernel;

typedef Exact_Kernel::Point_2                             Point_2;
typedef Exact_Kernel::Vector_2                            Vector_2;
typedef Exact_Kernel::Segment_2                           Segment_2;
typedef Exact_Kernel::Ray_2                               Ray_2;
typedef Exact_Kernel::Direction_2                         Direction_2;
typedef Exact_Kernel::Aff_transformation_2                Aff_transformation_2;
typedef CGAL::Polygon_2<Exact_Kernel>                     Polygon_2;
typedef Polygon_2::Vertex_const_iterator                  Vertex_const_iterator;
typedef CGAL::Polygon_with_holes_2<Exact_Kernel>          Polygon_with_holes_2;
typedef Polygon_with_holes_2::Hole_const_iterator         Hole_const_iterator;
typedef CGAL::Polygon_set_2<Exact_Kernel>                 Polygon_set_2;

// Min enclosing circle typedefs
typedef CGAL::Min_circle_2_traits_2<Exact_Kernel>  Min_Circle_Traits;
typedef CGAL::Min_circle_2<Min_Circle_Traits>      Min_circle;
typedef CGAL::Circle_2<Exact_Kernel> CGAL_Circle;


namespace
{
FT get_shortest_edge(const Polygon_set_2& polygon_set)
{
  FT shortest_edge = std::numeric_limits<double>::max();

  std::list<Polygon_with_holes_2> polygon_list;
  polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  // FIXME: For now it looks only at the outer boundary.
  // Should take holes into account as well
  for (const Polygon_with_holes_2& p : polygon_list)
  {
    const Polygon_2& bdr = p.outer_boundary();
    Point_2 prev = bdr.container().back();
    for (Polygon_2::Vertex_const_iterator vit = bdr.vertices_begin();
         vit != bdr.vertices_end(); vit++)
    {
      shortest_edge = std::min(shortest_edge, (*vit-prev).squared_length());
      prev = *vit;
    }
  }
  return shortest_edge;
}
//-----------------------------------------------------------------------------
Point_2 point_in_polygon(const Polygon_2& p)
{
  // Take the midpoint, s,  of a segment, then do a ray shooting in the inwards
  // normal direction and return the midpoint between s and the closest ray hit.

  Polygon_2::Edge_const_iterator eit = p.edges_begin();
  const Point_2 source = CGAL::midpoint(eit->source(), eit->target());
  const Aff_transformation_2 tr = Aff_transformation_2(CGAL::ROTATION,
                                                       p.orientation() == CGAL::COUNTERCLOCKWISE ? 1 : -1, 0);

  const Ray_2 r(source, tr(Direction_2(*eit)));
  Point_2 closest;
  FT min_squared_distance = std::numeric_limits<double>::max();

  eit++;
  for (; eit != p.edges_end(); eit++)
  {
    auto ii = CGAL::intersection(r, *eit);
    if (ii)
    {
      if (const Point_2* p = boost::get<Point_2>(&*ii))
      {
        const FT squared_distance = CGAL::squared_distance(source, *p);
        if (squared_distance < min_squared_distance)
        {
          closest = *p;
          min_squared_distance = squared_distance;
        }
      }
      else if (const Segment_2* s = boost::get<Segment_2>(&*ii))
      {
        // If the intersection is a segment, an edge in the polygon is paralell
        // to the ray. Simply check both the source and the target point of the
        // segment.

        {
          const FT squared_distance = CGAL::squared_distance(source, s->source());
          if (squared_distance < min_squared_distance)
          {
            closest = *p;
            min_squared_distance = squared_distance;
          }
        }

        {
          const FT squared_distance = CGAL::squared_distance(source, s->target());
          if (squared_distance < min_squared_distance)
          {
            closest = *p;
            min_squared_distance = squared_distance;
          }
        }
      } // end intersection is Segment_2
    }
  }
  return CGAL::midpoint(source, closest);
}
//-----------------------------------------------------------------------------
void truncate_short_edges(Polygon_set_2& polygon_set,
                          FT tolerance,
                          const std::set<Point_2> & collapsable_vertices = std::set<Point_2>())
{
  const FT tolerance_squared = tolerance*tolerance;

  Polygon_set_2 truncated_set;
  std::list<Polygon_with_holes_2> polygon_list;
  polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  // FIXME: For now it looks only at the outer boundary.
  // Should take holes into account as well
  for (const Polygon_with_holes_2& p : polygon_list)
  {
    Polygon_2 bdr = p.outer_boundary();

    bool changed = true;
    while(changed)
    {
      changed = false;

      Polygon_2::Vertex_iterator prev = bdr.vertices_end()-1;

      for (Polygon_2::Vertex_iterator vit = bdr.vertices_begin();
           vit != bdr.vertices_end(); vit++)
      {
        if ( (*vit-*prev).squared_length() < tolerance_squared)
        {
          if (collapsable_vertices.size() == 0 || collapsable_vertices.find(*vit) != collapsable_vertices.end())
          {
            bdr.erase(vit);
            changed = true;
            break;
          }
          else if (collapsable_vertices.size() == 0 || collapsable_vertices.find(*prev) != collapsable_vertices.end())
          {
            bdr.erase(prev);
            changed = true;
            break;
          }
          else
          {
            // std::cout << "Couldn't erase vertex" << std::endl;
          }
        } // end if short edge
        prev = vit;
      } // end vertex iteration
    } // end while(changed)

    Polygon_with_holes_2 pwh(bdr);
    for (Hole_const_iterator hit = p.holes_begin(); hit != p.holes_end(); ++hit)
    {
      pwh.add_hole(*hit);
    }
    truncated_set.insert(pwh);
  }
  polygon_set = truncated_set;
}
} // end anonymous namespace
//-----------------------------------------------------------------------------
namespace mshr
{

struct CSGCGALDomain2DImpl
{
  Polygon_set_2 polygon_set;

  CSGCGALDomain2DImpl(){}
  explicit CSGCGALDomain2DImpl(const Polygon_set_2& p)
    : polygon_set(p) {}
};
//-----------------------------------------------------------------------------
Polygon_2 make_circle(
    const Circle* c,
    const double segment_granularity
    )
{
  unsigned int num_segments = 0;
  if (c->segments() > 0)
  {
    num_segments = c->segments();
  } else {
    dolfin_assert(segment_granularity > 0.0);

    num_segments = std::round(
        (2 * DOLFIN_PI * c->radius()) / segment_granularity
        );
  }

  // Set the minimum number of segments to 5 (to be a bit circle-like)
  // (A polygon with less than 3 segments is degenerate)
  num_segments = std::max(num_segments, 5u);

  std::vector<Point_2> pts;
  pts.reserve(num_segments);

  for (std::size_t i = 0; i < num_segments; i++)
  {
    const double phi = (2*DOLFIN_PI*i) / num_segments;
    const double x = c->center().x() + c->radius()*cos(phi);
    const double y = c->center().y() + c->radius()*sin(phi);
    pts.push_back(Point_2(x, y));
  }

  return Polygon_2(pts.begin(), pts.end());
}
//-----------------------------------------------------------------------------
Polygon_2 make_ellipse(const Ellipse* e,
                       const double segment_granularity)
{
  unsigned int num_segments = 0;
  if (e->segments() > 0) {
    num_segments = e->segments();
  }
  else
  {
    dolfin_assert(segment_granularity > 0.0);
    // https://en.wikipedia.org/wiki/Ellipse#Circumference
    const double a_min_b = e->a() - e->b();
    const double a_plus_b = e->a() + e->b();
    const double h = (a_min_b*a_min_b) / (a_plus_b*a_plus_b);
    const double arc_length_ellipse_approx = DOLFIN_PI * a_plus_b * (
        1.0 + 3.0*h / (10 + sqrt(4.0 - 3.0*h))
        );
    num_segments = std::round(arc_length_ellipse_approx / segment_granularity);
  }

  // A polygon with less segments than 3 is degenerate
  // FIXME: Should this result in an empty polygon?
  num_segments = std::max(num_segments, 3u);

  std::vector<Point_2> pts;
  pts.reserve(num_segments);

  for (std::size_t i = 0; i < num_segments; i++)
  {
    const double phi = (2*DOLFIN_PI*i) / num_segments;
    const double x = e->center().x() + e->a()*cos(phi);
    const double y = e->center().y() + e->b()*sin(phi);
    pts.push_back(Point_2(x, y));
  }

  return Polygon_2(pts.begin(), pts.end());
}
//-----------------------------------------------------------------------------
Polygon_2 make_rectangle(const Rectangle* r)
{
  const double x0 = std::min(r->first_corner().x(), r->second_corner().x());
  const double y0 = std::min(r->first_corner().y(), r->second_corner().y());

  const double x1 = std::max(r->first_corner().x(), r->second_corner().x());
  const double y1 = std::max(r->first_corner().y(), r->second_corner().y());

  std::vector<Point_2> pts;
  pts.push_back(Point_2(x0, y0));
  pts.push_back(Point_2(x1, y0));
  pts.push_back(Point_2(x1, y1));
  pts.push_back(Point_2(x0, y1));

  Polygon_2 p(pts.begin(), pts.end());

  return p;
}
//-----------------------------------------------------------------------------
Polygon_2 make_polygon(const Polygon* p)
{
  std::vector<Point_2> pts;
  std::vector<dolfin::Point>::const_iterator v;
  for (v = p->vertices().begin(); v != p->vertices().end(); ++v)
    pts.push_back(Point_2(v->x(), v->y()));

  return Polygon_2(pts.begin(), pts.end());
}
//-----------------------------------------------------------------------------
std::unique_ptr<CSGCGALDomain2DImpl> do_transformation(const Polygon_set_2& p, Exact_Kernel::Aff_transformation_2 t)
{
  std::unique_ptr<CSGCGALDomain2DImpl> result(new CSGCGALDomain2DImpl);

  std::list<Polygon_with_holes_2> polygon_list;
  p.polygons_with_holes(std::back_inserter(polygon_list));

  std::list<Polygon_with_holes_2>::const_iterator pit;
  for (pit = polygon_list.begin(); pit != polygon_list.end(); ++pit)
  {
    const Polygon_with_holes_2& pwh = *pit;

    // Transform outer boundary
    Polygon_with_holes_2 transformed(CGAL::transform(t, pwh.outer_boundary()));

    // Transform holes
    for (Hole_const_iterator hit = pwh.holes_begin(); hit != pwh.holes_end(); hit++)
    {
      transformed.add_hole(CGAL::transform(t, *hit));
    }

    result->polygon_set.insert(transformed);
  }

  return result;
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::CSGCGALDomain2D()
  : impl(new CSGCGALDomain2DImpl)
{

}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::~CSGCGALDomain2D()
{
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::CSGCGALDomain2D(std::shared_ptr<const CSGGeometry> geometry,
                                 double segment_granularity)
 : impl(new CSGCGALDomain2DImpl)
{
  if (geometry->dim() != 2)
    dolfin::dolfin_error("CSGCGALDomain2D.cpp",
                         "Creating polygonal domain",
                         "Geometry has dimension %d, expected 2", geometry->dim());


  switch (geometry->getType())
  {
    case CSGGeometry::Union:
    {
      std::shared_ptr<const CSGUnion> u = std::dynamic_pointer_cast<const CSGUnion, const CSGGeometry>(geometry);
      dolfin_assert(u);

      CSGCGALDomain2D a(u->_g0, segment_granularity);
      CSGCGALDomain2D b(u->_g1, segment_granularity);

      impl.swap(a.impl);
      impl->polygon_set.join(b.impl->polygon_set);
      break;
    }
    case CSGGeometry::Intersection:
    {
      auto  u = std::dynamic_pointer_cast<const CSGIntersection, const CSGGeometry>(geometry);
      dolfin_assert(u);

      CSGCGALDomain2D a(u->_g0, segment_granularity);
      CSGCGALDomain2D b(u->_g1, segment_granularity);

      impl.swap(a.impl);
      impl->polygon_set.intersection(b.impl->polygon_set);
      break;
    }
    case CSGGeometry::Difference:
    {
      auto u = std::dynamic_pointer_cast<const CSGDifference, const CSGGeometry>(geometry);
      dolfin_assert(u);
      CSGCGALDomain2D a(u->_g0, segment_granularity);
      CSGCGALDomain2D b(u->_g1, segment_granularity);

      impl.swap(a.impl);
      impl->polygon_set.difference(b.impl->polygon_set);
      break;
    }
    case CSGGeometry::Translation :
    {
      auto t = std::dynamic_pointer_cast<const CSGTranslation, const CSGGeometry>(geometry);
      dolfin_assert(t);
      CSGCGALDomain2D a(t->g, segment_granularity);
      Exact_Kernel::Aff_transformation_2 translation(CGAL::TRANSLATION, Vector_2(t->t.x(), t->t.y()));
      std::unique_ptr<CSGCGALDomain2DImpl> transformed = do_transformation(a.impl->polygon_set, translation);
      impl.swap(transformed);
      break;
    }
    case CSGGeometry::Scaling :
    {
      auto t = std::dynamic_pointer_cast<const CSGScaling, const CSGGeometry>(geometry);
      dolfin_assert(t);
      CSGCGALDomain2D a(t->g, segment_granularity);
      Exact_Kernel::Aff_transformation_2 tr(CGAL::IDENTITY);

      // Translate if requested
      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2 (CGAL::TRANSLATION,
                                                 Vector_2(-t->c.x(), -t->c.y())) * tr;

      // Do the scaling
      tr = Exact_Kernel::Aff_transformation_2(CGAL::SCALING, t->s) * tr;

      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2(CGAL::TRANSLATION,
                                                Vector_2(t->c.x(), t->c.y())) * tr;

      std::unique_ptr<CSGCGALDomain2DImpl> transformed = do_transformation(a.impl->polygon_set,
                                                                           tr);
      impl.swap(transformed);
      break;
    }
    case CSGGeometry::Rotation :
    {
      auto t = std::dynamic_pointer_cast<const CSGRotation, const CSGGeometry>(geometry);
      dolfin_assert(t);
      CSGCGALDomain2D a(t->g, segment_granularity);
      Exact_Kernel::Aff_transformation_2 tr(CGAL::IDENTITY);

      // Translate if requested
      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2 (CGAL::TRANSLATION,
                                                 Vector_2(-t->c.x(), -t->c.y())) * tr;

      // Do the rotation
      tr = Exact_Kernel::Aff_transformation_2(CGAL::ROTATION, sin(t->theta), cos(t->theta)) * tr;

      if (t->translate)
        tr = Exact_Kernel::Aff_transformation_2(CGAL::TRANSLATION,
                                                Vector_2(t->c.x(), t->c.y())) * tr;

      std::unique_ptr<CSGCGALDomain2DImpl> transformed = do_transformation(a.impl->polygon_set,
                                                                           tr);
      impl.swap(transformed);
      break;
    }
    case CSGGeometry::Circle:
    {
      auto c = std::dynamic_pointer_cast<const Circle, const CSGGeometry>(geometry);
      dolfin_assert(c);
      impl->polygon_set.insert(make_circle(c.get(), segment_granularity));
      break;
    }
    case CSGGeometry::Ellipse:
    {
      auto c = std::dynamic_pointer_cast<const Ellipse, const CSGGeometry>(geometry);
      dolfin_assert(c);
      impl->polygon_set.insert(make_ellipse(c.get(), segment_granularity));
      break;
    }
    case CSGGeometry::Rectangle:
    {
      auto r = std::dynamic_pointer_cast<const Rectangle, const CSGGeometry>(geometry);
      dolfin_assert(r);
      impl->polygon_set.insert(make_rectangle(r.get()));
      break;
    }
    case CSGGeometry::Polygon:
    {
      auto p = std::dynamic_pointer_cast<const Polygon, const CSGGeometry>(geometry);
      dolfin_assert(p);
      impl->polygon_set.insert(make_polygon(p.get()));
      break;
    }
    default:
      dolfin::dolfin_error("CSGCGALMeshGenerator2D.cpp",
                           "converting geometry to cgal polyhedron",
                           "Unhandled primitive type");
  }

  // Truncate short edges
  // std::cout << "Creating domain, shortest edge: " << CGAL::to_double(get_shortest_edge(impl->polygon_set)) << std::endl;

  truncate_short_edges(impl->polygon_set, 1e-15);
  // std::cout << "Truncated, shortest edge: " << CGAL::to_double(get_shortest_edge(impl->polygon_set)) << std::endl;
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D::CSGCGALDomain2D(const CSGCGALDomain2D &other)
 : impl(new CSGCGALDomain2DImpl(other.impl->polygon_set))
{
}
//-----------------------------------------------------------------------------
CSGCGALDomain2D &CSGCGALDomain2D::operator=(const CSGCGALDomain2D &other)
{
  std::unique_ptr<CSGCGALDomain2DImpl> tmp(new CSGCGALDomain2DImpl(other.impl->polygon_set));

  impl.swap(tmp);

  return *this;
}
//-----------------------------------------------------------------------------
double CSGCGALDomain2D::compute_boundingcircle_radius() const
{
  std::list<Polygon_with_holes_2> polygon_list;
  impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  std::vector<Point_2> points;

  for (std::list<Polygon_with_holes_2>::const_iterator pit = polygon_list.begin();
       pit != polygon_list.end(); ++pit)
    for (Polygon_2::Vertex_const_iterator vit = pit->outer_boundary().vertices_begin();
         vit != pit->outer_boundary().vertices_end(); ++vit)
      points.push_back(*vit);

  Min_circle min_circle (points.begin(),
                         points.end(),
                         true); //randomize point order

  return sqrt(CGAL::to_double(min_circle.circle().squared_radius()));
}
//-----------------------------------------------------------------------------
std::size_t CSGCGALDomain2D::num_polygons() const
{
  return impl->polygon_set.number_of_polygons_with_holes();
}
//-----------------------------------------------------------------------------
std::vector<dolfin::Point> CSGCGALDomain2D::get_outer_polygon(std::size_t i) const
{
  std::vector<Polygon_with_holes_2> polygon_list;
  impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  const Polygon_2& polygon = polygon_list[i].outer_boundary();

  std::vector<dolfin::Point> res;
  res.reserve(polygon.size());
  for (std::size_t j = 0; j < polygon.size(); j++)
  {
    const Point_2& p = polygon.vertex(j);
    res.push_back(dolfin::Point(CGAL::to_double(p.x()), CGAL::to_double(p.y())));
  }

  return res;
}
//-----------------------------------------------------------------------------
void CSGCGALDomain2D::join_inplace(const CSGCGALDomain2D& other)
{
  impl->polygon_set.join(other.impl->polygon_set);
}
//-----------------------------------------------------------------------------
void CSGCGALDomain2D::difference_inplace(const CSGCGALDomain2D& other)
{
  impl->polygon_set.difference(other.impl->polygon_set);
}
//-----------------------------------------------------------------------------
void CSGCGALDomain2D::intersect_inplace(const CSGCGALDomain2D &other,
                                        double truncate_tolerance)
{
  if (truncate_tolerance > 0)
  {
    // store
    std::set<Point_2> original_vertices;
    {
      std::list<Polygon_with_holes_2> polygon_list;
      impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

      // FIXME: For now it looks only at the outer boundary.
      // Should take holes into account as well
      for (const Polygon_with_holes_2& p : polygon_list)
      {
        const Polygon_2& bdr = p.outer_boundary();
        for (Polygon_2::Vertex_const_iterator vit = bdr.vertices_begin();
             vit != bdr.vertices_end(); vit++)
        {
          original_vertices.insert(*vit);
        }

      }
    }
    // Do the actual intersection
    impl->polygon_set.intersection(other.impl->polygon_set);

    const FT shortest_edge = get_shortest_edge(impl->polygon_set);
    const FT tolerance_squared = truncate_tolerance*truncate_tolerance;
    if (shortest_edge < tolerance_squared)
    {
      truncate_short_edges(impl->polygon_set,
                           truncate_tolerance,
                           original_vertices);
    } // end
  }
  else
  {
    impl->polygon_set.intersection(other.impl->polygon_set);
  }
}
//-----------------------------------------------------------------------------
bool CSGCGALDomain2D::point_in_domain(dolfin::Point p) const
{
  const Point_2 p_(p.x(), p.y());
  return impl->polygon_set.oriented_side(p_) == CGAL::ON_POSITIVE_SIDE;
}
//-----------------------------------------------------------------------------
std::string CSGCGALDomain2D::str(bool verbose) const
{
  std::stringstream ss;
  const std::size_t num_polygons = impl->polygon_set.number_of_polygons_with_holes();
  ss << "<Polygonal domain with " << num_polygons
     << " outer polygon" << (num_polygons != 1 ? "s" : "") << std::endl;

  if (verbose)
  {
    std::list<Polygon_with_holes_2> polygon_list;
    impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));
    for (const Polygon_with_holes_2& p : polygon_list)
    {
      ss << "  Polygon ";
      const Polygon_2& bdr = p.outer_boundary();
      // ss << "[" << bdr.size() << "] Vertices: ";

      for (Polygon_2::Vertex_const_iterator vit = bdr.vertices_begin();
           vit != bdr.vertices_end(); vit++)
      {
        ss << CGAL::to_double(vit->x()) << " " << CGAL::to_double(vit->y()) << ", ";
      }
      ss << std::endl;

      Point_2 prev = bdr.container().back();
      ss << "Edge lengths:";
      for (Polygon_2::Vertex_const_iterator vit = bdr.vertices_begin();
           vit != bdr.vertices_end(); vit++)
      {
        ss << " "  << sqrt(CGAL::to_double((*vit-prev).squared_length()));
        prev = *vit;
      }
      ss << std::endl;



      ss << "[" << p.number_of_holes() << " holes]" << std::endl;
    }
  }

  ss << ">";
  return ss.str();
}
//-----------------------------------------------------------------------------
typedef std::map<Point_2, std::size_t>::iterator VertexMapIterator;
typedef std::set<std::pair<std::size_t, std::size_t>>::iterator SegmentIterator;
//-----------------------------------------------------------------------------
VertexMapIterator pslg_split_edge(SegmentIterator si,
                                  Point_2 p,
                                  std::map<Point_2, std::size_t>& vertex_map,
                                  std::vector<Point_2>& vertices,
                                  std::set<std::pair<std::size_t, std::size_t>>& segments)
{
  dolfin_assert(CGAL::do_intersect(Segment_2(vertices[si->first], vertices[si->second]), p));
  dolfin_assert(p != vertices[si->first]);
  dolfin_assert(p != vertices[si->second]);

  vertices.push_back(p);
  const std::size_t new_vertex_index = vertex_map.size();
  std::pair<VertexMapIterator, bool> result =
    vertex_map.insert(std::make_pair(p, new_vertex_index));
  dolfin_assert(result.second);

  const std::size_t source_vertex_index = si->first;
  const std::size_t target_vertex_index = si->second;

  segments.erase(si);
  segments.insert(std::make_pair(new_vertex_index, target_vertex_index));
  segments.insert(std::make_pair(source_vertex_index, new_vertex_index));
  return result.first;
}
//-----------------------------------------------------------------------------

// Check if point intersects any edge in the interior
static inline SegmentIterator
pslg_segment_intersection(Point_2 p,
                          const std::vector<Point_2>& vertices,
                          const std::set<std::pair<std::size_t, std::size_t>>& segments)
{
  for (SegmentIterator it = segments.begin(); it != segments.end(); ++it)
  {
    const Segment_2 current(vertices[it->first], vertices[it->second]);
    if (CGAL::do_intersect(p, current))
      return it;
  }
  return segments.end();
}

//-----------------------------------------------------------------------------
static inline void pslg_add_simple_polygon(std::map<Point_2, std::size_t>& vertex_map,
                                           std::vector<Point_2>& vertices,
                                           std::set<std::pair<std::size_t, std::size_t>>& segments,
                                           const Polygon_2& p)
{
  // Note that this function assumes that all edges are sufficiently long. Short
  // edges should be filtered out when preparing the domain (and subdomain)
  // polygon_set_2s.

  Point_2 prev = p.vertex(p.size()-1);

  for (std::size_t i = 0; i < p.size(); i++)
  {
    const Point_2 current = p.vertex(i);
    const Segment_2 s(prev, current);

    VertexMapIterator sit = vertex_map.find(s.source());
    if (sit == vertex_map.end())
    {
      // Did not find vertex
      SegmentIterator si = pslg_segment_intersection(s.source(), vertices, segments);
      if (si != segments.end())
      {
        sit = pslg_split_edge(si,
                              s.source(),
                              vertex_map,
                              vertices,
                              segments);
      }
      else
      {
        std::pair<VertexMapIterator, bool> i = vertex_map.insert(std::make_pair(s.source(), vertex_map.size()));
        vertices.push_back(s.source());
        dolfin_assert(i.second);
        sit = i.first;
      }
    }

    VertexMapIterator tit = vertex_map.find(s.target());
    if (tit == vertex_map.end())
    {
      SegmentIterator si = pslg_segment_intersection(s.target(), vertices, segments);
      if (si != segments.end())
      {
        tit = pslg_split_edge(si,
                              s.target(),
                              vertex_map,
                              vertices,
                              segments);
      }
      else
      {
        std::pair<VertexMapIterator, bool> i = vertex_map.insert(std::make_pair(s.target(), vertex_map.size()));
        vertices.push_back(s.target());
        dolfin_assert(i.second);
        tit = i.first;
      }
    }

    if (segments.count(std::make_pair(tit->second, sit->second)) == 0)
      segments.insert(std::make_pair(sit->second, tit->second));

    prev = current;
  }
}
//-----------------------------------------------------------------------------
std::pair<std::vector<dolfin::Point>, std::vector<std::pair<std::size_t, std::size_t>>>
  CSGCGALDomain2D::compute_pslg(const std::vector<std::pair<std::size_t, CSGCGALDomain2D>>& domains)
{
  std::vector<Point_2> vertices;
  std::map<Point_2, std::size_t> vertex_map;
  std::set<std::pair<std::size_t, std::size_t>> segments;
  // std::vector<Segment_2> inserted_segments;

  for (const std::pair<std::size_t, CSGCGALDomain2D>& domain : domains)
  {
    const Polygon_set_2& p = domain.second.impl->polygon_set;

    std::list<Polygon_with_holes_2> polygon_list;
    p.polygons_with_holes(std::back_inserter(polygon_list));

    for (const Polygon_with_holes_2& pwh : polygon_list)
    {
      pslg_add_simple_polygon(vertex_map,
                              vertices,
                              segments,
                              pwh.outer_boundary());

      // Add holes
      Hole_const_iterator hit;
      for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit)
      {
        pslg_add_simple_polygon(vertex_map,
                                vertices,
                                segments,
                                *hit);
      }
    }
  }

  std::vector<dolfin::Point> v(vertices.size());
  for (const std::pair<Point_2, std::size_t>& vertex : vertex_map)
  {
    v[vertex.second] = dolfin::Point(CGAL::to_double(vertex.first.x()), CGAL::to_double(vertex.first.y()));
  }

  std::vector<std::pair<std::size_t, std::size_t>> s(segments.begin(), segments.end());

  {
    double shortest_segment = std::numeric_limits<double>::max();
    for (std::pair<std::size_t, std::size_t> segment : s)
    {
      shortest_segment = std::min(shortest_segment, (v[segment.first]-v[segment.second]).norm());
    }

    double closest_points = std::numeric_limits<double>::max();
    for (auto i = v.begin(); i != v.end(); i++)
    {
      for (auto j = i+1; j != v.end(); j++)
      {
        closest_points = std::min(closest_points, (*i-*j).norm());
      }
    }

    // std::cout << "Num vertices: " << v.size() << ", num edges: " << s.size() << std::endl;
    // std::cout << "Shortest edge: " << shortest_segment << ", closest: " << closest_points << std::endl;
  }

  return std::make_pair(std::move(v), std::move(s));
}
//-----------------------------------------------------------------------------
double CSGCGALDomain2D::shortest_edge() const
{
  FT shortest_edge = std::numeric_limits<double>::max();

  std::list<Polygon_with_holes_2> polygon_list;
  impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  // FIXME: For now it looks only at the outer boundary.
  // Should take holes into account as well
  for (const Polygon_with_holes_2& p : polygon_list)
  {
    const Polygon_2& bdr = p.outer_boundary();
    Point_2 prev = bdr.container().back();
    for (Polygon_2::Vertex_const_iterator vit = bdr.vertices_begin();
         vit != bdr.vertices_end(); vit++)
    {
      shortest_edge = std::min(shortest_edge, (*vit-prev).squared_length());
      prev = *vit;
    }
  }

  return sqrt(CGAL::to_double(shortest_edge));
}
//-----------------------------------------------------------------------------
void CSGCGALDomain2D::get_points_in_holes(std::vector<dolfin::Point>& holes) const
{
  std::list<Polygon_with_holes_2> polygon_list;
  impl->polygon_set.polygons_with_holes(std::back_inserter(polygon_list));

  for (const Polygon_with_holes_2& pwh : polygon_list)
  {
    Hole_const_iterator hit;
    for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit)
    {
      const Point_2 p = point_in_polygon(*hit);
      holes.push_back(dolfin::Point(CGAL::to_double(p.x()), CGAL::to_double(p.y())));
    }
  }
}
}
