// Copyright (C) 2012 Anders Logg, 2014-2017 Benjamin Kehlet
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
//
// Modified by Johannes Ring, 2012

#include <mshr/CSGPrimitives2D.h>
#include <dolfin/math/basic.h>
#include <dolfin/log/LogStream.h>

#include <sstream>
#include <limits>
#include <algorithm>

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>

namespace mshr
{
//-----------------------------------------------------------------------------
CSGPrimitive2D::CSGPrimitive2D()
{}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Circle
//-----------------------------------------------------------------------------
Circle::Circle(dolfin::Point(c), double r, std::size_t segments)
  : c(c), _r(r), _segments(segments)
{
  if (_r < DOLFIN_EPS)
  {
    std::stringstream s;
    s << "Circle with center " << c.str() << " has zero or negative radius";
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                         "create circle",
                         s.str());
  }

  if (0 < _segments && _segments < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                         "create circle",
                         "Unable to create circle with fewer than 3 segments");
  }
}
//-----------------------------------------------------------------------------
std::string Circle::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Circle at (" << c.str() << ") with radius " << _r << ">";
  }
  else
  {
    s << "Circle(" << c.str() << ", " << _r << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> Circle::bounding_box() const
{
  return std::make_pair(dolfin::Point(c.x()-radius(), c.y()-radius()),
                        dolfin::Point(c.x()+radius(), c.y()+radius()));
}
//-----------------------------------------------------------------------------
bool Circle::inside(dolfin::Point p) const
{
  // dolfin::cout << "Inside: " << c << " r=" << _r << " ::: " << p << "(" << ((p-c).squared_norm() <= _r*_r ? "Inside" : "Outside") << dolfin::endl;
  return (p-c).squared_norm() <= _r*_r;
}


//-----------------------------------------------------------------------------
// Ellipse
//-----------------------------------------------------------------------------
Ellipse::Ellipse(dolfin::Point c, double a, double b,
                 std::size_t segments)
  : c(c), _a(a), _b(b), _segments(segments)
{
  if (_a < DOLFIN_EPS || _b < DOLFIN_EPS)
  {
    std::stringstream s;
    s << "Ellipse with center " << c.str() << " has invalid semi-axis";
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                         "create ellipse",
                         s.str());
  }

  if (0 < _segments && _segments < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create ellipse",
                 "Unable to create ellipse with fewer than 3 segments");
  }
}
//-----------------------------------------------------------------------------
std::string Ellipse::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Ellipse centered at (" << c.str() << ") with horizontal semi-axis "
      << _a << " and vertical semi-axis " << _b << ">";
  }
  else
  {
    s << "Ellipse(" << c.str() << ", " << _a << ", " << _b << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> Ellipse::bounding_box() const
{
  return std::make_pair(dolfin::Point(c.x()-a(), c.y()-b()),
                        dolfin::Point(c.x()+a(), c.y()+b()));
}
//-----------------------------------------------------------------------------
bool Ellipse::inside(dolfin::Point p) const
{
  const double x = (c.x()-p.x())/_a;
  const double y = (c.y()-p.y())/_b;
  return x*x + y*y <= 1;
}

//-----------------------------------------------------------------------------
// Rectangle
//-----------------------------------------------------------------------------
Rectangle::Rectangle(dolfin::Point a_, dolfin::Point b_)
  : a(dolfin::Point(std::min(a_.x(), b_.x()), std::min(a_.y(), b_.y()))),
    b(dolfin::Point(std::max(a_.x(), b_.x()), std::max(a_.y(), b_.y())))
{
  if (dolfin::near(a.x(), b.x()) || dolfin::near( a.y(), b.y()))
  {
    std::stringstream s;
    s << "Rectangle with corner " << a.str() << " and " << b.str() << " degenerated";
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create rectangle",
                         s.str());
  }
}
//-----------------------------------------------------------------------------
std::string Rectangle::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Rectangle with first corner at (" << a.str() << ") "
      << "and second corner at (" << b.str() << ")>";
  }
  else
  {
    s << "Rectangle( (" << a.str() << "), (" << b.str() << ") )";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> Rectangle::bounding_box() const
{
  return std::make_pair(a, b);
}
//-----------------------------------------------------------------------------
bool Rectangle::inside(dolfin::Point p) const
{
  return p.x() >= a.x() && p.x() <= b.x() && p.y() >= a.y() && p.y() <= b.y();
}


//-----------------------------------------------------------------------------
// Polygon
//-----------------------------------------------------------------------------
Polygon::Polygon(const std::vector<dolfin::Point>& vertices)
  : _vertices(vertices.begin(), vertices.end())
{
  if (_vertices.size() < 3)
  {
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create polygon",
                 "Polygon should have at least three vertices");
  }

  if (!ccw())
    dolfin::dolfin_error("CSGPrimitives2D.cpp",
                 "create polygon",
                 "Polygon vertices must be given in counter clockwise order");
}
//-----------------------------------------------------------------------------
std::string Polygon::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Polygon with vertices ";
    std::vector<dolfin::Point>::const_iterator p;
    for (p = _vertices.begin(); p != _vertices.end(); ++p)
    {
      s << "(" << p->x() << ", " << p->y() << ")";
      if ((p != _vertices.end()) && (p + 1 != _vertices.end()))
        s << ", ";
    }
    s << ">";
  }
  else
  {
    s << "Polygon (" << _vertices.size() << " vertices)";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
bool Polygon::ccw() const
{
  double signed_area = 0.0;

  dolfin::Point prev = _vertices.back();
  for (std::vector<dolfin::Point>::const_iterator it = _vertices.begin(),
	 v_end = _vertices.end();
       it != v_end;
       ++it)
  {
    signed_area += (prev.x()*it->y())-(it->x()*prev.y());
    prev = *it;
  }

  return signed_area > 0;
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> Polygon::bounding_box() const
{
  std::vector<dolfin::Point>::const_iterator it = _vertices.begin();
  dolfin::Point min = *it;
  dolfin::Point max = *it;
  it++;

  for (; it != _vertices.end(); it++)
  {
    const dolfin::Point& p = *it;
    min[0] = std::min(min[0], p[0]);
    min[1] = std::min(min[1], p[1]);
    min[2] = std::min(min[2], p[2]);

    max[0] = std::max(max[0], p[0]);
    max[1] = std::max(max[1], p[1]);
    max[2] = std::max(max[2], p[2]);
  }

  return std::make_pair(min, max);
}
//-----------------------------------------------------------------------------
bool Polygon::inside(dolfin::Point p) const
{
  typedef CGAL::Cartesian<double> K;
  typedef CGAL::Point_2<K> Point_2;
  typedef CGAL::Polygon_2<K> Polygon_2;

  Polygon_2 polygon;
  for (const dolfin::Point& vertex : _vertices)
    polygon.push_back(Point_2(vertex.x(), vertex.y()));

  return polygon.has_on_bounded_side(Point_2(p.x(), p.y()));
}

}
