// Copyright (C) 2012 Anders Logg, Benjamin Kehlet 2013-2017
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

#include <mshr/CSGOperators.h>

#include <dolfin/common/utils.h>
#include <dolfin/log/log.h>
#include <dolfin/common/constants.h>

#include <sstream>

namespace mshr
{

//-----------------------------------------------------------------------------
// CSGUnion
//-----------------------------------------------------------------------------
CSGOperator::CSGOperator()
{}

//-----------------------------------------------------------------------------
// CSGUnion
//-----------------------------------------------------------------------------
CSGUnion::CSGUnion(std::shared_ptr<CSGGeometry> g0,
                   std::shared_ptr<CSGGeometry> g1)
  : _g0(g0), _g1(g1)
{
  assert(g0);
  assert(g1);

  // Check dimensions
  if (g0->dim() != g1->dim())
  {
    dolfin::dolfin_error("CSGOperators.cpp",
                         "create union of CSG geometries",
                         "Dimensions of geometries don't match (%d vs %d)",
                         g0->dim(), g1->dim());
  }

  dim_ = g0->dim();
}
//-----------------------------------------------------------------------------
std::string CSGUnion::str(bool verbose) const
{
  assert(_g0);
  assert(_g1);

  std::stringstream s;

  if (verbose)
  {
    s << "<Union>\n"
      << "{\n"
      << dolfin::indent(_g0->str(true))
      << "\n"
      << dolfin::indent(_g1->str(true))
      << "\n}";
  }
  else
  {
    s << "(" << _g0->str(false) << " + " << _g1->str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> CSGUnion::bounding_box() const
{
  std::pair<dolfin::Point, dolfin::Point> a = _g0->bounding_box();
  std::pair<dolfin::Point, dolfin::Point> b = _g1->bounding_box();

  return std::make_pair(dolfin::Point(std::min(a.first.x(), b.first.x()),
                                      std::min(a.first.y(), b.first.y()),
                                      std::min(a.first.z(), b.first.z())),
                        dolfin::Point(std::max(a.second.x(), b.second.x()),
                                      std::max(a.second.y(), b.second.y()),
                                      std::max(a.second.z(), b.second.z())));
}
//-----------------------------------------------------------------------------
bool CSGUnion::inside(dolfin::Point p) const
{
  return _g0->inside(p) || _g1->inside(p);
}

//-----------------------------------------------------------------------------
// CSGDifference
//-----------------------------------------------------------------------------
CSGDifference::CSGDifference(std::shared_ptr<CSGGeometry> g0,
			     std::shared_ptr<CSGGeometry> g1)
  : _g0(g0), _g1(g1)
{
  assert(g0);
  assert(g1);

  // Check dimensions
  if (g0->dim() != g1->dim())
  {
    dolfin::dolfin_error("CSGOperators.cpp",
                         "create difference of CSG geometries",
                         "Dimensions of geomestries don't match (%d vs %d)",
                         g0->dim(), g1->dim());
  }

  dim_ = g0->dim();
}
//-----------------------------------------------------------------------------
std::string CSGDifference::str(bool verbose) const
{
  assert(_g0);
  assert(_g1);

  std::stringstream s;

  if (verbose)
  {
    s << "<Difference>\n"
      << "{\n"
      << dolfin::indent(_g0->str(true))
      << "\n"
      << dolfin::indent(_g1->str(true))
      << "\n}";
  }
  else
  {
    s << "(" << _g0->str(false) << " - " << _g1->str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> CSGDifference::bounding_box() const
{
  // FIXME: This is where the bounding_box implementation may overshoot
  return _g0->bounding_box();
}
//-----------------------------------------------------------------------------
bool CSGDifference::inside(dolfin::Point p) const
{
  return _g0->inside(p) && !_g1->inside(p);
}

//-----------------------------------------------------------------------------
// CSGIntersection
//-----------------------------------------------------------------------------
CSGIntersection::CSGIntersection(std::shared_ptr<CSGGeometry> g0,
                                 std::shared_ptr<CSGGeometry> g1)
  : _g0(g0), _g1(g1)
{
  assert(g0);
  assert(g1);

  // Check dimensions
  if (g0->dim() != g1->dim())
  {
    dolfin::dolfin_error("CSGOperators.cpp",
                         "create intersection of CSG geometries",
                         "Dimensions of geomestries don't match (%d vs %d)",
                         g0->dim(), g1->dim());
  }

  dim_ = g0->dim();
}
//-----------------------------------------------------------------------------
std::string CSGIntersection::str(bool verbose) const
{
  assert(_g0);
  assert(_g1);

  std::stringstream s;

  if (verbose)
  {
    s << "<Intersection>\n"
      << "{\n"
      << dolfin::indent(_g0->str(true))
      << "\n"
      << dolfin::indent(_g1->str(true))
      << "\n}";
  }
  else
  {
    s << "(" << _g0->str(false) << " * " << _g1->str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> CSGIntersection::bounding_box() const
{
  std::pair<dolfin::Point, dolfin::Point> a = _g0->bounding_box();
  std::pair<dolfin::Point, dolfin::Point> b = _g1->bounding_box();

  return std::make_pair(dolfin::Point(std::min(a.first.x(), b.first.x()),
                                      std::min(a.first.y(), b.first.y()),
                                      std::min(a.first.z(), b.first.z())),
                        dolfin::Point(std::max(a.second.x(), b.second.x()),
                                      std::max(a.second.y(), b.second.y()),
                                      std::max(a.second.z(), b.second.z())));
}
//-----------------------------------------------------------------------------
bool CSGIntersection::inside(dolfin::Point p) const
{
  return _g0->inside(p) && _g1->inside(p);
}

//-----------------------------------------------------------------------------
// CSGTranslation
//-----------------------------------------------------------------------------
CSGTranslation::CSGTranslation(std::shared_ptr<CSGGeometry> g,
                               dolfin::Point t)
  : g(g), t(t)
{
  assert(g);

  dim_ = g->dim();
}
//-----------------------------------------------------------------------------
std::string CSGTranslation::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    s << "<Translation>\n"
      << "{\n"
      << dolfin::indent(g->str(true) + "\nby\n" + t.str(true))
      << "\n}";
  }
  else
  {
    s << "(" << g->str(false) << " + " << t.str(false) << ")";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> CSGTranslation::bounding_box() const
{
  std::pair<dolfin::Point, dolfin::Point> aabb = g->bounding_box();

  return std::make_pair(aabb.first+t, aabb.second+t);
}
//-----------------------------------------------------------------------------
bool CSGTranslation::inside(dolfin::Point p) const
{
  return g->inside(p-t);
}

//-----------------------------------------------------------------------------
// CSGScaling
//-----------------------------------------------------------------------------
CSGScaling::CSGScaling(std::shared_ptr<CSGGeometry> g,
                       dolfin::Point c,
                       double s)
  : g(g), c(c), s(s), translate(true)
{
  assert(g);

  dim_ = g->dim();
}
//-----------------------------------------------------------------------------
CSGScaling::CSGScaling(std::shared_ptr<CSGGeometry> g,
                       double s)
  : g(g), c(0,0,0), s(s), translate(false)
{
  assert(g);

  dim_ = g->dim();
}
//-----------------------------------------------------------------------------
std::string CSGScaling::str(bool verbose) const
{
  std::stringstream ss;

  if (verbose)
  {
    ss << "<Scaling>\n"
      << "{\n"
      << dolfin::indent(g->str(true) + "\nby\n" + std::to_string(s));

      if (translate)
        ss << "\naround " << c.str(true);

      ss << "\n}";
  }
  else
  {
    ss << "(" << g->str(false) << " * " << std::to_string(s);
    if (translate)
      ss << "(" << c.str(true) << ")";
    ss << ")";
  }

  return ss.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> CSGScaling::bounding_box() const
{
  std::pair<dolfin::Point, dolfin::Point> aabb = g->bounding_box();

  const dolfin::Point scaled_c = (translate ? (1-s)*c : dolfin::Point(0,0,0));

  if (translate)
    return std::make_pair((aabb.first-c)*s + c, (aabb.second-c)*s + c);
  else
    return std::make_pair(aabb.first*s, aabb.second*s);

}
//-----------------------------------------------------------------------------
bool CSGScaling::inside(dolfin::Point p) const
{
  const dolfin::Point p_scaled = (translate ? (p-c)*s + c : p*s);
  return g->inside(p_scaled);
}

//-----------------------------------------------------------------------------
// CSGRotation
//-----------------------------------------------------------------------------
CSGRotation::CSGRotation(std::shared_ptr<CSGGeometry> g,
                         double theta)
  : g(g), rot_axis(.0,.0), c(.0,.0), theta(theta), translate(false)
{

  dim_ = g->dim();

  if (dim_ > 2)
    dolfin::dolfin_error("CSGOperators.cpp",
                         "Constructing CSG rotation",
                         "Rotation axis must be given in 3D");
}
//-----------------------------------------------------------------------------
CSGRotation::CSGRotation(std::shared_ptr<CSGGeometry> g,
                         dolfin::Point v,
                         double theta)
  : g(g),
    rot_axis(v),
    c(v),
    theta(theta),
    translate(g->dim() == 2 ? true : false)
{
  assert(g);

  dim_ = g->dim();

  // if (dim_ == 2)
  //   translate = true;
  // else
  //   translate = false;
}
//-----------------------------------------------------------------------------
CSGRotation::CSGRotation(std::shared_ptr<CSGGeometry> g,
                         dolfin::Point rot_axis,
                         dolfin::Point rot_center,
                         double theta)
  : g(g),
    rot_axis(rot_axis),
    c(rot_center),
    theta(theta),
    translate(true)
{
  assert(g);

  dim_ = g->dim();

  if (dim_ < 3)
    dolfin::dolfin_error("CSGOperators.cpp",
                         "Constructing CSG rotation",
                         "Can't give rotation axis for dimension < 3");
}
//-----------------------------------------------------------------------------
std::string CSGRotation::str(bool verbose) const
{
  std::stringstream ss;

  if (verbose)
  {
    ss << "<Rotation>\n"
      << "{\n"
      << dolfin::indent(g->str(true)
                        + (translate ? "\naround "+rot_axis.str(true) : "")
                        + "\nby " + std::to_string(theta/DOLFIN_PI) + " PI");

      ss << "\n}";
  }
  else
  {
    ss << "rotate(" << g->str(false)
       << ", " << std::to_string(theta/DOLFIN_PI) << " PI";

    if (translate)
      ss << ", " << rot_axis.str(true);

    ss << ")";
  }

  return ss.str();
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, dolfin::Point> CSGRotation::bounding_box() const
{
  if (dim() != 2)
  {
    dolfin_not_implemented();
  }

  const std::pair<dolfin::Point, dolfin::Point> aabb = g->bounding_box();
  const dolfin::Point axis(0,0,1);

  // Rotate all corners
  const dolfin::Point a = aabb.first.rotate(axis, theta);
  const dolfin::Point b = dolfin::Point(aabb.first.x(), aabb.second.y()).rotate(axis, theta);
  const dolfin::Point c = dolfin::Point(aabb.second.x(), aabb.first.y()).rotate(axis, theta);
  const dolfin::Point d = aabb.second.rotate(axis, theta);

  return std::make_pair(dolfin::Point(std::min(a.x(), std::min(b.x(), std::min(c.x(), d.x()))),
                                      std::min(a.y(), std::min(b.y(), std::min(c.y(), d.y())))),
                        dolfin::Point(std::max(a.x(), std::max(b.x(), std::max(c.x(), d.x()))),
                                      std::max(a.y(), std::max(b.y(), std::max(c.y(), d.y())))));

}
//-----------------------------------------------------------------------------
bool CSGRotation::inside(dolfin::Point p) const
{
  const dolfin::Point p_rotated = translate ? (p-c).rotate(dolfin::Point(0,0,1), -theta) + c  : p.rotate(dolfin::Point(0,0,1), -theta);
  return g->inside(p_rotated);
}

}
