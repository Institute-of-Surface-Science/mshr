// Copyright (C) 2012 Anders Logg, 2012-2017 Benjamin Kehlet
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


#ifndef __MSHR_PRIMITIVES_2D_H
#define __MSHR_PRIMITIVES_2D_H

#include "CSGPrimitive.h"

#include <dolfin/geometry/Point.h>
#include <vector>


namespace mshr
{

  /// @brief Base class for 2D primitives
  class CSGPrimitive2D : public CSGPrimitive
  {
  protected:
    CSGPrimitive2D();

  public:

    /// @brief get the dimension of the geometry (2)
    std::size_t dim() const { return 2; }
  };

  /// @brief A 2D circle
  ///
  /// { "small-icon" : "circle-small.png" }
  class Circle : public CSGPrimitive2D
  {
  public:

    /// @brief Create circle centered at x with radius r.
    ///
    /// @param c center.
    /// @param r radius.
    /// @param segments number of segments when computing the polygonal approximation
    Circle(dolfin::Point c, double r, std::size_t segments=0);

    /// @brief get informal string representation
    /// @param verbose  Verbosity level
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Circle; }

    /// @brief get center of circle
    dolfin::Point center() const { return c; }

    /// @brief get radius of circle
    double radius() const { return _r; }

    /// @brief get number of segments used when computing polygonal approximation
    std::size_t segments() const { return _segments; }

    std::pair<dolfin::Point, dolfin::Point> bounding_box() const;
    bool inside(dolfin::Point p) const;

  private:
    dolfin::Point c;
    double _r;
    const std::size_t _segments;
  };

  /// @brief A 2D ellipse
  ///
  /// { "small-icon" : "ellipse-small.png" }
  class Ellipse : public CSGPrimitive2D
  {
  public:

    /// @brief Create ellipse centered at c with horizontal semi-axis a and
    /// vertical semi-axis b.
    ///
    /// @param c        the center
    /// @param a        the horizontal semi-axis
    /// @param b        the vertical semi-axis
    /// @param segments the resolution when computing polygonal approximation
    Ellipse(dolfin::Point c, double a, double b, std::size_t segments=0);

    /// @brief get informal string representation
    /// @param verbose  Verbosity level
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Ellipse; }

    /// @brief get center of ellipse
    dolfin::Point center() const { return c; }

    /// @brief get horizontal semi-axis
    double a() const { return _a; }

    /// @brief get vertical semi-axis
    double b() const { return _b; }

    /// @brief get resolution when computing polygonal approximation
    std::size_t segments() const { return _segments; }

    std::pair<dolfin::Point, dolfin::Point> bounding_box() const;
    bool inside(dolfin::Point p) const;


  private:
    dolfin::Point c;
    double _a, _b;
    const std::size_t _segments;
  };

  /// @brief A 2D axis aligned rectangle
  ///
  /// { "small-icon" : "rectangle-small.png" }
  class Rectangle : public CSGPrimitive2D
  {
  public:

    /// @brief Create rectangle defined by two opposite corners
    ///
    /// @param a first corner.
    /// @param b second corner.
    Rectangle(dolfin::Point a, dolfin::Point b);

    /// @brief get informal string representation
    /// @param verbose  Verbosity level
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Rectangle; }

    /// @brief get first corner
    dolfin::Point first_corner() const { return a; }

    /// @brief get second corner
    dolfin::Point second_corner() const { return b; }

    std::pair<dolfin::Point, dolfin::Point> bounding_box() const;
    bool inside(dolfin::Point p) const;

  private:
    const dolfin::Point a, b;
  };

  /// @brief A 2D polygon
  ///
  /// { 'large-icon' : 'polygon-large.png', 'small-icon' : 'polygon-small.png' }
  class Polygon : public CSGPrimitive2D
  {
  public:

    /// @brief Create polygon defined by the given vertices. Vertices must be in counter-clockwise order and free of self-intersections.
    ///
    /// @param vertices A vector of dolfin::Points. Vertices are copied into the object.
    Polygon(const std::vector<dolfin::Point>& vertices);

    /// @brief get informal string representation
    /// @param verbose  Verbosity level
    std::string str(bool verbose) const;

    Type getType() const { return CSGGeometry::Polygon; }

    /// @brief check if vertices are counter clockwise oriented.
    bool ccw() const;

    /// @brief get vertices in polygon.
    const std::vector<dolfin::Point>& vertices() const { return _vertices; }

    std::pair<dolfin::Point, dolfin::Point> bounding_box() const;
    bool inside(dolfin::Point p) const;

  private:
    const std::vector<dolfin::Point> _vertices;
  };
}

#endif
