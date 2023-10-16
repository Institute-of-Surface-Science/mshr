// Copyright (C) 2012 Anders Logg and 2013-2017 Benjamin Kehlet
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


#ifndef __MSHR_GEOMETRY_H
#define __MSHR_GEOMETRY_H

#include <memory>
#include <cstddef>
#include <vector>
#include <list>

#include <dolfin/common/Variable.h>
#include <dolfin/geometry/Point.h>

namespace mshr
{

  /// Geometry described by Constructive Solid Geometry (CSG)
  class CSGGeometry : public dolfin::Variable
  {
   protected:
    /// Constructor
    CSGGeometry();

   public:
    /// Destructor
    virtual ~CSGGeometry();

    /// Return dimension of geometry
    virtual std::size_t dim() const = 0;

    /// Informal string representation
    virtual std::string str(bool verbose) const = 0;

    /// Define subdomain. This feature is 2D only.
    /// The subdomain is itself a CSGGeometry and the corresponding
    /// cells in the resulting will be marked with i
    /// If subdomains overlap, the latest added will take precedence.
    void set_subdomain(std::size_t i, std::shared_ptr<CSGGeometry> s);

    /// Define subdomain. This feature is 2D only.
    /// The subdomain is itself a CSGGeometry and the corresponding
    /// cells in the resulting will be marked with i
    /// If subdomains overlap, the latest added will take precedence.
    void set_subdomain(std::size_t i, CSGGeometry& s);

    /// @brief Has subdomains been set
    bool has_subdomains() const;

    /// @brief Return const list of subdomain geometries
    const std::list<std::pair<std::size_t, std::shared_ptr<const CSGGeometry>>>& get_subdomains() const { return subdomains; }

    // These functions are (for now) implemented for 2D only.
    std::pair<dolfin::Point, double> estimate_bounding_sphere(std::size_t numSamples=100) const;

    // Compute axis aligned box that is guaranteed to bound the geometry. May overshoot.
    virtual std::pair<dolfin::Point, dolfin::Point> bounding_box() const = 0;

    // Test if given point is inside geometry. 2D only (for now)
    virtual bool inside(dolfin::Point p) const = 0;

    // Test if line segment is entirely inside geometry. 2D only (for now)
    virtual bool inside(dolfin::Point p1, dolfin::Point p2) const;

    enum Type { Box,
                Sphere,
                Cylinder,
                Tetrahedron,
                Ellipsoid,
                Surface3D,
                Extrude2D,
                Circle,
                Ellipse,
                Rectangle,
                Polygon,
                Union,
                Intersection,
                Difference,
                Translation,
                Scaling,
                Rotation,
                TriPolyhedron };

    virtual Type getType() const = 0;
    virtual bool is_operator() const = 0;

   private :
    std::list<std::pair<std::size_t, std::shared_ptr<const CSGGeometry> > > subdomains;
  };
}

#endif
