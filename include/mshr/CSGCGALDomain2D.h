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

#ifndef __MSHR_CSGCGAL_DOMAIN2D_H
#define __MSHR_CSGCGAL_DOMAIN2D_H

#include <mshr/CSGGeometry.h>

#include <dolfin/geometry/Point.h>
#include <memory>


namespace mshr
{
// Forward declaration
struct CSGCGALDomain2DImpl;

/// @brief Polygonal meshing domain
///
/// Represents the polygonal meshing domain to be sent to the mesh
/// generator. It uses CGAL's Polygon_2 as backend and utilizied CGAL's
/// 2D Regularized Boolean Set-Operations package.
class CSGCGALDomain2D : public dolfin::Variable
{
 public:
  /// @brief Create empty polygon
  CSGCGALDomain2D();

  /// @brief Construct polygon from Dolfin CSG geometry
  CSGCGALDomain2D(
      std::shared_ptr<const CSGGeometry> geometry,
      double segment_granularity
      );

  /// @brief Destructor
  ~CSGCGALDomain2D();

  /// @brief Copy constructor
  CSGCGALDomain2D(const CSGCGALDomain2D &other);

  /// @brief Assignment operator
  CSGCGALDomain2D &operator=(const CSGCGALDomain2D &other);

  /// @brief Boolean join
  void join_inplace(const CSGCGALDomain2D& other);

  /// @brief Boolean intersection
  void intersect_inplace(const CSGCGALDomain2D& other,
                         double truncate_tolerance=-1.);

  /// @brief Boolean difference
  void difference_inplace(const CSGCGALDomain2D& other);

  /// @brief Determine if point is inside the polygon.
  bool point_in_domain(dolfin::Point p) const;

  /// @brief Compute the radius of the minimum bounding circle.
  double compute_boundingcircle_radius() const;

  std::size_t num_polygons() const;
  double shortest_edge() const;

  std::vector<dolfin::Point> get_outer_polygon(std::size_t i) const;

  /// @brief get one point per hole, strictly inside the hole.
  /// @param holes the returned points
  void get_points_in_holes(std::vector<dolfin::Point>& holes) const;

  /// @brief Informal string representation
  /// @param verbose Verbosity level
  std::string str(bool verbose) const;

  static
    std::pair<std::vector<dolfin::Point>,
              std::vector<std::pair<std::size_t, std::size_t>>>
    compute_pslg(const std::vector<std::pair<std::size_t, CSGCGALDomain2D>>& domains);

 private:
  std::unique_ptr<CSGCGALDomain2DImpl> impl;
};

}

#endif
