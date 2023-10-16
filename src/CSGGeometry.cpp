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

#include <dolfin/common/NoDeleter.h>
#include <dolfin/log/LogStream.h>
#include <mshr/CSGGeometry.h>

// Bounding sphere computation
#include <CGAL/Cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>

namespace mshr
{

//-----------------------------------------------------------------------------
CSGGeometry::CSGGeometry()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
CSGGeometry::~CSGGeometry()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void CSGGeometry::set_subdomain(std::size_t i, std::shared_ptr<CSGGeometry> s)
{
  if (dim() != 2)
    dolfin::dolfin_error("CSGGeometry.cpp",
                 "setting subdomain",
                 "Subdomains are currently supported only in 2D");

  if (s->dim() != dim())
    dolfin::dolfin_error("CSGGeometry.cpp",
                 "setting subdomain",
                 "Subdomain and domain must be of same dimension. Domain was dimension %d and subdomain was %d", dim(), s->dim());

  if (i == 0)
  {
    dolfin::dolfin_error("CSGGeometry.cpp",
                         "Setting reserved CSG subdomain (0)",
                         " Subdomain 0 is reserved and cannot be set by user");
  }

  // Check if i already used
  std::list<std::pair<std::size_t, std::shared_ptr<const CSGGeometry> > >::iterator it = subdomains.begin();
  while (it != subdomains.end())
  {
    if (it->first == i)
    {
      dolfin::warning("Double declaration of CSG subdomain with index %u.", i);

       // Remove existing declaration
       it = subdomains.erase(it);
    }
    else
      ++it;
  }

  subdomains.push_back(std::make_pair(i, s));
}
//-----------------------------------------------------------------------------
void CSGGeometry::set_subdomain(std::size_t i, CSGGeometry& s)
{
  set_subdomain(i, reference_to_no_delete_pointer(s));
}
//-----------------------------------------------------------------------------
bool CSGGeometry::has_subdomains() const
{
  return subdomains.size() > 0;
}
//-----------------------------------------------------------------------------
bool CSGGeometry::inside(dolfin::Point p1, dolfin::Point p2) const
{
  dolfin_not_implemented();
  return false;
}
//-----------------------------------------------------------------------------
std::pair<dolfin::Point, double> CSGGeometry::estimate_bounding_sphere(std::size_t numSamples) const
{
  std::pair<dolfin::Point, dolfin::Point> aabb = bounding_box();
  const dolfin::Point min = aabb.first;
  const dolfin::Point max = aabb.second;

  //typedef CGAL::Exact_predicates_exact_constructions_kernel K;
  typedef CGAL::Cartesian<double> K;
  typedef K::Point_2                      Point_2;
  typedef CGAL::Min_sphere_of_spheres_d_traits_2<K, K::FT> MinSphereTraits;
  typedef CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
  typedef MinSphereTraits::Sphere Sphere;

  std::vector<Sphere> S;
  S.reserve(numSamples);
  const std::size_t N = static_cast<std::size_t>(std::sqrt(static_cast<double>(numSamples)+.5));
  const double dX = (max.x()-min.x())/N;
  const double dY = (max.y()-min.y())/N;
  // const double dTheta = 2*DOLFIN_PI/N;
  for (std::size_t i = 0; i < N; i++)
  {
    for (std::size_t j = 0; j < N; j++)
    {
      const dolfin::Point p(min.x() + i*dX, min.y() + j*dY);
      if (inside(p))
        S.push_back(Sphere(Point_2(p.x(), p.y()), 0.0));
    }
  }

  Min_sphere ms(S.begin(), S.end());
  CGAL_assertion(ms.is_valid());

  auto it = ms.center_cartesian_begin();
  const double center_x = CGAL::to_double(*it);
  it++;
  const double center_y = CGAL::to_double(*it);
  const double radius = CGAL::to_double(ms.radius());

  return std::make_pair(dolfin::Point(center_x, center_y),
                        radius);
}

}
