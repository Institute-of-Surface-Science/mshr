// Copyright (C) 2017 Benjamin Kehlet
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


#ifndef POLYGON_UTILS_H__
#define POLYGON_UTILS_H__

class PolygonUtils
{
 public:

  // Computes the orientation, assuming that it is well defined
  static bool ccw(const std::vector<dolfin::Point>& vertices)
  {
    double signed_area = 0.0;

    dolfin::Point prev = vertices.back();
    for (std::vector<dolfin::Point>::const_iterator it = vertices.begin(),
	   v_end = vertices.end();
	 it != v_end;
	 ++it)
    {
      signed_area += (prev.x()*it->y())-(it->x()*prev.y());
      prev = *it;
    }

    return signed_area > 0; 
  }
};


#endif
