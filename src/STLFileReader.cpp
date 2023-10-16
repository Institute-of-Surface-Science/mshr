// Copyright (C) 2014 Benjamin Kehlet
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

#include <mshr/STLFileReader.h>
#include "FuzzyPointLocator.h"

#include <dolfin/geometry/Point.h>
#include <dolfin/common/constants.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/log/log.h>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <string>
#include <map>
#include <math.h> 
namespace
{
inline double strToDouble(const std::string& s)
{
  std::istringstream is(s);
  double val;
  is >> val;

  return val;
}

// get next line of file and trim away whitespace
inline void get_next_line(std::ifstream& file, std::string& line, std::size_t &lineno)
{
  // Skip white
  do
  {
    std::getline(file, line);
    boost::algorithm::trim(line);
    lineno++;
  } while (file.good() && line.size() == 0);
}

} // end anonymous namespace
//-----------------------------------------------------------------------------
namespace mshr
{

void STLFileReader::read(const std::string filename,
                         std::vector<std::array<double, 3> >& vertices,
                         std::vector<std::array<std::size_t, 3> >& facets)
{

  dolfin::log(dolfin:: TRACE, "Reading surface from %s ", filename.c_str());

  // vertices.clear();
  facets.clear();

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  std::ifstream file(filename.c_str());
  if (!file.is_open())
  {
    dolfin::dolfin_error("STLFileReader.cpp",
                         "open .stl file to read 3D surface",
                         "Failed to open file");
  }

  FuzzyPointMap vertex_map(1e-10);

  std::string line;
  std::size_t lineno = 0;
  const boost::char_separator<char> sep(" ");

  // Read the first line and trim away whitespaces
  get_next_line(file, line, lineno);

  if (line.substr(0, 5) != "solid")
  {
    dolfin::dolfin_error("STLFileReader.cpp",
                         "open .stl file to read 3D surface",
                         "File does not start with \"solid\" (line %u", lineno);
  }

  // TODO: Read name of solid

  do
  {
    // Some files contain color information before the vertex information
    get_next_line(file, line, lineno);
  } while (line.substr(0, 5) != "facet");

  while (file.good())
  {
    bool has_normal = false;
    dolfin::Point normal;

    // Read the line "facet normal n1 n2 n3"
    {
      tokenizer tokens(line, sep);
      tokenizer::iterator tok_iter = tokens.begin();

      if (*tok_iter != "facet")
        dolfin::dolfin_error("STLFileReader.cpp",
                             "open .stl file to read 3D surface",
                             "Expected keyword \"facet\" (line %u)", lineno);
      ++tok_iter;

      // Check if a normal different from zero is given
      if (tok_iter != tokens.end())
      {
        if  (*tok_iter != "normal")
          dolfin::dolfin_error("STLFileReader.cpp",
                               "open .stl file to read 3D surface",
                               "Expected keyword \"normal\"(line %u)", lineno);
        ++tok_iter;

        //dolfin::cout << "Read line: " << line << dolfin::endl;

        for (std::size_t i = 0; i < 3; ++i)
        {
          normal[i] = strToDouble(*tok_iter);
          ++tok_iter;
        }

        if (normal.norm() > DOLFIN_EPS)
          has_normal = true;

        if (tok_iter != tokens.end())
          dolfin::dolfin_error("STLFileReader.cpp",
                               "open .stl file to read 3D surface",
                               "Expected end of line (line %u)", lineno);
      }
    }

    // if (has_normal)
    // {
    //   dolfin::cout << "Has normal" << dolfin::endl;
    //   dolfin::cout << normal << dolfin::endl;
    // }

    // Read "outer loop" line
    get_next_line(file, line, lineno);

    if (line != "outer loop")
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword 'outer loop' (line %u)", lineno);

    std::array<std::size_t, 3> v_indices;

    get_next_line(file, line, lineno);

    tokenizer tokens(line, sep);
    tokenizer::iterator tok_iter = tokens.begin();

    if (*tok_iter != "vertex")
    {
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword vertex (line %u)", lineno);
    }

    int counter = 0;

    // Read lines with vertices
    do
    {
      // Only support for triangulated surfaces for now
      dolfin_assert(counter < 3);

      // Advance to next
      ++tok_iter;

      const double x = strToDouble(*tok_iter); ++tok_iter;
      const double y = strToDouble(*tok_iter); ++tok_iter;
      const double z = strToDouble(*tok_iter); ++tok_iter;

      const std::array<double, 3> vertex = {{x, y, z}};

      // Insert point and get index. vertex_map takes care of merging close
      // vertices
      v_indices[counter] = vertex_map.insert_point(vertex);

      // Get next line
      get_next_line(file, line, lineno);

      tokens = tokenizer(line, sep);
      tok_iter = tokens.begin();

      counter++;
    } while (*tok_iter == "vertex");

    // Read 'endloop' line
    if (line != "endloop")
    {
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword endloop (line %u)", lineno);
    }

    get_next_line(file, line, lineno);
    if (line != "endfacet")
    {
      dolfin::dolfin_error("STLFileReader.cpp",
                           "open .stl file to read 3D surface",
                           "Expected keyword endfacet (line %u)", lineno);
    }

    // Add facet to output
    facets.push_back(v_indices);

    // Get orientation right if normal is given
    if (has_normal)
    {
      // Compute normal
      const dolfin::Point v1(3, vertex_map[v_indices[0]].data());
      const dolfin::Point v2(3, vertex_map[v_indices[1]].data());
      const dolfin::Point v3(3, vertex_map[v_indices[2]].data());

      const dolfin::Point a = v2-v1;
      const dolfin::Point b = v3-v1;

      dolfin::Point n = a.cross(b);
      n /= n.norm();

      // dolfin::cout << "Normal: " << n << dolfin::endl;
      // if ( (n - normal).norm() > 1e-5 )
      //   dolfin::cout << "Diff: " << (n - normal).norm() << dolfin::endl;
    }

    // Get next line
    // either start of next facet or endsolid
    get_next_line(file, line, lineno);

    if (line.substr(0, 5) != "facet")
      break;
  }

  // Read the 'endsolid' line
  tokenizer tokens(line, sep);
  tokenizer::iterator tok_iter = tokens.begin();

  if (*tok_iter != "endsolid")
  {
    dolfin::dolfin_error("STLFileReader.cpp",
                         "open .stl file to read 3D surface",
                         "Expected keyword endsolid at line %u", lineno);
  }
  ++tok_iter;

  vertices = vertex_map.get_points();

  // TODO: Check name of solid
  dolfin::log(dolfin::TRACE, "Done reading surface");

  // closest_vertices(vertex_map);
}
void STLFileReader::write(const std::string filename,
                         std::vector<std::array<double, 3> >& vertices,
                         std::vector<std::array<std::size_t, 3> >& facets)
{

    std::ofstream file(filename);
    file.precision(6);

    if (!file.is_open())
    {
       dolfin::dolfin_error("STLFileReader.cpp",
                         "open file to write stl data",
                         "Failed to open file");
    }

    file << "solid "<< filename << std::endl;
    

    
    
    for (const std::array<std::size_t, 3>& f : facets)
    {

       std::array<double,3> n;
       float vec1x,vec1y,vec1z;
       float vec2x,vec2y,vec2z;
       float size;
      
       std::array<double, 3>& v0 = vertices[f[0]];
       std::array<double, 3>& v1 = vertices[f[1]];
       std::array<double, 3>& v2 = vertices[f[2]];
       
       vec1x = v1[0] - v0[0];
       vec1y = v1[1] - v0[1];
       vec1z = v1[2] - v0[2];

       vec2x = v2[0] - v0[0];
       vec2y = v2[1] - v0[1];
       vec2z = v2[2] - v0[2];
       
       n[0] = vec1y*vec2z-vec1z*vec2y;
       n[1] = vec1z*vec2x-vec1x*vec2z;
       n[2] = vec1x*vec2y-vec1y*vec2x; 		
       size = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

       file << "facet normal " << n[0]/size << " " << n[1]/size  << " " << n[2]/size <<std::endl;
       file << "outer loop"<< std::endl;
       file << "\t" << "vertex " << v0[0] << " " << v0[1]<< " " << v0[2] << std::endl;
       file << "\t" << "vertex " << v1[0] << " " << v1[1]<< " " << v1[2] << std::endl;
       file << "\t" << "vertex " << v2[0] << " " << v2[1]<< " " << v2[2] << std::endl;
       file <<"endloop" << std::endl;
       file <<"endfacet"<< std::endl;
    }
    file << "endsolid"  << std::endl;

}

}
