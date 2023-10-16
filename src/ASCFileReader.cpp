// Copyright (C) 2016 Lars Magnus Valnes
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

#include <mshr/ASCFileReader.h>

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

namespace
{
template<typename T>
inline double convert_string(const std::string& s)
{
  std::istringstream is(s);
  T val;
  is >> val;

  return val;
}

// get next line of file and trim away whitespace
inline void get_next_line(std::ifstream& file, std::string& line, std::size_t &lineno)
{
  do
  {
    std::getline(file, line);
    boost::algorithm::trim(line);
    lineno++;
  } while ( !file.eof() && line == "");
}
} // end anonymous namespace
//-----------------------------------------------------------------------------
namespace mshr
{

void ASCFileReader::read(const std::string filename,
                         std::vector<std::array<double, 3> >& vertices,
                         std::vector<std::array<std::size_t, 3> >& facets)
{
  dolfin::log(dolfin:: TRACE, "Reading surface from %s ", filename.c_str());

  vertices.clear();
  facets.clear();
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  std::ifstream file(filename.c_str());
  if (!file.is_open())
  {
    dolfin::dolfin_error("ASCFileReader.cpp",
                         "open .asc file to read 3D surface",
                         "Failed to open file");
  }

  std::string line;
  std::size_t lineno = 0;
  const boost::char_separator<char> sep(" ");
  const boost::char_separator<char> sep2("  ");
  
  get_next_line(file, line, lineno);  
  // First line is not important
  get_next_line(file, line, lineno);
  tokenizer tokens(line, sep); 
  tokenizer::iterator tok_iter = tokens.begin();
  const std::size_t num_vertices = convert_string<std::size_t>(*tok_iter);
  tok_iter++;
  const std::size_t num_facets   = convert_string<std::size_t>(*tok_iter);
  tok_iter++;
  dolfin_assert(tok_iter == tokens.end());

  vertices.reserve(num_vertices);
  facets.reserve(num_facets);

  get_next_line(file, line, lineno);   
  for (std::size_t i = 0; i < num_vertices; i++)
  {
    tokenizer tokens(line, sep2);
    tokenizer::iterator tok_iter = tokens.begin();

    std::array<double, 3> vertex;
    vertex[0] = convert_string<double>(*tok_iter);
    tok_iter++;
    vertex[1] = convert_string<double>(*tok_iter);
    tok_iter++;
    vertex[2] = convert_string<double>(*tok_iter);
    tok_iter++;

    tok_iter++;

    dolfin_assert(tok_iter == tokens.end());
    vertices.push_back(vertex);
    get_next_line(file, line, lineno);
  }
  for (std::size_t i = 0; i < num_facets; i++)
  {
    tokenizer tokens(line, sep);
    tokenizer::iterator tok_iter = tokens.begin();

    std::array<std::size_t, 3> facet;
    facet[0] = convert_string<std::size_t>(*tok_iter);
    tok_iter++;
    facet[1] = convert_string<std::size_t>(*tok_iter);
    tok_iter++;
    facet[2] = convert_string<std::size_t>(*tok_iter);
    tok_iter++;
    tok_iter++;

    dolfin_assert(tok_iter == tokens.end());

    facets.push_back(facet);
    get_next_line(file, line, lineno); 

  }
  
}
//-----------------------------------------------------------------------------
void ASCFileReader::write(const std::string filename,
                          const std::vector<std::array<double, 3> >& vertices,
                          const std::vector<std::array<std::size_t, 3> >& facets)
{
  std::ofstream file(filename);
  file.precision(6);

  if (!file.is_open())
  {
    dolfin::dolfin_error("ASCFileReader.cpp",
                         "open file to write asc data",
                         "Failed to open file");
  }

  file << "#!" <<  std::endl;
  file << vertices.size() << " " << facets.size() << std::endl;
	// Orded vertices
  for (const std::array<double, 3>& v : vertices)
  {
    file << v[0] << "  " << v[1] << "  " << v[2] << "  " << 0 << std::endl;
  }

  for (const std::array<std::size_t, 3>& f : facets)
  {
    file << f[0] << " " << f[1] << " " << f[2] << " " << 0 <<  std::endl;
  }
}
}
