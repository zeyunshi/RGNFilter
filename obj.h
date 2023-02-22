#ifndef OBJ_H
#define OBJ_H

//== INCLUDES ===================================================================

//STANDARD
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdio>
#include<string>
#include<vector>
#include<assert.h>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

#define REPORT_LOCATION \
  __FILE__ << " " << __LINE__ << " "

  
//=== IMPLEMENTATION ============================================================

//! @brief if this file line carries no information
inline bool is_line_invalid(const std::string &line) {
  return (line.empty() || 13 == line[0]);
}

//-------------------------------------------------------------------------------

template <typename OS, typename FLOAT_MATRIX, typename INT_MATRIX>
int obj2tri(OS &os, FLOAT_MATRIX &node, INT_MATRIX &tri) {
  if (!os) return -1;
  
  //get node and face number
  size_t node_num = 0, tri_num = 0;
  std::string line, word;
  while (getline(os, line)) {
    if (is_line_invalid(line)) continue;
    std::istringstream instream(line);
    instream >> word;

    if ("v" == word || "V" == word) {
      ++node_num;
    } else if ('f' == word[0] || 'F' == word[0]) {
      ++tri_num;
    }
    word.clear();
  }
  os.clear();
  os.seekg(0, std::ios::beg);

  //allocate space for triangle mesh
  if (0 == node_num || 0 == tri_num) {
    std::cerr << REPORT_LOCATION
              << "no vertex or face data find" << std::endl;
    return -1;
  }
  node.resize(3, node_num);
  tri.resize(3, tri_num);

  int process_state = 0;
  //get node and face data
  size_t node_cnt = 0, tri_cnt = 0;
  while (getline(os, line)) {
    if (is_line_invalid(line)) continue;
    std::istringstream instream(line);
    instream >> word;

    if ("v" == word || "V" == word) {
      double w = 1.0;
      instream >> node(0, node_cnt) >> node(1, node_cnt) >> node(2, node_cnt);
      if (instream >> w) {
        if (w < 1e-6) {
          std::cerr << REPORT_LOCATION
                    << "error occured when read vertex coordinates" << std::endl;
          process_state = -1;
          break;
        }
        for (size_t i = 0; i < 3; ++i)
          node(i, node_cnt) /= w;
      }
      ++node_cnt;
    } else if ('f' == word[0] || 'F' == word[0]) {
      //only triangle mesh supported
      std::string pair[3], test;
      instream >> pair[0] >> pair[1] >> pair[2] >> test;
      if (!test.empty() || !instream.eof()) {
        std::cerr << REPORT_LOCATION
                  << "only triangle mesh supported" << std::endl;
        process_state = -1;
        break;
      }

      //get vertex id in the this triangle face
      for (size_t i = 0; i < 3; ++i) {
        std::string vertex_str = pair[i].substr(0, pair[i].find('/'));

        long unsigned int tmp;
        sscanf(vertex_str.c_str(), "%lu", &tmp);
        tri(i, tri_cnt) = tmp;
        --tri(i, tri_cnt);
        if (tri(i, tri_cnt) >= node_num) {
          std::cerr << REPORT_LOCATION
                    << "vertex index in face exceed limit" << std::endl;
          process_state = -1;
          break;
        }
      }
      if (process_state != 0)
        break;
      ++tri_cnt;
    }
    word.clear();
  }
  // os.close();

  return process_state;
}

//-------------------------------------------------------------------------------

template <typename FLOAT_MATRIX, typename INT_MATRIX>
int load_obj(const char *file_name, FLOAT_MATRIX &pts, INT_MATRIX &tris) {
  std::ifstream os(file_name);
  if (!os) return __LINE__;

  return obj2tri(os, pts, tris);
}

//-------------------------------------------------------------------------------

template <typename OS, typename FLOAT_MATRIX>
int obj2point(OS &os, FLOAT_MATRIX &node, FLOAT_MATRIX &normal) {
  //get node number
  size_t node_num = 0;
  size_t normal_num = 0;
  std::string line, word;
  while (getline(os, line)) {
    if (is_line_invalid(line)) continue;
    std::istringstream instream(line);
    instream >> word;

    if ("v" == word || "V" == word)
      ++node_num;
    if ("vn" == word || "VN" == word)
      ++normal_num;
    word.clear();
  }
  os.clear();
  os.seekg(0, std::ios::beg);

  //allocate space for triangle mesh
  if (0 == node_num) {
    std::cerr << REPORT_LOCATION
              << "no vertex data find" << std::endl;
    return -1;
  }
  if (node_num != normal_num) {
    std::cerr << REPORT_LOCATION
              << "vertex number != normal number" << std::endl;
    return -1;
  }
  node.resize(3, node_num);
  normal.resize(3, node_num);

  int process_state = 0;
  //get node and face data
  size_t node_cnt = 0;
  size_t normal_cnt = 0;
  while (getline(os, line)) {
    if (is_line_invalid(line)) continue;
    std::istringstream instream(line);
    instream >> word;

    if ("v" == word || "V" == word) {
      double w = 1.0;
      instream >> node(0, node_cnt) >> node(1, node_cnt) >> node(2, node_cnt);
      if (instream >> w) {
        if (w < 1e-6) {
          std::cerr << REPORT_LOCATION
                    << "error occured when read vertex coordinates" << std::endl;
          process_state = -1;
          break;
        }
        for (size_t i = 0; i < 3; ++i)
          node(i, node_cnt) /= w;
      }
      ++node_cnt;
    } else if ("vn" == word || "VN" == word) {
      instream >> normal(0, normal_cnt)
               >> normal(1, normal_cnt) >> normal(2, normal_cnt);
      ++normal_cnt;
      
    }
    word.clear();
  }

  return process_state;
}

template <typename OS, typename FLOAT_MATRIX>
int obj2point(OS &os, FLOAT_MATRIX &node) {
  FLOAT_MATRIX normal;
  obj2point(os, node, normal);
}

//-------------------------------------------------------------------------------

template <typename OS, typename FLOAT, typename INT>
int tri2obj(OS &os,
             const FLOAT *node, size_t node_num,
             const INT *tri, size_t tri_num) {
  if (!os) return -1;
  for(size_t i = 0; i < node_num; ++i)
    os << "v " << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";
  for(size_t i = 0; i < tri_num; ++i)
    os << "f " << tri[i*3+0]+1 << " " << tri[i*3+1]+1 << " " << tri[i*3+2]+1 << "\n";
  return 0;
}

//-------------------------------------------------------------------------------

template <typename OS, typename FLOAT, typename INT>
int line2obj(OS &os,
             const FLOAT *node, size_t node_num,
             const INT *line, size_t line_num) {
  if (!os) return __LINE__;
  for(size_t i = 0; i < node_num; ++i)
    os << "v " << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";
  for(size_t i = 0; i < line_num; ++i)
    os << "l " << line[i*2+0]+1 << " " << line[i*2+1]+1 << "\n";
  return 0;
}

//-------------------------------------------------------------------------------

template <typename OS, typename FLOAT_MATRIX, typename INT_MATRIX>
int obj2line(OS &os, FLOAT_MATRIX &node, INT_MATRIX &lines) {
  if (!os) return -1;
  
  //get node and line number
  size_t node_num = 0, line_num = 0;
  std::string line, word;
  while (getline(os, line)) {
    if (is_line_invalid(line)) continue;
    std::istringstream instream(line);
    instream >> word;

    if ("v" == word || "V" == word) {
      ++node_num;
    } else if ('l' == word[0] || 'L' == word[0]) {
      ++line_num;
    }
    word.clear();
  }
  os.clear();
  os.seekg(0, std::ios::beg);

  //allocate space for line mesh
  if (0 == node_num || 0 == line_num) {
    std::cerr << REPORT_LOCATION
              << "no vertex or line data find" << std::endl;
    return -1;
  }
  node.resize(3, node_num);
  lines.resize(2, line_num);

  int process_state = 0;
  //get node and line data
  size_t node_cnt = 0, line_cnt = 0;
  while (getline(os, line)) {
    if (is_line_invalid(line)) continue;
    std::istringstream instream(line);
    instream >> word;

    if ("v" == word || "V" == word) {
      double w = 1.0;
      instream >> node(0, node_cnt) >> node(1, node_cnt) >> node(2, node_cnt);
      if (instream >> w) {
        if (w < 1e-6) {
          std::cerr << REPORT_LOCATION
                    << "error occured when read vertex coordinates" << std::endl;
          process_state = -1;
          break;
        }
        for (size_t i = 0; i < 3; ++i)
          node(i, node_cnt) /= w;
      }
      ++node_cnt;
    } else if ('l' == word[0] || 'L' == word[0]) {
      //only line mesh supported
      std::string pair[2], test;
      instream >> pair[0] >> pair[1] >> test;
      if (!test.empty() || !instream.eof()) {
        std::cerr << REPORT_LOCATION
                  << "only lines supported" << std::endl;
        process_state = -1;
        break;
      }

      //get vertex id in the this line 
      for (size_t i = 0; i < 2; ++i) {
        std::string vertex_str = pair[i].substr(0, pair[i].find('/'));

        long unsigned int tmp;
        sscanf(vertex_str.c_str(), "%lu", &tmp);
        lines(i, line_cnt) = tmp;
        --lines(i, line_cnt);
        if (lines(i, line_cnt) >= node_num) {
          std::cerr << REPORT_LOCATION
                    << "vertex index in line exceed limit" << std::endl;
          process_state = -1;
          break;
        }
      }
      if (process_state != 0)
        break;
      ++line_cnt;
    }
    word.clear();
  }
  // os.close();

  return process_state;
}


//-------------------------------------------------------------------------------

template <typename FLOAT_MATRIX, typename INT_MATRIX>
int save_obj(const char *file_name, FLOAT_MATRIX &pts, INT_MATRIX &tris) {
  std::ofstream os(file_name);
  if (!os) return __LINE__;

  return tri2obj(os, pts, tris);
}


//-------------------------------------------------------------------------------

template <typename OS, typename FLOAT, typename INT>
int tri2obj(OS &os,
            const FLOAT *node, const size_t node_num,
            const FLOAT *normal, const size_t normal_num,
            const INT *tri, const size_t tri_num) {
  for(size_t i = 0; i < node_num; ++i)
    os << "v " << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";
  for (size_t i = 0; i < normal_num; ++i)
    os << "vt " << normal[i*3+0] << " " << normal[i*3+1] << " " << normal[i*3+2] << "\n";
  if (node_num == normal_num) {
    for(size_t i = 0; i < tri_num; ++i)
      os << "f "
         << tri[i*3+0]+1 << "\\" << tri[i*3+0]+1 << " "
         << tri[i*3+1]+1 << "\\" << tri[i*3+1]+1 << " " 
         << tri[i*3+2]+1 << "\\" << tri[i*3+2]+1 << "\n";
  } else if (tri_num == normal_num) {
    for(size_t i = 0; i < tri_num; ++i)
      os << "f "
         << tri[i*3+0]+1 << "\\" << i << " "
         << tri[i*3+1]+1 << "\\" << i << " "
         << tri[i*3+2]+1 << "\\" << i << "\n";
  } else {
    std::cout << REPORT_LOCATION << "normal info error!" << std::endl;
  }
  return 0;
}

//-------------------------------------------------------------------------------

template <typename OS, typename FLOAT>
int point2obj(OS &os,
               const FLOAT *node, size_t node_num) {
  for(size_t i = 0; i < node_num; ++i)
    os << "v " << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";
  return 0;
}


//===============================================================================

#endif //OBJ_H

//===============================================================================
