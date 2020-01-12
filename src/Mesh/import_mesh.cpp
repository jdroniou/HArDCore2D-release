// Class to read a typ2 mesh file
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include <import_mesh.hpp>
#include <iterator>
#include <iostream>
#include <ctype.h>  //ispace
#include <stdio.h>

using namespace HArDCore2D;
MeshReaderTyp2::MeshReaderTyp2(std::string file_name)
    : _file_name(file_name) {}

bool MeshReaderTyp2::read_mesh(std::vector<std::vector<double> >& vertices,
                                 std::vector<std::vector<size_t> >& cells,
                                 std::vector<std::vector<double> >& centers) {
    std::ifstream in_file(
        _file_name.c_str());  // check here to see if file name has extension?
    if (in_file.good() == false) {
        return false;
    }
    // this assumes the file is formated
    //"Vertices"
    // int nv
    // double double
    //... ...
    // cells
    // int nc
    // double double
    //.... ....
//    std::cout << "Reading Typ2 Mesh" << std::endl;
    std::string string;
    std::string vertices_name("Vertices");
    std::string cells_name("cells");
    std::string centers_name("centers");
    // flag to know whether importing vert or cells
    Flag flag;
    flag = header_;
    size_t nV = 0;
    size_t nC = 0;
    while (in_file.good()) {
        std::getline(in_file, string);
        if (string.empty()) {
            std::cout << string << std::endl;
            continue;
        }
        // clean the file
        while (std::isspace(string.at(0))) {
            string.erase(0, 1);
        }  // str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
        while (std::isspace(*(string.end() - 1))) {
            string.erase(string.end() - 1);
        }
        if (string.compare(vertices_name) == 0) {
            flag = nv_;
            continue;
        }
        if (string.compare(cells_name) == 0) {
            flag = nc_;
            continue;
        }
        if (string.compare(centers_name) == 0) {
            flag = centers_;
            continue;
        }
        std::istringstream buffer(string);

        switch (flag) {
            case nv_: {
                nV = std::stoi(string);
                flag = vertices_;
                continue;
            }
            case vertices_: {
                std::vector<double> line(
                    (std::istream_iterator<double>(buffer)),
                    std::istream_iterator<double>());

                std::vector<double> tmp;
                for (std::vector<double>::iterator it = line.begin(); it != line.end(); it++) {
                    tmp.push_back(*it);
                }
                vertices.push_back(tmp);
                continue;
            }
            case nc_: {
                nC = std::stoi(string);
                flag = cells_;
                continue;
            }
            case cells_: {
                std::vector<size_t> line((std::istream_iterator<size_t>(buffer)),
                                      std::istream_iterator<size_t>());
                std::vector<size_t> tmp;
                for (std::vector<size_t>::iterator it = line.begin(); it != line.end(); it++) {
                    tmp.push_back(*it - 1);
                }
                cells.push_back(tmp);  // adjust from matlab/fotran indexes to c++
                continue;
            }
            case centers_: {
                std::vector<double> line(
                    (std::istream_iterator<double>(buffer)),
                    std::istream_iterator<double>());

                std::vector<double> tmp;
                for (std::vector<double>::iterator it = line.begin(); it != line.end(); it++) {
                    tmp.push_back(*it);
                }
                centers.push_back(tmp);
                continue;
            }
      default:
    break;
        }
        in_file.close();
    }
    size_t found = _file_name.find_last_of("/\\");
    std::string file_name_noroot = _file_name.substr(found+1);
    printf("%s file read: %lu/%lu cells, %lu/%lu vertices\n",file_name_noroot.c_str(),cells.size(),nC,vertices.size(),nV);
    return true;
}
