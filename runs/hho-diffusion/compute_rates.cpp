#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <bits/stdc++.h>

int main()
{

     std::ifstream inFile;

     inFile.open("outputs/data_rates.dat");
     if (!inFile)
     {
          std::cerr << "Unable to open file\n";
          exit(1);
     }

     size_t num_objects = 7;
     std::string item;
     std::vector<std::string> all_items;
     std::vector<double> each_items[num_objects];

     while (inFile >> item) 
     {    
          all_items.push_back(item);
     }

     inFile.close();
     size_t num_items = all_items.size();
     for (size_t i = num_objects; i < num_items; i++) {
          double temp = std::stod(all_items[i]);
          for (size_t j = 0; j < num_objects; j++) {
               if(i % num_objects == j) {
                    each_items[j].push_back(temp);
                    break;
               }
          } 
     }

     num_items /= num_objects;
     num_items--;

     std::vector<double> meshsizes = each_items[0];
     std::vector<double> l2errors = each_items[1];
     std::vector<double> h1errors = each_items[2];
     std::vector<double> energyerrors = each_items[3];
     std::vector<double> nbedgedofs = each_items[4];
     std::vector<double> meshreg = each_items[5];
     std::vector<double> meshskewness = each_items[6];

     char l2rate[10];
     char h1rate[10];
     char energyrate[10];

     const int precision = 2;
     const int width = 12;

     std::cout.precision(precision);
     std::cout << std::scientific;

     std::cout << "\n ---- Compute convergence rates ----\n";
     std::cout << "L2 error:\n";
     std::cout << std::setw(width) << "Mesh size";
     std::cout << std::setw(width) << "Error";
     std::cout << std::setw(width) << "Rate";
     std::cout << std::setw(width) << "\t(mesh reg / skew)";
     std::cout << "\n";
     for (int i = 0; i < num_items; i++)
     {
          if (i >= 1)
          {
               sprintf(l2rate, "%5.3f",(log(l2errors[i] / l2errors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }else{
               sprintf(l2rate, " ");
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << l2errors[i];
          std::cout << std::setw(width) << l2rate;
          char meshregskew[100];
          sprintf(meshregskew, "\t(%5.2e / %5.2e)", meshreg[i], meshskewness[i]);
          std::cout << std::setw(width) << meshregskew;
          std::cout << "\n";
     }
     std::cout << "----------------------------------------\n";
     std::cout << "H1 error:\n";
     std::cout << std::setw(width) << "Mesh size";
     std::cout << std::setw(width) << "Error";
     std::cout << std::setw(width) << "Rate";
     std::cout << "\n";
     for (int i = 0; i < num_items; i++)
     {
          if (i >= 1)
          {
               sprintf(h1rate, "%5.3f",(log(h1errors[i] / h1errors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }else{
               sprintf(h1rate, " ");
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << h1errors[i];
          std::cout << std::setw(width) << h1rate;
          std::cout << "\n";
     }
     std::cout << "----------------------------------------\n";
     std::cout << "Energy error:\n";
     std::cout << std::setw(width) << "Mesh size";
     std::cout << std::setw(width) << "Error";
     std::cout << std::setw(width) << "Rate";
     std::cout << "\n";
     for (int i = 0; i < num_items; i++)
     {
          if (i >= 1)
          {
               sprintf(energyrate, "%5.3f",(log(energyerrors[i] / energyerrors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }else{
               sprintf(energyrate, " ");
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << energyerrors[i];
          std::cout << std::setw(width) << energyrate;
          std::cout << "\n";
     }

     return 0;
}
