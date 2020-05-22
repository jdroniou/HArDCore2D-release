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

     size_t num_objects = 4;
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

     std::string l2rate;
     std::string h1rate;
     std::string energyrate;

     const int precision = 4;
     const int width = 12;

     std::cout.precision(precision);

     std::cout << "\n ---- Compute convergence rates ----\n";
     std::cout << "L2 error:\n";
     std::cout << std::setw(width) << "Mesh size";
     std::cout << std::setw(width) << "Error";
     std::cout << std::setw(width) << "Rate";
     std::cout << "\n";
     for (int i = 0; i < num_items; i++)
     {
          if (i >= 1)
          {
               l2rate = std::to_string((log(l2errors[i] / l2errors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << l2errors[i];
          std::cout << std::setw(width) << l2rate;
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
               h1rate = std::to_string((log(h1errors[i] / h1errors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
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
               energyrate = std::to_string((log(energyerrors[i] / energyerrors[i - 1])) / (log(meshsizes[i] / meshsizes[i - 1])));
          }
          std::cout << std::setw(width) << meshsizes[i];
          std::cout << std::setw(width) << energyerrors[i];
          std::cout << std::setw(width) << energyrate;
          std::cout << "\n";
     }

     return 0;
}