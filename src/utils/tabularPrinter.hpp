#ifndef tabular_printer_hpp_
#define tabular_printer_hpp_
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>

void printTable(std::map<int, std::vector<std::string>> table, std::vector<std::string> headers, std::string delimiter)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(std::cout);
  std::vector<int> columnWidths(table.size(), 0);

  // for each column in the table, compute the maximum length string
  for(auto [columnNumber, entries] : table)
  {
    for(auto entry : entries)
    {
      if(entry.length() > columnWidths[columnNumber])
      {
        columnWidths[columnNumber] = entry.length();
      }
    }
  }

  // print the headers, padding each entry to the maximum length
  for(int i = 0; i < headers.size(); i++)
  {
    if(i != 0) std::cout << delimiter;
    std::cout << std::left
              << std::setw(columnWidths[i])
              << std::setfill(' ')
              << headers[i];
  }

  std::cout << "\n";

  // print contents of table
  for(int i = 0; i < table[0].size(); i++)
  {
    for(int j = 0; j < table.size(); j++)
    {
      if(j != 0) std::cout << delimiter;
      std::cout << std::left
                << std::setw(columnWidths[j])
                << std::setfill(' ')
                << table[j][i];
    }
    std::cout << "\n";
  }

  std::cout.copyfmt(oldState);
}
#endif