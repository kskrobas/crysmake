#include "elements.h"
#include <iostream>

namespace  Elements {

 std::map<std::wstring, std::wstring> mass{

        {L"Al", L"26.982"},
        {L"H",  L"1.00794"},
        {L"C",  L"12.0107"},
        {L"Ca", L"40.078"},
        {L"Cd", L"112.411"},
        {L"Ce", L"140.116"},
        {L"Cr", L"51.996"},
        {L"Cu", L"63.546"},
        {L"Fe", L"55.845"},
        {L"F",  L"18.9984"},
        {L"Ga", L"69.723"},
        {L"Mn", L"54.94"},
        {L"Mo", L"95.94"},
        {L"N",  L"14.007"},
        {L"Na", L"22.9897"},
        {L"Ni", L"58.6934"},
        {L"Mg", L"24.305"},
        {L"O",  L"15.9994"},
        {L"P",  L"30.973762"},
        {L"S",  L"32.065"},
        {L"Se", L"78.96"},
        {L"Si", L"28.0855"},
        {L"U",  L"238.028913"},
        {L"W",  L"183.84"},
        {L"Te", L"127.6"},
        {L"Ti", L"47.867"},
        {L"Xe", L"131.29"},
        {L"Zn", L"65.32"},
         {L"Zr", L"91.22"},
        {L"H2O",L"18.015"}

 };

 void addElement(std::wstring & name__, std::wstring & mass__)
 {

     auto search = mass.find(name__);
            if (search != mass.end()) {
                std::wcout<<"WARNING: specified atom type is already defined  "<<name__<<std::endl;
                std::wcout<<"         value will be overwritten"<<std::endl;
            }

            mass.emplace(name__,mass__);
 }


}
