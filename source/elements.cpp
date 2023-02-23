#include "elements.h"
#include <iostream>

namespace  Elements {

 std::map<std::wstring, std::wstring> mass{

        {L"H",  L"1.00794"},
        {L"C",  L"12.0107"},
        {L"Cd", L"112.411"},
        {L"Ga", L"69.723"},
        {L"N",  L"14.007"},
        {L"Ni", L"58.694"},
        {L"P", L"30.973762"},
        {L"S",  L"32.065"},
        {L"Se", L"78.96"},
        {L"Si", L"28.0855"},
        {L"Te", L"127.6"},
        {L"Zn", L"65.32"},
        {L"H2O", L"18.015"}

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
