#ifndef ELEMENTS_H
#define ELEMENTS_H


#include <map>
#include <string>



namespace  Elements {

extern  std::map<std::wstring, std::wstring> mass;

void addElement(std::wstring & name__, std::wstring & mass__);

}

#endif // ELEMENTS_H
