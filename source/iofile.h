#ifndef IOFILE
#define IOFILE


#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>

#include"databuffer.h"

bool loadAtomsFromNDL(dataBuffer &db,bool verbose=false);
bool loadAtomsFromXYZ(dataBuffer &db,bool verbose=false);
bool loadAtomsFromLammps(dataBuffer &db,bool verbose=false);
bool saveDataToXYZ(const std::string &filename, vdb *ptr_vdb);
bool saveDataToDMP(const std::string &filename, vdb *ptr_vdb);
bool saveDataToCONFIG(const std::string &filename, vdb *ptr_vdb,bool verbose=false);
bool saveDataToCOORD(const std::string &filename, vdb *ptr_vdb,bool verbose=false);
bool saveDataToLammps(const std::string &filename, vdb *ptr_vdb,bool verbose=false);
bool saveDataToD2S(const std::string &filename, vdb *ptr_vdb,bool verbose=false);

#endif // IOFILE

