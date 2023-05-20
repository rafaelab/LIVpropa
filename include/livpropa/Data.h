#ifndef LIVPROPA_DATA_H
#define LIVPROPA_DATA_H

#include <cstdlib>
#include <cstring>
#include <fstream>


#include <kiss/path.h>
#include <kiss/logger.h>



namespace livpropa {


/*
 Returns the full path to the data files
*/
std::string getDataPath(std::string filename);

/* 
 Returns the location where LIVpropa was installed
*/
std::string  getInstallPrefix();


} // namespace livpropa


#endif // LIVPROPA_DATA_H
