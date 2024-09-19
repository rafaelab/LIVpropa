#ifndef LIVPROPA_DATA_H
#define LIVPROPA_DATA_H

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include <kiss/path.h>
#include <kiss/logger.h>

#include "livpropa/Common.h"


namespace livpropa {


/*
 Returns the full path to the data files
*/
string getDataPath(string filename);

/* 
 Returns the location where LIVpropa was installed
*/
string  getInstallPrefix();


} // namespace livpropa


#endif // LIVPROPA_DATA_H
