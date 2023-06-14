//
// Created by finge on 6/13/2023.
//

#ifndef PULSERECONSTRUCTION_H5_FUNCTIONS_H
#define PULSERECONSTRUCTION_H5_FUNCTIONS_H

#include "H5Cpp.h"
#include "hdf5.h"
#include <string>
#include <vector>

std::vector<std::vector<double>> h5readVectord(std::string file, std::string var);

#endif //PULSERECONSTRUCTION_H5_FUNCTIONS_H
