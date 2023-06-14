//
// Created by finge on 6/13/2023.
//

#include "h5_functions.h"

using namespace H5;

template <typename T> T h5read(std::string file, std::string var) {
    auto item_type = PredType::NATIVE_DOUBLE;
    auto mem_type = VarLenType(&item_type);

    H5File datafile(file.c_str(), H5F_ACC_READONLY);
    DataSet dataset = datafile.openDataSet(var.c_str());
    DataSpace dataspace = dataset.getSpace();

    hsize_t dims[2];
    hsize_t rank = dataspace.getSimpleExtentDims(dims, NULL);

    hsize_t dimsm[1];
    dimsm[0] = dims[1];
    DataSpace memspace(1, dimsm);


}