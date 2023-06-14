//
// Created by finge on 6/13/2023.
//

#include "h5_functions.h"

using namespace H5;

std::vector<std::vector<double>> h5readVectord(std::string file, std::string var) {
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

    std::vector<std::vector<T>> datavec;
    datavec.resize(dims[0]);
    for (hsize_t i = 0; i < dims[0]; i++) {
        datavec[i].resize(dims[1]);
    }

    // Initialize hyperslabs
    hsize_t dataCount[2] = {1, dims[1]};
    hsize_t dataOffset[2] = {0, 0};
    const hsize_t memCount[1] = {dims[1]};
    const hsize_t memOffset[1] = {0};
    memspace.selectHyperslab(H5S_SELECT_SET, memCount, memOffset);
    for (hsize_t i = 0; i < dims[0]; i++) {
        dataOffset[0] = i;
        dataspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
        dataset.read(data[i].data(), PredType, memspace, dataspace);
    }
    return datavec;
}