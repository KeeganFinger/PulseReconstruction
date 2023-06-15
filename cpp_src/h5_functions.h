//
// Created by finge on 6/13/2023.
//

#ifndef PULSERECONSTRUCTION_H5_FUNCTIONS_H
#define PULSERECONSTRUCTION_H5_FUNCTIONS_H

#include "H5Cpp.h"
#include "hdf5.h"

using namespace H5;

template <typename T>
void h5readMatrixD(H5FILE &datafile, std::string &var, T &datavec) {
	auto item_type = PredType::NATIVE_DOUBLE;
	auto mem_type = VarLenType(&item_type);

	DataSet dataset = datafile.openDataSet(var.c_str());
	DataSpace dataspace = dataset.getSpace();

	hsize_t dims[2];
	hsize_t rank = dataspace.getSimpleExtentDims(dims, NULL);
	
	hsize_t dimsm[1];
	dimsm[0] = dims[1];
	DataSpace memspace(1, dimsm);

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
		dataset.read(datavec[i].data(), item_type, memspace, dataspace);
	}
}

template <typename T>
void h5readVectorD(H5FILE &datafile, std::string &var, T &datavec) {
	auto item_type = PredType::NATIVE_DOUBLE;
        auto mem_type = VarLenType(&item_type);

        DataSet dataset = datafile.openDataSet(var.c_str());
        DataSpace dataspace = dataset.getSpace();

        hsize_t dims[2];
        hsize_t rank = dataspace.getSimpleExtentDims(dims, NULL);

        hsize_t dimsm[1];
        dimsm[0] = dims[0];
        DataSpace memspace(1,dimsm);

        datavec.resize(dims[0]);

        // Initialize hyperslabs
        hsize_t dataCount[1] = {dims[0]};
        hsize_t dataOffset[1] = {0};
        const hsize_t memCount[1] = {dims[0]};
        const hsize_t memOffset[1] = {0};
        memspace.selectHyperslab(H5S_SELECT_SET, memCount, memOffset);
        dataspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
        dataset.read(datavec.data(), item_type, memspace, dataspace);
}

template <typename T>
void h5readScalarI(H5FILE &datafile, std::string &var, T &data) {
	auto item_type = PredType::NATIVE_INT;
	auto mem_type = VarLenType(&item_type);

	DataSet dataset = datafile.openDataSet(var.c_str());
	DataSpace dataspace = dataset.getSpace();

	hsize_t dims[2];
	hsize_t rank = dataspace.getSimpleExtentDims(dims, NULL);

	hsize_t dimsm[1];
	dimsm[0] = dims[0];
	DataSpace memspace(1,dimsm);

	// Initialize hyperslabs
	hsize_t dataCount[1] = {dims[0]};
	hsize_t dataOffset[1] = {0};
	const hsize_t memCount[1] = {dims[0]};
	const hsize_t memOffset[1] = {0};
	memspace.selectHyperslab(H5S_SELECT_SET, memCount, memOffset);
	dataspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
	dataset.read(&data, item_type, memspace, dataspace);
}

template <typename T>
void h5readScalarD(H5FILE &datafile, std::string &var, T &data) {
    auto item_type = PredType::NATIVE_DOUBLE;
    auto mem_type = VarLenType(&item_type);

    DataSet dataset = datafile.openDataSet(var.c_str());
    DataSpace dataspace = dataset.getSpace();

    hsize_t dims[2];
    hsize_t rank = dataspace.getSimpleExtentDims(dims, NULL);

    hsize_t dimsm[1];
    dimsm[0] = dims[0];
    DataSpace memspace(1,dimsm);

    // Initialize hyperslabs
    hsize_t dataCount[1] = {dims[0]};
    hsize_t dataOffset[1] = {0};
    const hsize_t memCount[1] = {dims[0]};
    const hsize_t memOffset[1] = {0};
    memspace.selectHyperslab(H5S_SELECT_SET, memCount, memOffset);
    dataspace.selectHyperslab(H5S_SELECT_SET, dataCount, dataOffset);
    dataset.read(&data, item_type, memspace, dataspace);
}

#endif //PULSERECONSTRUCTION_H5_FUNCTIONS_H
