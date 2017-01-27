#include <iostream>
#include "biom.hpp"

using namespace H5;
using namespace su;

/* datasets defined by the BIOM 2.x spec */ 
const std::string OBS_INDPTR = std::string("/observation/matrix/indptr");
const std::string OBS_INDICES = std::string("/observation/matrix/indices");
const std::string OBS_DATA = std::string("/observation/matrix/data");
const std::string OBS_IDS = std::string("/observation/ids");

const std::string SAMPLE_INDPTR = std::string("/sample/matrix/indptr");
const std::string SAMPLE_INDICES = std::string("/sample/matrix/indices");
const std::string SAMPLE_DATA = std::string("/sample/matrix/data");
const std::string SAMPLE_IDS = std::string("/sample/ids");

biom::biom(std::string filename) {
    file = H5File(filename, H5F_ACC_RDONLY);

    /* establish the datasets */
    obs_indices = file.openDataSet(OBS_INDICES.c_str());
    obs_data = file.openDataSet(OBS_DATA.c_str());
    sample_indices = file.openDataSet(SAMPLE_INDICES.c_str());
    sample_data = file.openDataSet(SAMPLE_DATA.c_str());
    
    /* cache IDs and indptr */
    sample_ids = std::vector<std::string>();
    obs_ids = std::vector<std::string>();
    sample_indptr = std::vector<uint32_t>();
    obs_indptr = std::vector<uint32_t>();

    load_ids(OBS_IDS.c_str(), obs_ids);
    load_ids(SAMPLE_IDS.c_str(), sample_ids);
    load_indptr(OBS_INDPTR.c_str(), obs_indptr);    
    load_indptr(SAMPLE_INDPTR.c_str(), sample_indptr);    

    /* cache shape and nnz info */
    n_samples = sample_ids.size();
    n_obs = obs_ids.size();
    set_nnz();

    /* define reusable temporary arrays
     *
     * NOT THREAD SAFE
     */
    tmp_obs_indices = (uint32_t*)malloc(sizeof(uint32_t) * n_samples);
    tmp_sample_indices = (uint32_t*)malloc(sizeof(uint32_t) * n_obs);
    tmp_obs_data = (double*)malloc(sizeof(double) * n_samples);
    tmp_sample_data = (double*)malloc(sizeof(double) * n_obs);

    /* define a mapping between an ID and its corresponding offset */
    obs_id_index = std::unordered_map<std::string, uint32_t>();
    sample_id_index = std::unordered_map<std::string, uint32_t>();

    create_id_index(obs_ids, obs_id_index);
    create_id_index(sample_ids, sample_id_index);
}

biom::~biom() {
    /* free the reusable arrays */
    free(tmp_obs_indices);
    free(tmp_sample_indices);
    free(tmp_obs_data);
    free(tmp_sample_data);
}

void biom::set_nnz() {
    // should these be cached?
    DataType dtype = obs_data.getDataType();
    DataSpace dataspace = obs_data.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);
    nnz = dims[0];
}

void biom::load_ids(const char *path, std::vector<std::string> &ids) {
    DataSet ds_ids = file.openDataSet(path);
    DataType dtype = ds_ids.getDataType();
    DataSpace dataspace = ds_ids.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);

    /* the IDs are a dataset of variable length strings */
    char **dataout = (char**)malloc(sizeof(char*) * dims[0]);
    ds_ids.read((void*)dataout, dtype);

    ids.reserve(dims[0]);
    for(unsigned int i = 0; i < dims[0]; i++) {
        ids.push_back(dataout[i]);
    }
    
    for(unsigned int i = 0; i < dims[0]; i++) 
        free(dataout[i]);
    free(dataout);
}

void biom::load_indptr(const char *path, std::vector<uint32_t> &indptr) {
    DataSet ds = file.openDataSet(path);
    DataType dtype = ds.getDataType();
    DataSpace dataspace = ds.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);
    
    uint32_t *dataout = (uint32_t*)malloc(sizeof(uint32_t) * dims[0]);
    ds.read((void*)dataout, dtype);

    indptr.reserve(dims[0]);
    for(unsigned int i = 0; i < dims[0]; i++)
        indptr.push_back(dataout[i]);
    free(dataout);
}

void biom::create_id_index(std::vector<std::string> &ids, 
                           std::unordered_map<std::string, uint32_t> &map) {
    uint32_t count = 0;
    map.reserve(ids.size());
    for(auto i = ids.begin(); i != ids.end(); i++, count++) {
        map[*i] = count;
    }
}


void biom::get_obs_data(std::string id, double* out) {
    uint32_t idx = obs_id_index.at(id);
    uint32_t start = obs_indptr[idx];
    uint32_t end = obs_indptr[idx + 1];
 
    hsize_t count[1] = {end - start};
    hsize_t offset[1] = {start};

    DataType indices_dtype = obs_indices.getDataType();
    DataType data_dtype = obs_data.getDataType();

    DataSpace indices_dataspace = obs_indices.getSpace();
    DataSpace data_dataspace = obs_data.getSpace();
    
    DataSpace indices_memspace(1, count, NULL);
    DataSpace data_memspace(1, count, NULL);

    indices_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 
    data_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 

    // reset our output buffer
    for(unsigned int i = 0; i < n_samples; i++)
        out[i] = 0.0;
    
    // NOT THREAD SAFE due to use of tmp_obs_*.
    obs_indices.read((void*)tmp_obs_indices, indices_dtype, indices_memspace, indices_dataspace);
    obs_data.read((void*)tmp_obs_data, data_dtype, data_memspace, data_dataspace);

    for(unsigned int i = 0; i < count[0]; i++) {
        out[tmp_obs_indices[i]] = tmp_obs_data[i];
    }
}

void biom::get_sample_data(std::string id, double* out) {
    uint32_t idx = sample_id_index.at(id);
    uint32_t start = sample_indptr[idx];
    uint32_t end = sample_indptr[idx + 1];

    hsize_t count[1] = {end - start};
    hsize_t offset[1] = {start};

    DataType indices_dtype = sample_indices.getDataType();
    DataType data_dtype = sample_data.getDataType();

    DataSpace indices_dataspace = sample_indices.getSpace();
    DataSpace data_dataspace = sample_data.getSpace();
    
    DataSpace indices_memspace(1, count, NULL);
    DataSpace data_memspace(1, count, NULL);

    indices_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 
    data_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 

    // reset our output buffer
    for(unsigned int i = 0; i < n_obs; i++)
        out[i] = 0.0;
    
    // NOT THREAD SAFE due to use of tmp_sample_*. 
    sample_indices.read((void*)tmp_sample_indices, indices_dtype, indices_memspace, indices_dataspace);
    sample_data.read((void*)tmp_sample_data, data_dtype, data_memspace, data_dataspace);

    for(unsigned int i = 0; i < count[0]; i++) {
        out[tmp_sample_indices[i]] = tmp_sample_data[i];
    }
}
