#include <H5Cpp.h>
#include <H5Dpublic.h>
#include <vector>
#include <unordered_map>

namespace su {
    class biom {
        public:
            std::vector<std::string> sample_ids;
            std::vector<std::string> obs_ids;

            std::vector<uint32_t> sample_indptr;
            std::vector<uint32_t> obs_indptr;

            uint32_t n_samples;
            uint32_t n_obs;
            uint32_t nnz;

            biom(std::string filename);
            ~biom();

            // need to store dataset attributes
            // need to store sample counts to support convert to proportions
            void get_obs_data(std::string id, double* out);
            void get_sample_data(std::string id, double* out);
        private:
            H5::DataSet obs_indices;
            H5::DataSet sample_indices;
            H5::DataSet obs_data;
            H5::DataSet sample_data;
            H5::H5File file;
            
            // not thread safe
            uint32_t *tmp_obs_indices;
            uint32_t *tmp_sample_indices;
            double *tmp_obs_data;
            double *tmp_sample_data;
            // end not thread safe

            std::unordered_map<std::string, uint32_t> obs_id_index;
            std::unordered_map<std::string, uint32_t> sample_id_index;
            
            void load_ids(const char *path, std::vector<std::string> &ids);
            void load_indptr(const char *path, std::vector<uint32_t> &indptr);
            void set_nnz();
            void create_id_index(std::vector<std::string> &ids, 
                                 std::unordered_map<std::string, uint32_t> &map);
    };
}
