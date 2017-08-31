#include <H5Cpp.h>
#include <H5Dpublic.h>
#include <vector>
#include <unordered_map>

namespace su {
    class biom {
        public:
            // cache the IDs contained within the table
            std::vector<std::string> sample_ids;
            std::vector<std::string> obs_ids;

            // cache both index pointers into both CSC and CSR representations
            std::vector<uint32_t> sample_indptr;
            std::vector<uint32_t> obs_indptr;

            uint32_t n_samples;  // the number of samples
            uint32_t n_obs;      // the number of observations
            uint32_t nnz;        // the total number of nonzero entries
            double *sample_counts;

            /* default constructor
             *
             * @param filename The path to the BIOM table to read
             */
            biom(std::string filename);

            /* copy constructor */
            biom(const biom &obj);

            /* default destructor
             *
             * Temporary arrays are freed
             */
            ~biom();

            /* get a dense vector of observation data
             *
             * @param id The observation ID to fetch
             * @param out An allocated array of at least size n_samples. 
             *      Values of an index position [0, n_samples) which do not
             *      have data will be zero'd.
             */
            void get_obs_data(std::string id, double* out);
        private:
            /* retain DataSet handles within the HDF5 file */
            H5::DataSet obs_indices;
            H5::DataSet sample_indices;
            H5::DataSet obs_data;
            H5::DataSet sample_data;
            H5::H5File file;
            uint32_t **obs_indices_resident;
            double **obs_data_resident;
            unsigned int *obs_counts_resident;

            unsigned int get_obs_data_direct(std::string id, uint32_t *& current_indices_out, double *& current_data_out);
            unsigned int get_sample_data_direct(std::string id, uint32_t *& current_indices_out, double *& current_data_out);
            double* get_sample_counts();

            /* At construction, lookups mapping IDs -> index position within an
             * axis are defined
             */
            std::unordered_map<std::string, uint32_t> obs_id_index;
            std::unordered_map<std::string, uint32_t> sample_id_index;
 
            /* load ids from an axis
             *
             * @param path The dataset path to the ID dataset to load
             * @param ids The variable representing the IDs to load into
             */          
            void load_ids(const char *path, std::vector<std::string> &ids);

            /* load the index pointer for an axis
             *
             * @param path The dataset path to the index pointer to load
             * @param indptr The vector to load the data into
             */
            void load_indptr(const char *path, std::vector<uint32_t> &indptr);

            /* count the number of nonzero values and set nnz */
            void set_nnz();

            /* create an index mapping an ID to its corresponding index 
             * position.
             *
             * @param ids A vector of IDs to index
             * @param map A hash table to populate
             */
            void create_id_index(std::vector<std::string> &ids, 
                                 std::unordered_map<std::string, uint32_t> &map);
    };
}
