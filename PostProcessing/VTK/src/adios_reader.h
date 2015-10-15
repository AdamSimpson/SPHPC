#ifndef VTK_SRC_ADIOS_READER_H_
#define VTK_SRC_ADIOS_READER_H_

#include <vector>
#include <string>
#include <boost/mpi.hpp>

extern "C" {
  #include "adios_read.h"
}

class AdiosReader {
  public:
    AdiosReader(const std::string& name, const boost::mpi::communicator& comm);
    ~AdiosReader();

    template<typename T>
    std::vector<T> FetchValue(const std::string& value_name);

    AdiosReader(const AdiosReader&)            =delete;
    AdiosReader& operator=(const AdiosReader&) =delete;
    AdiosReader(AdiosReader&&) noexcept        =delete;
    AdiosReader& operator=(AdiosReader&&)      =delete;

  private:
    std::string file_name;
    ADIOS_FILE *adios_file;
    ADIOS_SELECTION *adios_selection = NULL;
    boost::mpi::communicator communicator;
};

#endif
