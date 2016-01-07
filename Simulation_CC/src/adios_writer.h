#pragma once

#include  <stdexcept>
#include "dimension.h"
#include "distributor.h"
#include "particles.h"

extern "C" {
  #include "adios.h"
  #include "adios_error.h"
}

/////////////////////////////////////////////////
// Adios file/stream writer class
/////////////////////////////////////////////////

template <typename Real, Dimension Dim>
class AdiosWriter {
public:
  /**
    @brief Constructor: initilize adios environment and required members
    @param[in] adios_writer_xml string to adios writer .xml file
    @param[in] distributor reference
    @param[in] particles reference
    @param[in] parameters reference
  **/
  AdiosWriter(const std::string& adios_writer_xml,
              const Distributor<Real,Dim>& distributor,
              const Particles<Real,Dim>& particles): distributor_{distributor},
                                                     particles_{particles}
  {
    int err = adios_init(adios_writer_xml.data(), static_cast<MPI_Comm>(distributor.compute_comm()));
    if(err) {
      // Not all ranks return err on failure and so
      // MPI_Finalize will deadlock when distributor is destructed
      // So we let adios print the error message and then exit instead of
      // throwing exception
      exit(1);
    }
  }

  /**
    @brief Destructor: finalize adios environment
  **/
  ~AdiosWriter(){
    adios_finalize(distributor_.compute_rank());
  }

  AdiosWriter(const AdiosWriter&)            =delete;
  AdiosWriter& operator=(const AdiosWriter&) =delete;
  AdiosWriter(AdiosWriter&&) noexcept        =delete;
  AdiosWriter& operator=(AdiosWriter&&)      =delete;

  /**
    @brief Write particle coordinates as described in adios_writer_xml
    As custom types are no supported by Adios vectors are written as raw bytes
  **/
  void write_particles() {
    int err;
    err = adios_open(&adios_handle_, "particles", "sim-output.bp", "a",
                         static_cast<MPI_Comm>(distributor_.compute_comm()));
    if(err)
      throw std::runtime_error(adios_get_last_errmsg());

    // Set ADIOS particle "group" size
    uint64_t group_bytes = 3 * sizeof(int64_t)                                  // global, local, offset count
                         + distributor_.resident_count() * sizeof(Vec<Real,Dim>);  // (x, y {,z})
    uint64_t total_bytes;
    err = adios_group_size(adios_handle_, group_bytes, &total_bytes);
    if(err)
      throw std::runtime_error(adios_get_last_errmsg());

    // Compute offset in global output for current rank(sum of ranks to the "left")
    uint64_t local_bytes = distributor_.resident_count() * sizeof(Vec<Real,Dim>);
    std::size_t offset_bytes = boost::mpi::scan(distributor_.compute_comm(), local_bytes, std::plus<std::size_t>());
    offset_bytes -= local_bytes;

    uint64_t global_bytes = distributor_.global_resident_count() * sizeof(Vec<Real,Dim>);

    err  = adios_write(adios_handle_, "global_bytes", &global_bytes);
    err |= adios_write(adios_handle_, "local_bytes", &local_bytes);
    err |= adios_write(adios_handle_, "offset_bytes", &offset_bytes);
    // adios_write takes a non-const pointer so we unsafely cast it away
    void* positions_ptr = static_cast<void*>(const_cast<Vec<Real,Dim>*>
                                            (particles_.positions().data()));

    err |= adios_write(adios_handle_, "positions", positions_ptr);

    if(err)
      throw std::runtime_error(adios_get_last_errmsg());

    err = adios_close(adios_handle_);
    if(err)
      throw std::runtime_error(adios_get_last_errmsg());
  }

private:
    int64_t adios_handle_;
    const Distributor<Real,Dim>& distributor_;
    const Particles<Real,Dim>& particles_;
};
