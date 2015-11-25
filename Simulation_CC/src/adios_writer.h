#pragma once

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
              const Distributor<Real,Dim>& distributor_,
              const Particles<Real,Dim>& particles_,
              const Parameters<Real,Dim>& parameters_): m_distributor{distributor_},
                                                        m_particles{particles_},
                                                        m_parameters{parameters_}
  {
    int err = adios_init(adios_writer_xml.data(), static_cast<MPI_Comm>(m_distributor.GetComputeComm()));
    if(err)
      throw std::runtime_error(adios_get_last_errmsg());
  }

  /**
    @brief Destructor: finalize adios environment
  **/
  ~AdiosWriter(){
    adios_finalize(m_distributor.GetComputeRank());
  }

  AdiosWriter(const AdiosWriter&)            =delete;
  AdiosWriter& operator=(const AdiosWriter&) =delete;
  AdiosWriter(AdiosWriter&&) noexcept        =delete;
  AdiosWriter& operator=(AdiosWriter&&)      =delete;

  /**
    @brief Write particle coordinates as described in adios_writer_xml
    As custom types are no supported by Adios vectors are written as raw bytes
  **/
  void WriteParticles() {
    int err;
    err = adios_open(&m_adios_handle, "particles", "sim-output.bp", "a",
                         static_cast<MPI_Comm>(m_distributor.GetComputeComm()));
    if(err)
      throw std::runtime_error(adios_get_last_errmsg());

    // Set ADIOS particle "group" size
    uint64_t group_bytes = 3 * sizeof(int64_t)                                  // global, local, offset count
                         + m_particles.GetLocalCount() * sizeof(Vec<Real,Dim>);  // (x, y {,z})
    uint64_t total_bytes;
    err = adios_group_size(m_adios_handle, group_bytes, &total_bytes);
    if(err)
      throw std::runtime_error(adios_get_last_errmsg());

    // Compute offset in global output for current rank(sum of ranks to the "left")
    uint64_t local_bytes = m_particles.GetLocalCount() * sizeof(Vec<Real,Dim>);
    std::size_t offset_bytes = boost::mpi::scan(m_distributor.GetComputeComm(), local_bytes, std::plus<std::size_t>());
    offset_bytes -= local_bytes;

    uint64_t global_bytes = m_parameters.GetGlobalParticleCount() * sizeof(Vec<Real,Dim>);

    err  = adios_write(m_adios_handle, "global_bytes", &global_bytes);
    err |= adios_write(m_adios_handle, "local_bytes", &local_bytes);
    err |= adios_write(m_adios_handle, "offset_bytes", &offset_bytes);
    // adios_write takes a non-const pointer so we unsafely cast it away
    void* pos_pointer = static_cast<void*>(const_cast<Vec<Real,Dim>*>
                          (m_particles.GetPositionsPointer()));
    err |= adios_write(m_adios_handle, "xyz", pos_pointer);

    if(err)
      throw std::runtime_error(adios_get_last_errmsg());

    err = adios_close(m_adios_handle);
    if(err)
      throw std::runtime_error(adios_get_last_errmsg());
  }

private:
    int64_t m_adios_handle;
    const Distributor<Real,Dim>& m_distributor;
    const Particles<Real,Dim>& m_particles;
    const Parameters<Real,Dim>& m_parameters;
};
