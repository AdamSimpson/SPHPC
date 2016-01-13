#pragma once

#include <cmath>
#include "dimension.h"
#include "array.h"
#include "vec.h"
#include "parameters.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/for_each.h"
#include "thrust/iterator/zip_iterator.h"
#include "thrust/sort.h"
#include "thrust/binary_search.h"
#include "thrust/execution_policy.h"

#define MAX_NEIGHBORS 60

struct NeighborBin {
  std::size_t neighbor_indices[MAX_NEIGHBORS];
  int count;
};

// Iterator for range based for loops over neighbor indices
const std::size_t* begin(const NeighborBin& bin) {
  return bin.neighbor_indices;
}

const std::size_t* end(const NeighborBin& bin) {
  return bin.neighbor_indices + bin.count;
}

template<typename Real, Dimension Dim>
class Neighbors {
public:
  Neighbors(const Parameters<Real,Dim>& parameters): parameters_{parameters},
                                                     bin_spacing_{parameters.smoothing_radius()},
                                                     bin_dimensions_{static_cast<Vec<std::size_t,Dim>>(ceil(
                                                                                                            (parameters.boundary().extent())
                                                                                                            /bin_spacing_) + static_cast<Real>(2))},
                                                     begin_indices_{product(bin_dimensions_)},
                                                     end_indices_{product(bin_dimensions_)},
                                                     bin_ids_{parameters.max_particles_local()},
                                                     particle_ids_{parameters.max_particles_local()},
                                                     neighbor_bins_{parameters.max_particles_local()} {};

// Access neighbor bins with subscript operator
const NeighborBin& operator[] (const std::size_t index) const {
  return neighbor_bins_[index];
}

/**
  return linearized bin value of Vec<Real,2>
  Each point is shifted by bin_spacing in each direction
  To account for a 1 block grid boundary. This boundary allows neighbors to easily
  be searched for without worrying about boundary conditions
 
**/
std::size_t calculate_bin_id(const Vec<Real,2>& point) const {
  const auto point_shifted = point + bin_spacing_;
  const auto bin_location = static_cast< Vec<std::size_t, two_dimensional> >(floor(point_shifted/bin_spacing_));
  return (bin_location.y * bin_dimensions_.x + bin_location.x);
}

/**
  return linearized bin value of Vec<Real,3>
  Each point is shifted by bin_spacing in each direction
  To account for a 1 block grid boundary. This boundary allows neighbors to easily
  be searched for without worrying about boundary conditions
 
**/
std::size_t calculate_bin_id(const Vec<Real,3>& point) const {
  const auto point_shifted = point + bin_spacing_;
  const auto bin_location = static_cast< Vec<std::size_t, three_dimensional> >(floor(point_shifted/bin_spacing_));
  return (bin_dimensions_.x * bin_dimensions_.y * bin_location.z) +
    (bin_location.y * bin_dimensions_.x + bin_location.x);
}

  /**
     Calclate bin value for each particle
  **/
  void calculate_bins(const IndexSpan& span, const Vec<Real,Dim>* position_stars) {
    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(thrust::device, begin, end, [=] (std::size_t i) {
      bin_ids_[i] = this->calculate_bin_id(position_stars[i]);
      particle_ids_[i] = i;
      });
  }

  /**
     Sort particle_ids array into bin order
   **/
  void sort_bins(std::size_t particle_count) {
    thrust::sort_by_key(thrust::device, bin_ids_.data(), bin_ids_.data() + particle_count , particle_ids_.data());
  }

  /**
     Find begin/end bounds of bins
   **/
  void find_bin_bounds(std::size_t particle_count) {
    auto counting_begin = thrust::make_counting_iterator((std::size_t)0); 
    auto counting_end = thrust::make_counting_iterator(product(bin_dimensions_)); 

    thrust::lower_bound(thrust::device, bin_ids_.data(), bin_ids_.data() + particle_count,
                        counting_begin, counting_end,
                        begin_indices_.data());
    thrust::upper_bound(thrust::device, bin_ids_.data(), bin_ids_.data() + particle_count,
                        counting_begin, counting_end,
                        end_indices_.data());
  }

  /**
     Calculate neighbor bin indices based upon particle coordinate
     @todo use array_span
   **/
  void calculate_neighbor_indices(const Vec<Real, 2>& coord, std::size_t* neighbor_indices) {
    int index = 0;
    for(int i=-1; i<2; i++) {
      for(int j=-1; j<2; j++) {
        const Vec<Real,2> neighbor_coord{coord.x + i*bin_spacing_ ,
                                         coord.y + j*bin_spacing_};
        neighbor_indices[index] = calculate_bin_id(neighbor_coord);
        ++index;
      }
    }
  }

  /**
     Calculate neighbor bin indices based upon particle coordinate
     @todo use array_span
  **/
  void calculate_neighbor_indices(const Vec<Real, 3>& coord, std::size_t* neighbor_indices) {
    int index = 0;
    for(int i=-1; i<2; i++) {
      for(int j=-1; j<2; j++) {
        for(int k=-1; k<2; k++) {
          const Vec<Real,3> neighbor_coord{coord.x + i*bin_spacing_ ,
                                           coord.y + j*bin_spacing_ ,
                                           coord.z + k*bin_spacing_};
          neighbor_indices[index] = calculate_bin_id(neighbor_coord);
          ++index;
        }
      }
    }
  }

  constexpr std::size_t neighbor_count() const {
    return static_cast<std::size_t>(pow(3, Dim));
  }

  /**
     Fill the neighbor bins in the specified particle span
   **/
  void fill_neighbors(IndexSpan span, const Vec<Real,Dim>* position_stars) {
    //    auto index = calculate_bin_id(Vec<Real,3>{55.0,55.0,55.0});
    //    std::cout<<"index "<<index<<"bins: "<<neighbor_bins_[index].count<<std::endl;
    
    const Real smoothing_radius_squared = parameters_.smoothing_radius() *  parameters_.smoothing_radius();

    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(thrust::device, begin, end, [=] (std::size_t particle_index) {
        const auto position_star = position_stars[particle_index];
        auto& bin = neighbor_bins_[particle_index];
        std::size_t neighbor_bin_indices[neighbor_count()];

        // Zero out neighbor count for current bin
        neighbor_bins_[particle_index].count = 0;

        calculate_neighbor_indices(position_stars[particle_index], neighbor_bin_indices);

        for(auto neighbor_bin_index : neighbor_bin_indices) {
          const auto begin_index = begin_indices_[neighbor_bin_index];
          const auto end_index   = end_indices_[neighbor_bin_index];
          for(auto j = begin_index; j < end_index; ++j) {
            const std::size_t neighbor_particle_index = particle_ids_[j];
            if(particle_index == neighbor_particle_index)
              continue;

            const auto neighbor_position_star = position_stars[neighbor_particle_index];
            const Real distance_squared = magnitude_squared(position_star - neighbor_position_star);
            if(distance_squared < smoothing_radius_squared && bin.count < MAX_NEIGHBORS) {
              bin.neighbor_indices[bin.count] = neighbor_particle_index;
              ++bin.count;
            }
          }
        }
    });
  }

  /**
    Find all particle neighbors.
    particles_to_bin_span is the span of particles to bin/sort
    particles_to_fill_span is the span of particles to fill in neighbors for
    This is seperate as we need to bin/sort resident + halo but
    Only need to fill neighbors for resident
  **/
  void find(const IndexSpan& particles_to_bin_span, const IndexSpan& particles_to_fill_span,
            const Vec<Real,Dim>* coords) {
    const auto particles_to_bin_count = particles_to_bin_span.end - particles_to_bin_span.begin;
    this->calculate_bins(particles_to_bin_span, coords);
    this->sort_bins(particles_to_bin_count);
    this->find_bin_bounds(particles_to_bin_count);
    this->fill_neighbors(particles_to_bin_span, coords);
  }

  ~Neighbors()                           = default;
  Neighbors(const Neighbors&)            = default;
  Neighbors& operator=(const Neighbors&) = default;
  Neighbors(Neighbors&&) noexcept        = default;
  Neighbors& operator=(Neighbors&&)      = default;

  sim::Array<NeighborBin> neighbor_bins_;

private:
  const Parameters<Real,Dim>& parameters_;

  Real bin_spacing_;
  Vec<std::size_t, Dim> bin_dimensions_;

  sim::Array<std::size_t> begin_indices_; // Begin indices for bin ids
  sim::Array<std::size_t> end_indices_;   // End indices for bin ids
  sim::Array<std::size_t> bin_ids_;   // Array of bin ids
  sim::Array<std::size_t> particle_ids_;  // Array of particle ids
};
