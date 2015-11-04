#pragma once

#include "vec.h"
#include "hemi/array.h"
#include "hemi.parallel_for.h"
#include "parameters.h"

template< typename T, Dimension D >
Particles(const Parameters<T>& params, const int max_local);

template< typename T, Dimension D >
~Particles() {};

template< typename T, Dimension D >
void ApplyGravity();

template< typename T, Dimension D >
void ComputeDensities();

template< typename T, Dimension D >
void PredictPositions();

template< typename T, Dimension D >
void UpdatePositions();

template< typename T, Dimension D >
void ComputeLambdas();

template< typename T, Dimension D >
void UpdatePositionStars();

template< typename T, Dimension D >
void UpdateDeltaPositions();

template< typename T, Dimension D >
void UpdateVelocities();

template< typename T, Dimension D >
void ApplyBoundaryConditions();

template< typename T, Dimension D >
void ApplySurfaceTension();

template< typename T, Dimension D >
void ApplyViscosity();

template< typename T, Dimension D >
void ComputeVorticity();

template< typename T, Dimension D >
void ApplyVorticityConfinement();
