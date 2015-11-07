#pragma once

#include "vec.h"
#include "parameters.h"

template<typename Real, typename Integer, Dimension Dim>
Particles<Real,Integer,Dim>::Particles(const Parameters<Real,Integer,Dim>& params): parameters{params},
                             pos{hemi::Array< Vec<Real,Dim> >(parameters.GetMaxParticlesLocal())},
                             pos_star{hemi::Array< Vec<Real,Dim> >(parameters.GetMaxParticlesLocal())},
                             v{hemi::Array< Vec<Real,Dim> >(parameters.GetMaxParticlesLocal())},
                             delta_pos{hemi::Array< Vec<Real,Dim> >(parameters.GetMaxParticlesLocal())},
                             vorticity{hemi::Array< Vec<Real,Dim> >(parameters.GetMaxParticlesLocal())},
                             color{hemi::Array< Vec<Real,Dim> >(parameters.GetMaxParticlesLocal())},
                             density{hemi::Array< Real >(parameters.GetMaxParticlesLocal())},
                             lambda{hemi::Array< Real >(parameters.GetMaxParticlesLocal())},
                             id{hemi::Array< Integer >(parameters.GetMaxParticlesLocal())}
{};

template<typename Real, typename Integer, Dimension Dim>
void ApplyGravity();

template< typename Real, typename Integer, Dimension Dim>
void ComputeDensities();

template< typename Real, typename Integer, Dimension Dim>
void PredictPositions();

template< typename Real, typename Integer, Dimension Dim>
void UpdatePositions();

template< typename Real, typename Integer, Dimension Dim>
void ComputeLambdas();

template< typename Real, typename Integer, Dimension Dim>
void UpdatePositionStars();

template< typename Real, typename Integer, Dimension Dim>
void UpdateDeltaPositions();

template< typename Real, typename Integer, Dimension Dim>
void UpdateVelocities();

template< typename Real, typename Integer, Dimension Dim>
void ApplyBoundaryConditions();

template< typename Real, typename Integer, Dimension Dim>
void ApplySurfaceTension();

template< typename Real, typename Integer, Dimension Dim>
void ApplyViscosity();

template< typename Real, typename Integer, Dimension Dim>
void ComputeVorticity();

template< typename Real, typename Integer, Dimension Dim>
void ApplyVorticityConfinement();
