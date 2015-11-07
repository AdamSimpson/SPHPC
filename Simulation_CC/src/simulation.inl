#include <string>

template<typename Real, typename Integer, Dimension Dim>
Simulation<Real,Integer,Dim>::Simulation(const std::string& input_file_name):
  distributor{Distributor<Real,Integer,Dim>()},
  parameters{Parameters<Real,Integer,Dim>(input_file_name)},
  particles{Particles<Real,Integer,Dim>(parameters)}
{
//  distributor.SetParameters(parameters);
}
