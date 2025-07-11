#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  amrex::ParmParse pp("prob");
  pp.query("inject_fuel", PeleC::h_prob_parm_device->inject_fuel);
  pp.query("centx", PeleC::h_prob_parm_device->centx);
  pp.query("centz", PeleC::h_prob_parm_device->centz);
  pp.query("r_hole", PeleC::h_prob_parm_device->r_hole);
  pp.query("init_type", PeleC::h_prob_parm_device->init_type);
  pp.query("cavity_depth", PeleC::h_prob_parm_device->cavity_depth);
  
  //Read inflow EB jet data
  amrex::Vector<amrex::Real> jet_ct(2, 0);
  pp.queryarr("jet_center_EB", jet_ct, 0, 2);
  for (int idim = 0; idim < 2; idim++) {
    PeleC::h_prob_parm_device->jet_center_EB[idim] = jet_ct[idim];
  }
  pp.query("jet_radius_EB", PeleC::h_prob_parm_device->jet_radius_EB);
  pp.query("jet_vel_EB", PeleC::h_prob_parm_device->jet_vel_EB);
  pp.query("jet_T_EB", PeleC::h_prob_parm_device->jet_T_EB);

}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
}
