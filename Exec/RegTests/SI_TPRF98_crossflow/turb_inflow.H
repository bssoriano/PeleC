#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Forcing.H"
#include "Utilities.H"

void
read_turb_input_file(
  const std::string iname,
  int z_coord,
  int nx,
  int ny,
  int nz,
  int nscal,
  int ix_dir_limit[2],
  amrex::Vector<amrex::Real>& gridx_input,
  amrex::Vector<amrex::Real>& gridy_input,
  amrex::Vector<amrex::Real>& gridz_input,
  amrex::Vector<amrex::Real>& data);


void get_inputs_turb(
  const std::string& iname,
  int z_coord);


void read_turbulent_bc(int ix_limit[2], amrex::Real v_mean);