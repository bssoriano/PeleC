#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include <AMReX_ErrorList.H>

struct TaggingParm
{
  amrex::Real denerr = 1.0e10;
  amrex::Real dengrad = 1.0e10;
  amrex::Real denratio = 1.0e10;
  amrex::Real presserr = 1.0e10;
  amrex::Real pressgrad = 1.0e10;
  amrex::Real velerr = 1.0e10;
  amrex::Real velgrad = 1.0e10;
  amrex::Real vorterr = 1.0e10;
  amrex::Real temperr = 1.0e10;
  amrex::Real lotemperr = -1.0e10;
  amrex::Real tempgrad = 1.0e10;
  amrex::Real ftracerr = 1.0e10;
  amrex::Real ftracgrad = 1.0e10;
  amrex::Real vfracerr = 1.0e10;
  amrex::Real detag_eb_factor = 3.0;

  std::string eb_refine_type = "static";

  int max_denerr_lev = 10;
  int max_dengrad_lev = 10;
  int max_denratio_lev = 10;
  int max_presserr_lev = 10;
  int max_pressgrad_lev = 10;
  int max_velerr_lev = 10;
  int max_velgrad_lev = 10;
  int max_vorterr_lev = 10;
  int max_temperr_lev = 10;
  int max_lotemperr_lev = 10;
  int max_tempgrad_lev = 10;
  int max_ftracerr_lev = 10;
  int max_ftracgrad_lev = 10;
  int max_vfracerr_lev = 10;


  amrex::Real dom_lbnd_x = 1.0e10;
  amrex::Real dom_hbnd_x = 1.0e10;
  amrex::Real dom_lbnd_y = 1.0e10;
  amrex::Real dom_hbnd_y = 1.0e10;
  amrex::Real dom_lbnd_z = 1.0e10;
  amrex::Real dom_hbnd_z = 1.0e10;
  int max_dom_bnd_lev = 10;

  int max_eb_refine_lev = 10;
  int min_eb_refine_lev = -1;
  int adapt_eb_refined_lev = -1;

  amrex::Vector<amrex::AMRErrorTag> err_tags;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error_region(
  const int i,
  const int j,
  const int k,
  amrex::GeometryData const& geomdata,
  amrex::Real x_lo,
  amrex::Real x_hi,
  amrex::Real y_lo,
  amrex::Real y_hi,
  amrex::Real z_lo,
  amrex::Real z_hi,
  amrex::Array4<char> const& tag,
  char tagval) noexcept
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  // const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];


  // Tag on regions of high field values
  if ((x_lo < x) && (x < x_hi) && (y_lo < y) && (y < y_hi)
#if AMREX_SPACEDIM == 3
      && (z_lo < z) && (z < z_hi)
  //    (dom_lbnd[2] < z) && (z < dom_hbnd[2])) {
#endif
     ){
        tag(i, j, k) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if ((field(iv) >= fielderr) && (!flags(iv).isCovered())) {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_loerror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of low field values
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if ((field(iv) <= fielderr) && (!flags(iv).isCovered())) {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_graderror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fieldgrad,
  char tagval) noexcept
{
  // Tag on regions of high field gradient
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> as = {{0.0}};

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    const auto ivp = iv + amrex::IntVect::TheDimensionVector(idir);
    const auto ivm = iv - amrex::IntVect::TheDimensionVector(idir);
    const bool cp = flags(iv).isCovered() || flags(ivp).isCovered();
    const bool cm = flags(iv).isCovered() || flags(ivm).isCovered();
    const amrex::Real ap = !cp ? std::abs(field(ivp) - field(iv)) : 0.0;
    const amrex::Real am = !cm ? std::abs(field(iv) - field(ivm)) : 0.0;
    as[idir] = amrex::max<amrex::Real>(ap, am);
  }

#if AMREX_SPACEDIM > 1
  if (amrex::max<amrex::Real>(AMREX_D_DECL(as[0], as[1], as[2])) >= fieldgrad)
#elif AMREX_SPACEDIM == 1
  if (as[0] >= fieldgrad)
#endif
  {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_ratioerror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fieldratio,
  char tagval) noexcept
{
  // Tag on regions of high ratio of adjacent cell values in field
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> as = {{0.0}};

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    const auto ivp = iv + amrex::IntVect::TheDimensionVector(idir);
    const auto ivm = iv - amrex::IntVect::TheDimensionVector(idir);
    const bool cp = flags(iv).isCovered() || flags(ivp).isCovered();
    const bool cm = flags(iv).isCovered() || flags(ivm).isCovered();
    amrex::Real ap = !cp ? std::abs(field(ivp) / field(iv)) : 0.0;
    ap = !cp ? amrex::max<amrex::Real>(ap, 1.0 / ap) : 0.0;
    amrex::Real am = !cm ? std::abs(field(ivm) / field(iv)) : 0.0;
    am = !cm ? amrex::max<amrex::Real>(am, 1.0 / am) : 0.0;
    as[idir] = amrex::max<amrex::Real>(am, ap);
  }

#if AMREX_SPACEDIM > 1
  if (amrex::max<amrex::Real>(AMREX_D_DECL(as[0], as[1], as[2])) >= fieldratio)
#elif AMREX_SPACEDIM == 1
  if (as[0] >= fieldratio)
#endif
  {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_abserror(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real fielderr,
  char tagval) noexcept
{
  // Tag on regions of high field values
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if ((std::abs(field(iv)) >= fielderr) && (!flags(iv).isCovered())) {
    tag(iv) = tagval;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
tag_error_bounds(
  const int i,
  const int j,
  const int k,
  amrex::Array4<amrex::EBCellFlag const> const& flags,
  amrex::Array4<char> const& tag,
  amrex::Array4<amrex::Real const> const& field,
  const amrex::Real lbnd,
  const amrex::Real ubnd,
  char tagval) noexcept
{
  // Tag on regions inside bounds
  amrex::ignore_unused(k);
  const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
  if (((lbnd < field(iv)) && (field(iv) < ubnd)) && (!flags(iv).isCovered())) {
    tag(iv) = tagval;
  }
}

#endif
