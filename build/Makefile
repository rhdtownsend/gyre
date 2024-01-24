# File     : Makefile
# Purpose  : makefile for build subdir

# Build tool executables
TOOLS ?= yes

# Enable debugging (with a performance penalty)
DEBUG ?= no

# Enable OpenMP parallelization
OMP ?= yes

# Enable FPE checks
FPE ?= yes

# Enable correctly-rounded math
CRMATH?=yes

# Enable portable math (for bit-for-bit reproducibility; setting to
# no may give a performance boost)
PORTABLE ?= yes

# Use IEEE fortran features
IEEE ?= yes

# Enable experimental features
EXPERIMENTAL ?= no

############ DO NOT EDIT BELOW THIS LINE ############
### (unless you think you know what you're doing) ###
#####################################################

# Targets

BIN_TARGETS = gyre gyre_tides
LIB_TARGETS = libgyre

ifeq ($(TOOLS),yes)
  BIN_TARGETS := $(BIN_TARGETS) eval_lambda
endif

# Source file lists & related variables
#
# xxx_SRCS = source files required to compile xxx
# xxx_USES = source files providing modules USE'd by xxx
# xxx_LIBS = library files required to link xxx

gyre_SRCS = gyre.fypp
gyre_USES = gyre_m.fypp
gyre_LIBS = libgyre

gyre_tides_SRCS = gyre_tides.fypp
gyre_tides_USES = gyre_m.fypp
gyre_tides_LIBS = libgyre

eval_lambda_SRCS = eval_lambda.fypp
eval_lambda_USES = gyre_m.fypp
eval_lambda_LIBS = libgyre

libgyre_SRCS = ad_bound_m.fypp ad_bvp_m.fypp ad_diff_m.fypp		\
   ad_eqns_m.fypp ad_match_m.fypp ad_trans_m.fypp amdl_file_m.fypp	\
   atmos_m.fypp b3_file_m.fypp band_sysmtx_m.fypp block_sysmtx_m.fypp	\
   bound_m.fypp bracket_search_m.fypp bvp_m.fypp cfunc_m.fypp		\
   cheb_fit_m.fypp clinalg_m.fypp colloc_diff_m.fypp constants_m.fypp	\
   context_m.fypp contour_map_m.fypp contour_path_m.fypp		\
   contour_search_m.fypp deriv_m.fypp detail_m.fypp			\
   diff_factory_m.fypp diff_m.fypp discrim_m.fypp eqns_m.fypp		\
   evol_model_m.fypp ext_complex_m.fypp ext_m.fypp ext_real_m.fypp	\
   famdl_file_m.fypp fgong_file_m.fypp freq_context_m.fypp		\
   freq_m.fypp func_m.fypp gyre_m.fypp grid_factory_m.fypp		\
   grid_m.fypp grid_par_m.fypp grid_spec_m.fypp grid_util_m.fypp	\
   grid_weights_m.fypp gsm_file_m.fypp hdf_writer_m.fypp		\
   hgroup_m.fypp hom_model_m.fypp interp_m.fypp kinds_m.fypp		\
   linalg_m.fypp losc_file_m.fypp math_m.fypp magnus_diff_m.fypp	\
   memory_m.fypp mesa_file_m.fypp minim_m.fypp minmod_search_m.fypp	\
   mirk_diff_m.fypp mode_m.fypp mode_par_m.fypp model_factory_m.fypp	\
   model_m.fypp model_par_m.fypp model_util_m.fypp nad_bound_m.fypp	\
   nad_bvp_m.fypp nad_diff_m.fypp nad_eqns_m.fypp nad_match_m.fypp	\
   nad_trans_m.fypp null_rot_m.fypp num_par_m.fypp orbit_par_m.fypp	\
   order_m.fypp osc_file_m.fypp osc_par_m.fypp out_par_m.fypp		\
   out_util_m.fypp parfait_file_m.fypp parfait_model_m.fypp		\
   parfaitd_file_m.fypp parfaitd_model_m.fypp point_m.fypp		\
   poly_file_m.fypp poly_model_m.fypp prox_search_m.fypp		\
   qad_eval_m.fypp rad_bound_m.fypp rad_bvp_m.fypp rad_diff_m.fypp	\
   rad_eqns_m.fypp rad_match_m.fypp rad_trans_m.fypp resp_m.fypp	\
   root_m.fypp rot_factory_m.fypp rot_m.fypp rot_par_m.fypp		\
   sad_bound_m.fypp sad_bvp_m.fypp sad_diff_m.fypp sad_eqns_m.fypp	\
   sad_match_m.fypp sad_trans_m.fypp scan_m.fypp scan_par_m.fypp	\
   search_m.fypp state_m.fypp status_m.fypp string_m.fypp		\
   summary_m.fypp sysmtx_factory_m.fypp sysmtx_m.fypp system_m.fypp	\
   tar_eigen_m.fypp tar_fit_m.fypp tar_rot_m.fypp tidal_coeff_m.fypp	\
   tidal_resp_m.fypp tidal_context_m.fypp tide_par_m.fypp		\
   tnad_bvp_m.fypp tnad_diff_m.fypp tnad_eqns_m.fypp			\
   trapz_diff_m.fypp txt_writer_m.fypp util_m.fypp version_m.fypp	\
   wave_m.fypp wdec_file_m.fypp writer_m.fypp

# Directories

SRC_DIRS = $(addprefix ../src/, ad bvp cad common context core diff	\
   ext extern frontend/gyre frontend/tools grid interp lib math		\
   math/intrinsic matrix mode model nad output par poly rad rot sad	\
   search tar tide tnad)

# Include standard rules

include Make.inc
