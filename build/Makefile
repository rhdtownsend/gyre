# File     : Makefile
# Purpose  : makefile for build subdir

# Default rule in case make is invoked in current dir

all:
	make -C ..

# Targets

LIB_TARGETS = libgyre

INC_FILES = gyre_m.mod gyre.inc

ifeq ($(FRONTENDS),yes)
   BIN_TARGETS += gyre gyre_tides
endif

ifeq ($(TOOLS),yes)
   BIN_TARGETS += eval_lambda eval_tidal_coeff eval_love build_poly poly_to_fgong poly_to_txt build_tar_fit
endif

ifeq ($(IFACES),yes)
   BIN_TARGETS += test_gyre_mesa
   LIB_TARGETS += libgyre_mesa
   INC_FILES += gyre_mesa_m.mod
endif

# Source file lists & related variables
#
# xxx_SRCS = source files required to compile xxx
# xxx_USES = source files providing modules USE'd by xxx
# xxx_LIBS = library files required to link xxx

gyre_SRCS = gyre.fypp
gyre_USES = $(libgyre_SRCS)
gyre_LIBS = libgyre

gyre_tides_SRCS = gyre_tides.fypp
gyre_tides_USES = $(libgyre_SRCS)
gyre_tides_LIBS = libgyre

eval_lambda_SRCS = eval_lambda.fypp
eval_lambda_USES = $(libgyre_SRCS)
eval_lambda_LIBS = libgyre

eval_tidal_coeff_SRCS = eval_tidal_coeff.fypp
eval_tidal_coeff_USES = $(libgyre_SRCS)
eval_tidal_coeff_LIBS = libgyre

build_poly_SRCS = build_poly.fypp lane_emden_m.fypp
build_poly_USES = $(libgyre_SRCS)
build_poly_LIBS = libgyre

poly_to_fgong_SRCS = poly_to_fgong.fypp
poly_to_fgong_USES = $(libgyre_SRCS)
poly_to_fgong_LIBS = libgyre

poly_to_txt_SRCS = poly_to_txt.fypp
poly_to_txt_USES = $(libgyre_SRCS)
poly_to_txt_LIBS = libgyre

build_tar_fit_SRCS = build_tar_fit.fypp
build_tar_fit_USES = $(libgyre_SRCS)
build_tar_fit_LIBS = libgyre

eval_love_SRCS = eval_love.fypp
eval_love_USES = $(libgyre_SRCS)
eval_love_LIBS = libgyre

test_gyre_mesa_SRCS = test_gyre_mesa.fypp
test_gyre_mesa_USES = $(libgyre_mesa_SRCS) $(libgyre_SRCS)
test_gyre_mesa_LIBS = libgyre_mesa libgyre

libgyre_SRCS = ad_bound_m.fypp ad_bvp_m.fypp ad_diff_m.fypp		\
   ad_eqns_m.fypp ad_match_m.fypp ad_trans_m.fypp amdl_file_m.fypp	\
   anapoly_0_model_m.fypp anapoly_1_model_m.fypp			\
   anapoly_5_model_m.fypp anapoly_5_1_model_m.fypp atmos_m.fypp		\
   b3_file_m.fypp band_sysmtx_m.fypp blas_lapack_m.fypp			\
   block_sysmtx_m.fypp bound_m.fypp bracket_search_m.fypp bvp_m.fypp	\
   cheb_fit_m.fypp colloc_diff_m.fypp constants_m.fypp context_m.fypp	\
   contour_map_m.fypp contour_path_m.fypp contour_search_m.fypp		\
   cubic_interp_m.fypp deriv_m.fypp detail_m.fypp diff_factory_m.fypp	\
   diff_m.fypp discrim_m.fypp eqns_m.fypp evol_model_m.fypp		\
   ext_complex_m.fypp ext_m.fypp ext_real_m.fypp famdl_file_m.fypp	\
   fgong_file_m.fypp freq_context_m.fypp freq_m.fypp gyre_m.fypp	\
   grid_m.fypp grid_par_m.fypp grid_refine_m.fypp			\
   grid_scaffold_m.fypp grid_spec_m.fypp grid_util_m.fypp		\
   grid_weights_m.fypp gsm_file_m.fypp hdf_writer_m.fypp		\
   hom_model_m.fypp interp_m.fypp integ_m.fypp linalg_m.fypp		\
   linear_interp_m.fypp losc_file_m.fypp math_m.fypp			\
   magnus_diff_m.fypp mesa_file_m.fypp minim_m.fypp			\
   minmod_search_m.fypp mirk_diff_m.fypp mode_m.fypp mode_par_m.fypp	\
   model_factory_m.fypp model_m.fypp model_par_m.fypp			\
   model_util_m.fypp nad_bound_m.fypp nad_bvp_m.fypp nad_diff_m.fypp	\
   nad_eqns_m.fypp nad_match_m.fypp nad_trans_m.fypp null_rot_m.fypp	\
   num_par_m.fypp orbit_par_m.fypp osc_file_m.fypp osc_par_m.fypp	\
   out_par_m.fypp out_util_m.fypp parfait_file_m.fypp			\
   parfait_model_m.fypp parfaitd_file_m.fypp parfaitd_model_m.fypp	\
   point_m.fypp poly_file_m.fypp poly_model_m.fypp prox_search_m.fypp	\
   qad_eval_m.fypp rad_bound_m.fypp rad_bvp_m.fypp rad_diff_m.fypp	\
   rad_eqns_m.fypp rad_match_m.fypp rad_trans_m.fypp resp_m.fypp	\
   root_m.fypp rot_factory_m.fypp rot_m.fypp rot_par_m.fypp		\
   sad_bound_m.fypp sad_bvp_m.fypp sad_diff_m.fypp sad_eqns_m.fypp	\
   sad_match_m.fypp sad_trans_m.fypp scan_m.fypp scan_par_m.fypp	\
   search_m.fypp special_m.fypp state_m.fypp status_m.fypp		\
   summary_m.fypp sysmtx_factory_m.fypp sysmtx_m.fypp			\
   tar_eigen_m.fypp tar_fit_m.fypp tar_rot_m.fypp tidal_coeff_m.fypp	\
   tidal_resp_m.fypp tidal_context_m.fypp tide_par_m.fypp		\
   tnad_bvp_m.fypp tnad_diff_m.fypp tnad_eqns_m.fypp			\
   trapz_diff_m.fypp txt_writer_m.fypp util_m.fypp version_m.fypp	\
   wave_m.fypp wdec_file_m.fypp writer_m.fypp

libgyre_mesa_SRCS = gyre_mesa_m.fypp
libgyre_mesa_USES = $(libgyre_SRCS)
libgyre_mesa_LIBS = libgyre

# Build flags

FFLAGS = -I$(FORUM_INC_DIR)

FPPFLAGS = -I$(FORUM_INC_DIR)

LDFLAGS = `lapack_link` `hdf5_link` `lapack95_link` `odepack_link` -L$(FORUM_LIB_DIR) -lforum

# Include standard rules

include Make.inc

# Target for linting

linter :
	fortitude check --file-extensions=f90

# Targets for source checking

check_src:
	@for src in $(ALL_SRCS); do \
             find $(SRC_DIRS) -path '*/forum/*' -prune -or -name $$src -print | grep -q  . || echo $$src MISSING; \
         done

check_src : SRC_DIRS ?= ../src

.PHONY: all linter check_src
