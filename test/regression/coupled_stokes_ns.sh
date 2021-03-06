#!/bin/bash

PROG="${GRINS_TEST_DIR}/grins_flow_regression"

INPUT="${GRINS_TEST_INPUT_DIR}/coupled_stokes_ns.in ${GRINS_TEST_DATA_DIR}/coupled_stokes_ns.xdr 1.0e-11"

PETSC_OPTIONS="-pc_factor_levels 4 -sub_pc_factor_levels 4"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS 
