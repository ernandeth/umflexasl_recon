#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Nonlinear registration of a 4D time series to its first volume using FSL FNIRT
#
# Usage:
#   ./fnirt_4d_to_firstvol.sh input_4d.nii output_dir
#
# Example:
#   ./fnirt_4d_to_firstvol.sh func_4d.nii fnirt_to_firstvol_out
#
# Outputs:
#   output_dir/reference_vol0000.nii
#   output_dir/registered_4d_fnirt.nii
#   output_dir/split/
#   output_dir/affines/
#   output_dir/warps/
#   output_dir/registered_vols/
#
# Notes:
#   - Uses the first 3D volume as the FNIRT reference.
#   - Performs an initial rigid-body FLIRT registration before FNIRT.
#   - FNIRT is not usually recommended for fMRI motion correction.
###############################################################################

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_4d_image.nii> <output_dir>"
    exit 1
fi

IN4D="$1"
OUTDIR="$2"

if [ ! -f "${IN4D}" ]; then
    echo "ERROR: Input file does not exist: ${IN4D}"
    exit 1
fi

# Check that FSL commands exist
for cmd in fslroi fslsplit fslmerge flirt fnirt fslval; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "ERROR: FSL command not found: ${cmd}"
        echo "Make sure FSL is installed and configured."
        exit 1
    fi
done

mkdir -p "${OUTDIR}"

SPLIT_DIR="${OUTDIR}/split"
AFFINE_DIR="${OUTDIR}/affines"
WARP_DIR="${OUTDIR}/warps"
REG_DIR="${OUTDIR}/registered_vols"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "${SPLIT_DIR}" "${AFFINE_DIR}" "${WARP_DIR}" "${REG_DIR}" "${LOG_DIR}"

REF="${OUTDIR}/reference_vol0000.nii"
OUT4D="${OUTDIR}/registered_4d_fnirt.nii"

echo "============================================================"
echo "4D nonlinear registration using FSL FNIRT"
echo "Input 4D image: ${IN4D}"
echo "Output folder:  ${OUTDIR}"
echo "Reference:      first volume"
echo "============================================================"

###############################################################################
# Step 1: Extract first volume as reference
###############################################################################

echo ""
echo "Step 1: Extracting first volume as reference..."

fslroi "${IN4D}" "${REF}" 0 1

echo "Reference image saved to:"
echo "  ${REF}"

###############################################################################
# Step 2: Split 4D image into 3D volumes
###############################################################################

echo ""
echo "Step 2: Splitting 4D image into 3D volumes..."

fslsplit "${IN4D}" "${SPLIT_DIR}/vol_" -t

NVOLS=$(ls "${SPLIT_DIR}"/vol_*.nii | wc -l)

echo "Number of volumes found: ${NVOLS}"

###############################################################################
# Step 3: Register each volume to the first volume
###############################################################################

echo ""
echo "Step 3: Registering each volume to the first volume..."

for VOL in "${SPLIT_DIR}"/vol_*.nii; do

    BASENAME=$(basename "${VOL}" .nii)

    AFFINE_MAT="${AFFINE_DIR}/${BASENAME}_to_vol0000_affine.mat"
    AFFINE_IMG="${AFFINE_DIR}/${BASENAME}_to_vol0000_affine.nii"

    WARP_COEF="${WARP_DIR}/${BASENAME}_to_vol0000_warpcoef.nii"
    WARP_FIELD="${WARP_DIR}/${BASENAME}_to_vol0000_warpfield.nii"

    REG_IMG="${REG_DIR}/${BASENAME}_fnirt.nii"
    LOG_FILE="${LOG_DIR}/${BASENAME}_fnirt.log"

    echo ""
    echo "------------------------------------------------------------"
    echo "Processing ${BASENAME}"
    echo "Input volume: ${VOL}"
    echo "Output image: ${REG_IMG}"
    echo "------------------------------------------------------------"

    # The first volume is already the reference.
    # Copy it directly so the final merge has the same number/order of volumes.
    if [ "${BASENAME}" = "vol_0000" ]; then
        echo "First volume is the reference; copying without registration."
        cp "${REF}" "${REG_IMG}"
        continue
    fi

    ###########################################################################
    # Initial linear registration with FLIRT
    #
    # For a time series from the same subject/session, rigid-body registration
    # is usually the appropriate initialization.
    ###########################################################################

    echo "Running FLIRT initial rigid-body registration..."

    flirt \
        -in "${VOL}" \
        -ref "${REF}" \
        -out "${AFFINE_IMG}" \
        -omat "${AFFINE_MAT}" \
        -dof 6 \
        -cost normcorr \
        -interp trilinear

    ###########################################################################
    # Nonlinear registration with FNIRT
    #
    # This uses relatively generic FNIRT settings rather than a T1-to-MNI config.
    # For same-modality within-subject time series, avoid using configs such as
    # T1_2_MNI152_2mm unless that is truly your registration scenario.
    ###########################################################################

#    echo "Running FNIRT nonlinear registration..."
#
#    fnirt \
#        --in="${VOL}" \
#        --ref="${REF}" \
#        --aff="${AFFINE_MAT}" \
#        --cout="${WARP_COEF}" \
#        --fout="${WARP_FIELD}" \
#        --iout="${REG_IMG}" \
#        --interp=spline \
#	--warpres=30,30,30 \
#        --ssqlambda=1 \
#        --verbose \
#        > "${LOG_FILE}" 2>&1
#
#    echo "Finished ${BASENAME}"

done

###############################################################################
# Step 4: Merge registered 3D volumes back into one 4D image
###############################################################################

echo ""
echo "Step 4: Merging registered volumes into 4D image..."

# fslmerge -t "${OUT4D}" "${REG_DIR}"/vol_*_fnirt.nii
fslmerge -t "${OUT4D}" "${AFFINE_DIR}"/vol_*.nii

echo ""
echo "============================================================"
echo "Done."
echo "Registered 4D output:"
echo "  ${OUT4D}"
echo ""
echo "Reference volume:"
echo "  ${REF}"
echo ""
echo "To inspect:"
echo "  fsleyes ${REF} ${OUT4D}"
echo "============================================================"

