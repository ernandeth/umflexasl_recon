#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# Estimate FLIRT affine transformations from a 4D magnitude time series,
# then apply those transformations to a corresponding 4D phase time series.
#
# Usage:
#   ./flirt_realign_phase_from_magnitude.sh magnitude_4d.nii phase_4d.nii output_dir [dof] [phase_interp]
#
# Example:
#   ./flirt_realign_phase_from_magnitude.sh mag.nii phase.nii flirt_realign_out
#
# Example with 12 DOF and nearest-neighbor phase interpolation:
#   ./flirt_realign_phase_from_magnitude.sh mag.nii phase.nii flirt_realign_out 12 nearestneighbour
#
# Defaults:
#   dof=12
#   phase_interp=trilinear
#
# Main outputs:
#   output_dir/magnitude_flirt_realigned.nii
#   output_dir/phase_flirt_realigned.nii
#   output_dir/mats/mag_vol_XXXX_to_ref.mat
###############################################################################

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <magnitude_4d.nii> <phase_4d.nii> <output_dir> [dof] [phase_interp]"
    exit 1
fi

MAG4D="$1"
REAL4D="$2"
IMAG4D="$3"
OUTDIR="$4"
DOF="${5:-12}"
REAL_INTERP="${6:-trilinear}"

# Interpolation for magnitude registration output
MAG_INTERP="trilinear"

mkdir -p "${OUTDIR}"

MAG_SPLIT_DIR="${OUTDIR}/mag_split"
CPLX_SPLIT_DIR="${OUTDIR}/cplx_split"
MAG_REG_DIR="${OUTDIR}/mag_registered_vols"
CPLX_REG_DIR="${OUTDIR}/cplx_registered_vols"
MAT_DIR="${OUTDIR}/mats"
LOG_DIR="${OUTDIR}/logs"

mkdir -p \
    "${MAG_SPLIT_DIR}" \
    "${CPLX_SPLIT_DIR}" \
    "${MAG_REG_DIR}" \
    "${CPLX_REG_DIR}" \
    "${MAT_DIR}" \
    "${LOG_DIR}"

MAG_REF="${OUTDIR}/magnitude_ref_vol0000.nii"
REAL_REF="${OUTDIR}/real_ref_vol0000.nii"
IMAG_REF="${OUTDIR}/imag_ref_vol0000.nii"

MAG_OUT4D="${OUTDIR}/r_im_mag.nii"
REAL_OUT4D="${OUTDIR}/r_im_real.nii"
IMAG_OUT4D="${OUTDIR}/r_im_imag.nii"

echo "============================================================"
echo "FLIRT-based realignment of magnitude and phase time series"
echo "Magnitude 4D:       ${MAG4D}"
echo "Real 4D:            ${REAL4D}"
echo "Imaginary 4D:       ${REAL4D}"
echo "Output directory:   ${OUTDIR}"
echo "Reference volume:   first magnitude volume"
echo "FLIRT DOF:          ${DOF}"
echo "interpolation:${REAL_INTERP}"
echo "============================================================"

###############################################################################
# Check required commands
###############################################################################

for cmd in fslval fslroi fslsplit fslmerge flirt fslinfo; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "ERROR: FSL command not found: ${cmd}"
        echo "Make sure FSL is installed and configured."
        exit 1
    fi
done

#############################################################
# make sure we write nii and not compressed nii
###########################################################
export FSLOUTPUTTYPE=NIFTI 
echo "output is ...${FSLOUTPUTTYPE}"

###############################################################################
# Check input files
###############################################################################

if [ ! -f "${MAG4D}" ]; then
    echo "ERROR: Magnitude file not found: ${MAG4D}"
    exit 1
fi

if [ ! -f "${REAL4D}" ]; then
    echo "ERROR: real file not found: ${REAL4D}"
    exit 1
fi

if [ ! -f "${IMAG4D}" ]; then
    echo "ERROR:  Imaginary file not found: ${IMAG4D}"
    exit 1
fi


NMAG=$(fslval "${MAG4D}" dim4)
NREAL=$(fslval "${REAL4D}" dim4)
NIMAG=$(fslval "${IMAG4D}" dim4)

echo ""
echo "Number of magnitude volumes: ${NMAG}"
echo "Number of real volumes:     ${NREAL}"
echo "Number of imag volumes:     ${NIMAG}"

if [ "${NMAG}" -ne "${NREAL}" ]; then
    echo "ERROR: Magnitude and real images have different numbers of volumes."
    exit 1
fi

###############################################################################
# Step 1: Extract reference volumes
###############################################################################

echo ""
echo "Step 1: Extracting first magnitude and phase volumes..."

fslroi "${MAG4D}" "${MAG_REF}" 0 1
fslroi "${REAL4D}" "${REAL_REF}" 0 1
fslroi "${IMAG4D}" "${IMAG_REF}" 0 1

echo "Magnitude reference:"
echo "  ${MAG_REF}"
echo "Complex reference:"
echo "  ${REAL_REF}"
echo "  ${IMAG_REF}"

###############################################################################
# Step 2: Split 4D magnitude and phase images
###############################################################################

echo ""
echo "Step 2: Splitting 4D magnitude and phase time series..."

fslsplit "${MAG4D}" "${MAG_SPLIT_DIR}/mag_vol_" -t
fslsplit "${REAL4D}" "${CPLX_SPLIT_DIR}/real_vol_" -t
fslsplit "${IMAG4D}" "${CPLX_SPLIT_DIR}/imag_vol_" -t

###############################################################################
# Step 3: Estimate FLIRT transforms from magnitude and apply to phase
###############################################################################

echo ""
echo "Step 3: Estimating FLIRT transforms from magnitude and applying to phase..."

for ((i=0; i<${NMAG}; i++)); do

    IDX=$(printf "%04d" "${i}")

    MAG_VOL="${MAG_SPLIT_DIR}/mag_vol_${IDX}.nii"
    REAL_VOL="${CPLX_SPLIT_DIR}/real_vol_${IDX}.nii"
    IMAG_VOL="${CPLX_SPLIT_DIR}/imag_vol_${IDX}.nii"

    MAT="${MAT_DIR}/mag_vol_${IDX}_to_ref.mat"

    MAG_REG="${MAG_REG_DIR}/mag_vol_${IDX}_flirt.nii"
    REAL_REG="${CPLX_REG_DIR}/real_vol_${IDX}_flirt.nii"
    IMAG_REG="${CPLX_REG_DIR}/imag_vol_${IDX}_flirt.nii"

    MAG_LOG="${LOG_DIR}/mag_vol_${IDX}_flirt.log"
    REAL_LOG="${LOG_DIR}/real_vol_${IDX}_applyxfm.log"
    IMAG_LOG="${LOG_DIR}/imag_vol_${IDX}_applyxfm.log"

    echo "Processing volume ${IDX}"

    if [ ! -f "${MAG_VOL}" ]; then
        echo "ERROR: Missing magnitude volume: ${MAG_VOL}"
        exit 1
    fi

    if [ ! -f "${REAL_VOL}" ]; then
        echo "ERROR: Missing Real volume: ${REAL_VOL}"
        exit 1
    fi

    if [ ! -f "${IMAG_VOL}" ]; then
        echo "ERROR: Missing Imaginary volume: ${IMAG_VOL}"
        exit 1
    fi

    if [ "${IDX}" = "0000" ]; then

        echo "  Volume 0000 is the reference; using identity transform."

        cat > "${MAT}" <<EOF
1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1
EOF

        cp "${MAG_REF}" "${MAG_REG}"
        cp "${REAL_REF}" "${REAL_REG}"
        cp "${IMAG_REF}" "${IMAG_REG}"

    else

        #######################################################################
        # Estimate transform from current magnitude volume to reference magnitude
        #######################################################################

        flirt \
            -in "${MAG_VOL}" \
            -ref "${MAG_REF}" \
            -out "${MAG_REG}" \
            -omat "${MAT}" \
            -dof "${DOF}" \
            -cost normcorr \
            -interp "${MAG_INTERP}" \
            > "${MAG_LOG}" 2>&1

        #######################################################################
        # Apply same transform to corresponding phase volume
        #
        # This assumes:
        #   real_vol_XXXX is geometrically aligned with mag_vol_XXXX
        #   real_vol_0000 is geometrically aligned with mag_vol_0000
        #######################################################################

        flirt \
            -in "${IMAG_VOL}" \
            -ref "${IMAG_REF}" \
            -applyxfm \
            -init "${MAT}" \
            -out "${IMAG_REG}" \
            -interp "${REAL_INTERP}" \
            > "${IMAG_LOG}" 2>&1

        flirt \
            -in "${REAL_VOL}" \
            -ref "${REAL_REF}" \
            -applyxfm \
            -init "${MAT}" \
            -out "${REAL_REG}" \
            -interp "${REAL_INTERP}" \
            > "${REAL_LOG}" 2>&1

    fi

done

###############################################################################
# Step 4: Merge registered volumes back into 4D images
###############################################################################

echo ""
echo "Step 4: Merging registered volumes into 4D images..."

fslmerge -t "${MAG_OUT4D}" "${MAG_REG_DIR}"/mag_vol_*_flirt.nii
fslmerge -t "${REAL_OUT4D}" "${CPLX_REG_DIR}"/real_vol_*_flirt.nii
fslmerge -t "${IMAG_OUT4D}" "${CPLX_REG_DIR}"/imag_vol_*_flirt.nii

echo ""
echo "============================================================"
echo "Done."
echo ""
echo "Realigned magnitude 4D:"
echo "  ${MAG_OUT4D}"
echo ""
echo "Realigned complex 4D:"
echo "  ${REAL_OUT4D}"
echo "  ${IMAG_OUT4D}"
echo ""
echo "Affine matrices:"
echo "  ${MAT_DIR}/mag_vol_XXXX_to_ref.mat"
echo ""
echo "Magnitude reference:"
echo "  ${MAG_REF}"
echo ""
echo "complex reference:"
echo "  ${REAL_REF}"
echo "  ${IMAG_REF}"
echo ""
echo "To inspect:"
echo "  fsleyes ${MAG_REF} ${MAG_OUT4D} ${REAL_OUT4D} ${IMAG_OUT4D}"
echo ""
echo "To Clean up:"
echo "rm -rf realigned"
echo "============================================================"
