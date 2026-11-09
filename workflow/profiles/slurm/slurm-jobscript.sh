#!/bin/bash
# properties = {properties}
# Keep the scheduler GPU assignment when containers use --cleanenv.
if [[ "${{CUDA_VISIBLE_DEVICES+x}}" == x ]]; then
    export SINGULARITYENV_CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES"
    export APPTAINERENV_CUDA_VISIBLE_DEVICES="$CUDA_VISIBLE_DEVICES"
fi
{exec_job}
