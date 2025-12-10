#!/bin/bash
env -i PATH=/lustre09/project/6019267/shared/tools/variants/annotation/openCravat_env/bin:/usr/bin:/bin
module --force purge
module load StdEnv/2023 python/3.10
source /lustre09/project/6019267/shared/tools/variants/annotation/openCravat_env/bin/activate

# Run openCravat
{0}
