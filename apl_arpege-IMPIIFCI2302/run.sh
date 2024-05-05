#!/bin/bash

. /home/gmap/mrpm/khatib/public/bin/env_intel-mpi-2023.2

set -x




# -g -O2 -march=core-avx2 -finline-functions -finline-limit=500 -Winline -qopt-prefetch=4 -fast-transcendentals -fimf-use-svml -no-fma \
# -convert big_endian -assume byterecl -align array64byte,all -traceback \
# -fpic -qopenmp -qopenmp-threadprivate compat -fp-model source \
# -ftz -diag-disable=remark,cpu-dispatch \
# -qopt-report=5 -qopt-report-phase=vec -free -DLINUX -DLITTLE_ENDIAN -DLITTLE -DADDRESS64 \
# -list -show nomap \

exec /home/gmap/mrpm/khatib/public/bin/mpiifort_wrapper \
  -c \
  -g -O0 \
  apl_arpege.F90
