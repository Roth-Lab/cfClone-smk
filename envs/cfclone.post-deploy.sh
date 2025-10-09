#!env bash
set -o pipefail
set -e


export LD_PRELOAD=${CONDA_PREFIX}/lib/libstdc++.so.6
export JULIAUP_DEPOT_PATH_BACKUP=${JULIAUP_DEPOT_PATH:-}
export JULIAUP_DEPOT_PATH="$CONDA_PREFIX/depots/juliaup"
export JULIA_DEPOT_PATH_BACKUP=${JULIA_DEPOT_PATH:-}
export JULIA_DEPOT_PATH="$CONDA_PREFIX/depots/julia"
export BRIDGESTAN="$CONDA_PREFIX/bin/bridgestan"

mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d
mkdir -p ${CONDA_PREFIX}/etc/conda/deactivate.d
mkdir -p ${CONDA_PREFIX}/depots/juliaup

activate_fp=${CONDA_PREFIX}/etc/conda/activate.d/cfclone_activate.sh
deactivate_fp=${CONDA_PREFIX}/etc/conda/deactivate.d/cfclone_deactivate.sh


mkdir -p ${CONDA_PREFIX}/bin/bridgestan
cd ${CONDA_PREFIX}/bin

git clone --single-branch --depth 1 ssh://git@github.com/Roth-Lab/pt-cfclone cfclone

rm -rf ./cfclone/.git

cd cfclone

sed -i'.bak' 's/n_local_mpi_processes = n_threads/n_threads = n_threads/' lib.jl
sed -i 's/target,/target, multithreaded = true,/' lib.jl

juliaup add 1.11.5

juliaup default 1.11.5

julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'

pkg_version=$(julia -e 'using Pkg; Pkg.activate("."); include("lib.jl"); println(BridgeStan.get_version())')

cd ${CONDA_PREFIX}/bin/bridgestan

wget -O bridgestan.tar.gz "https://github.com/roualdes/bridgestan/releases/download/v${pkg_version}/bridgestan-${pkg_version}.tar.gz"

tar -xzf bridgestan.tar.gz --strip-components=1

rm bridgestan.tar.gz

mkdir -p make

echo "TBB_CXX_TYPE=gcc"  >> make/local
echo "TBB_INTERFACE_NEW=true" >> make/local
echo "TBB_INC=${CONDA_PREFIX}/include/" >> make/local
echo "TBB_LIB=${CONDA_PREFIX}/lib/" >> make/local
echo "PRECOMPILED_HEADERS=false" >> make/local
echo "STAN_THREADS=true" >> make/local

cd ${CONDA_PREFIX}/bin/cfclone

julia -e 'using Pkg; Pkg.activate("."); include("lib.jl"); infer(cfdna_tsv = "test/data/VOA10055P_cfdna.tsv", clone_cn_tsv = "test/data/VOA10055P_clone_cn.tsv", n_threads = 1, n_rounds=1, subsampling = 1000)'

rm -r ./results

echo 'export PATH="$CONDA_PREFIX/bin/cfclone:$PATH"' > $activate_fp
echo 'export JULIAUP_DEPOT_PATH_BACKUP="${JULIAUP_DEPOT_PATH:-}"' >> $activate_fp
echo 'export JULIAUP_DEPOT_PATH="$CONDA_PREFIX/depots/juliaup"' >> $activate_fp
echo 'export JULIA_DEPOT_PATH_BACKUP="${JULIA_DEPOT_PATH:-}"' >> $activate_fp
echo 'export JULIA_DEPOT_PATH="$CONDA_PREFIX/depots/julia"' >> $activate_fp
echo 'export BRIDGESTAN_BACKUP="${BRIDGESTAN:-}"' >> $activate_fp
echo 'export BRIDGESTAN="$CONDA_PREFIX/bin/bridgestan"' >> $activate_fp
echo 'export LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6"' >> $activate_fp

echo 'export JULIAUP_DEPOT_PATH="${JULIAUP_DEPOT_PATH_BACKUP}"' > $deactivate_fp
echo "unset JULIAUP_DEPOT_PATH_BACKUP" >> $deactivate_fp
echo 'export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_BACKUP}"' >> $deactivate_fp
echo "unset JULIA_DEPOT_PATH_BACKUP" >> $deactivate_fp
echo 'export BRIDGESTAN="${BRIDGESTAN_BACKUP}"' >> $deactivate_fp
echo "unset BRIDGESTAN_BACKUP" >> $deactivate_fp
echo "unset LD_PRELOAD" >> $deactivate_fp