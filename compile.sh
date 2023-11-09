cmake .. \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ARCH_TURING75=ON \
  -DKokkos_ENABLE_CUDA_LAMBDA=ON \
  -DCMAKE_BUILD_TYPE=Release 
