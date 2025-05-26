#pragma once

//***************************************************************************
//* Based on the wrapper from the SPAdes project:
//* https://github.com/ablab/spades/blob/main/src/common/utils/parallel/openmp_wrapper.h
//***************************************************************************

#ifndef __OMP_WRAPPER_H__
#define __OMP_WRAPPER_H__

#include <cstdlib>

#ifdef _OPENMP
# include <omp.h>
#else
/* Provide single-threaded stubs */
# define omp_set_num_threads(x)  ((void)(x))
# define omp_get_max_threads()   1
# define omp_get_thread_num()    0
# define omp_get_num_threads()   1
# define omp_lock_t              size_t
# define omp_init_lock(x)        ((void)(x))
# define omp_destroy_lock(x)     ((void)(x))
# define omp_set_lock(x)         ((void)(x))
# define omp_unset_lock(x)       ((void)(x))
#endif

static inline unsigned spades_set_omp_threads(unsigned max_threads) {
    // Fix number of threads according to OMP capabilities.
    unsigned omp_threads = (unsigned)omp_get_max_threads();
    if (max_threads > omp_threads)
        max_threads = omp_threads;

    // Inform OpenMP runtime about this :)
    omp_set_num_threads((int) max_threads);

    return max_threads;
}

#endif /* __OMP_WRAPPER_H__ */