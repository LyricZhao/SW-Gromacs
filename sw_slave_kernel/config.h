/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Include file for configuration macros from the build system.
 *
 * This header is not installed, so headers must not reference macros defined
 * here.
 *
 * \inlibraryapi
 */
#ifndef GMX_CONFIG_H
#define GMX_CONFIG_H

/* TODO: For now, disable Doxygen warnings from here */
/*! \cond */

/* IEEE754 floating-point format. Memory layout is defined by macros
 * GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER and GMX_IEEE754_BIG_ENDIAN_WORD_ORDER. 
 */
#define GMX_FLOAT_FORMAT_IEEE754

/* Work around broken calloc() */
/* #undef GMX_BROKEN_CALLOC */

/* Do not optimize FFTW setups (not needed with SSE FFT kernels) */
/* #undef GMX_DISABLE_FFTW_MEASURE */

/* Use Built-in FFTPACK FFT library */
#define GMX_FFT_FFTPACK

/* Use FFTW3 FFT library */
/* #undef GMX_FFT_FFTW3 */

/* Use Intel MKL FFT library */
/* #undef GMX_FFT_MKL */

/* Target platform is x86 or x86_64 */
/* #undef GMX_TARGET_X86 */

/* Target platform is BlueGene/Q */
/* #undef GMX_TARGET_BGQ */

/** Define if we are building natively on Windows */
/* #undef GMX_NATIVE_WINDOWS */

/** Define if we are building for Cygwin */
/* #undef GMX_CYGWIN */

/** Define if we have sufficient C++11 support */
/* #undef GMX_CXX11 */

/* GCC bug in AVX maskload/maskstore arguments - worked around internally */
/* #undef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG */

/* SSE2 was selected for SIMD instruction set level */
/* #undef GMX_SIMD_X86_SSE2 */

/* SSE4.1 was selected as SIMD instructions */
/* #undef GMX_SIMD_X86_SSE4_1 */

/* AVX 128-bit FMA was selected as SIMD instructions */
/* #undef GMX_SIMD_X86_AVX_128_FMA */

/* AVX 256-bit was selected as SIMD instructions */
/* #undef GMX_SIMD_X86_AVX_256 */

/* AVX2 256-bit SIMD instruction set level was selected */
/* #undef GMX_SIMD_X86_AVX2_256 */

/* MIC (Xeon Phi) SIMD instruction set level was selected */
/* #undef GMX_SIMD_X86_MIC */

/* AVX-512F foundation level instruction SIMD */
/* #undef GMX_SIMD_X86_AVX_512F */

/* AVX-512ER foundation level instruction SIMD */
/* #undef GMX_SIMD_X86_AVX_512ER */

/* 32-bit ARM NEON SIMD instruction set level was selected */
/* #undef GMX_SIMD_ARM_NEON */

/* ARM (AArch64) NEON Advanced SIMD instruction set level was selected */
/* #undef GMX_SIMD_ARM_NEON_ASIMD */

/* IBM QPX was selected as SIMD instructions (e.g. BlueGene/Q) */
/* #undef GMX_SIMD_IBM_QPX */

/* IBM VMX was selected as SIMD instructions (Power 6 and later) */
/* #undef GMX_SIMD_IBM_VMX */

/* IBM VSX was selected as SIMD instructions (Power 7 and later) */
/* #undef GMX_SIMD_IBM_VSX */
 
/* Fujitsu Sparc64 HPC-ACE SIMD acceleration */
/* #undef GMX_SIMD_SPARC64_HPC_ACE */

/* Reference SIMD implementation for testing */
/* #undef GMX_SIMD_REFERENCE */

/* String for SIMD instruction choice (for writing to log files and stdout) */
#define GMX_SIMD_STRING "NONE"

/* Calling convention string (if any) for routines with SIMD variable args */
#define gmx_simdcall 

/* Target mantissa accuracy for SIMD single precision math */
#define GMX_SIMD_ACCURACY_BITS_SINGLE 22

/* Target mantissa accuracy for SIMD double precision math */
#define GMX_SIMD_ACCURACY_BITS_DOUBLE 44

/* Integer byte order is big endian. */
/* #undef GMX_INTEGER_BIG_ENDIAN */

/* Use our own instead of system XDR libraries */
/* #undef GMX_INTERNAL_XDR */

/* Compile to use TNG library */
#define GMX_USE_TNG

/* Add support for tracing using Extrae */
/* #undef HAVE_EXTRAE */

/* Use MPI (with mpicc) for parallelization */
#define GMX_LIB_MPI

/* Use threads_mpi for parallelization */
/* #undef GMX_THREAD_MPI */

#if defined GMX_LIB_MPI || defined GMX_THREAD_MPI
/* Make a parallel version of GROMACS using message passing
   (MPI or thread_mpi) */
#define GMX_MPI
#endif

/* MPI_IN_PLACE exists for collective operations */
#define MPI_IN_PLACE_EXISTS

/* Use OpenMP multithreading */
/* #undef GMX_OPENMP */

/* Can and should use nice(3) to set priority */
#define GMX_USE_NICE

/* Maximum number of OpenMP threads supported */
#define GMX_OPENMP_MAX_THREADS 32

/* Use if can't rename checkpoints */
/* #undef GMX_NO_RENAME */

/* Use (modified) Gamess-UK for QM-MM calculations */
/* #undef GMX_QMMM_GAMESS */

/* Use (modified) Gaussian0x for QM-MM calculations */
/* #undef GMX_QMMM_GAUSSIAN */

/* Use (modified) Mopac 7 for QM-MM calculations */
/* #undef GMX_QMMM_MOPAC */

/* Use ORCA for QM-MM calculations */
/* #undef GMX_QMMM_ORCA */

/* Use the GROMACS software 1/sqrt(x) */
#define GMX_SOFTWARE_INVSQRT

/* Use sub-counters */
/* #undef GMX_CYCLE_SUBCOUNTERS */

/* Compile with plugin support */
/* #undef GMX_USE_PLUGINS */

/* Fallback path for VMD plug-ins */
#define GMX_VMD_PLUGIN_PATH ""

/* Define when pthreads are used */
#define THREAD_PTHREADS

/* Define when Windows threads are used */
/* #undef THREAD_WINDOWS */

/* Define native atomic operations are found */
#define TMPI_ATOMICS

/* Define for busy wait option  */
/* See gmxpre-config.h.cmakein for explanation for the #ifdef */
#ifndef TMPI_WAIT_FOR_NO_ONE
/* #undef TMPI_WAIT_FOR_NO_ONE */
#endif

/* Define for copy buffer option */
/* #undef TMPI_COPY_BUFFER */

/* Define for tmpi warnings option */
/* #undef TMPI_WARNINGS */

/* Define for profiling option */
/* #undef TMPI_PROFILE */

/* Define for Linux pthread_setaffinity_np */
#define HAVE_PTHREAD_SETAFFINITY

/* Define for X-Windows */
/* #undef GMX_X11 */

/* Enable x86 gcc inline assembly */
/* #undef GMX_X86_GCC_INLINE_ASM */

/* Use GPU native acceleration */
/* #undef GMX_GPU */

/* CUDA runtime API version (identical to CUDART_VERSION from cuda_runtime_api.h) */
/* #undef GMX_CUDA_VERSION */

/* Defined if texture objects are supported by CUDA */
/* #undef HAVE_CUDA_TEXOBJ_SUPPORT */

/* Use NVML */
/* #undef HAVE_NVML */

/* Use OpenCL acceleators */
/* #undef GMX_USE_OPENCL */

/* Define relative path to OpenCL kernels */
#define OCL_INSTALL_DIR "share/gromacs/opencl"

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO

/* Define to 1 if _fseeki64 (and presumably _fseeki64) exists and is declared. */
/* #undef HAVE__FSEEKI64 */

/* Have io.h (windows)*/
/* #undef HAVE_IO_H */

/* Define to 1 if you have the posix_memalign() function. */
#define HAVE_POSIX_MEMALIGN

/* Define to 1 if you have the memalign() function. */
/* #undef HAVE_MEMALIGN */

/* Define to 1 if you have the MSVC _aligned_malloc() function. */
/* #undef HAVE__ALIGNED_MALLOC */

/* Define to 1 if you have the clock_gettime() function. */
#define HAVE_CLOCK_GETTIME

/* Define to 1 if you have the gettimeofday() function. */
#define HAVE_GETTIMEOFDAY

/* Define to 1 if you have the rdtscp instruction. */
/* #undef HAVE_RDTSCP */

/* Define to 1 if you have the isfinite() function. */
/* #undef HAVE_ISFINITE */

/* Define to 1 if you have the _isfinite() function. */
/* #undef HAVE__ISFINITE */

/* Define to 1 if you have the _finite() function. */
/* #undef HAVE__FINITE */

/* Define to 1 if you have the fsync() function. */
#define HAVE_FSYNC

/* Define to 1 if you have the Windows _commit() function. */
/* #undef HAVE__COMMIT */

/* Define to 1 if you have the fileno() function. */
#define HAVE_FILENO

/* Define to 1 if you have the _fileno() function. */
/* #undef HAVE__FILENO */

/* Define to 1 if you have the sigaction() function. */
#define HAVE_SIGACTION

/* Define to 1 if you have the rsqrt() function. */
/* #undef HAVE_RSQRT */

/* Define to 1 if you have the rsqrtf() function. */
/* #undef HAVE_RSQRTF */

/* Define to 1 if you have the sqrtf() function. */
#define HAVE_SQRTF

/* Define to 1 if yo have the <unistd.h> header file. */
#define HAVE_UNISTD_H

/* Define to 1 if yo have the <pwd.h> header file. */
#define HAVE_PWD_H

/* Define to 1 if yo have the <dirent.h> header file. */
#define HAVE_DIRENT_H

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H

/* Define to 1 if you have the <x86intrin.h> header file */
/* #undef HAVE_X86INTRIN_H */

/* Define to 1 if you have the <intrin.h> header file */
/* #undef HAVE_INTRIN_H */

/* Define to 1 if you have the <sched.h> header */
#define HAVE_SCHED_H

/* Define to 1 if you have the POSIX <regex.h> header file. */
#define HAVE_POSIX_REGEX

/* Define to 1 if you have the C++11 <regex> header file. */
/* #undef HAVE_CXX11_REGEX */

/* Define to 1 if you have the sysconf() function */
#define HAVE_SYSCONF

/* Define to 1 if you have the all the affinity functions in sched.h */
/* #undef HAVE_SCHED_AFFINITY */

/* Bytes in IEEE fp word are in big-endian order if set, little-endian if not.
   Only relevant when FLOAT_FORMAT_IEEE754 is defined. */
/* #undef GMX_IEEE754_BIG_ENDIAN_BYTE_ORDER */

/* The two words in a double precision variable are in b ig-endian order if
   set, little-endian if not. Do NOT assume this is the same as the byte
   order! Only relevant when FLOAT_FORMAT_IEEE754 is defined. */
/* #undef GMX_IEEE754_BIG_ENDIAN_WORD_ORDER */

/* Define if SIGUSR1 is present */
#define HAVE_SIGUSR1

/* Enable gromacs quotes */
#define GMX_COOL_QUOTES

/* default name mangling maybe wrong on exotic plattforms */
#define F77_FUNC(name,NAME) name ## _

/* Define if we have pipes */
#define HAVE_PIPES

/* Define if we have feenableexcept */
#define HAVE_FEENABLEEXCEPT

/* Define if we have zlib */
#define HAVE_ZLIB

/*! \endcond */

#endif
