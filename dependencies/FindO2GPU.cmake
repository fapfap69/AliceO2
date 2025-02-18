# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

if(NOT DEFINED ENABLE_CUDA)
  set(ENABLE_CUDA "AUTO")
endif()
if(NOT DEFINED ENABLE_OPENCL1)
  set(ENABLE_OPENCL1 "AUTO")
endif()
if(NOT DEFINED ENABLE_OPENCL2)
  set(ENABLE_OPENCL2 "AUTO")
endif()
if(NOT DEFINED ENABLE_HIP)
  set(ENABLE_HIP "AUTO")
endif()
string(TOUPPER "${ENABLE_CUDA}" ENABLE_CUDA)
string(TOUPPER "${ENABLE_OPENCL1}" ENABLE_OPENCL1)
string(TOUPPER "${ENABLE_OPENCL2}" ENABLE_OPENCL2)
string(TOUPPER "${ENABLE_HIP}" ENABLE_HIP)

if(CUDA_COMPUTETARGET AND CUDA_COMPUTETARGET STREQUAL "default")
  set(CUDA_COMPUTETARGET 86 89)
endif()

if(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET STREQUAL "default")
  set(HIP_AMDGPUTARGET gfx906;gfx908)
endif()

function(set_target_cuda_arch target)
  if(CUDA_COMPUTETARGET AND CUDA_COMPUTETARGET STREQUAL "86")
    message(STATUS "Using optimized CUDA settings for Ampere GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_AMPERE)
  elseif(CUDA_COMPUTETARGET AND CUDA_COMPUTETARGET STREQUAL "75")
    message(STATUS "Using optimized CUDA settings for Turing GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_TURING)
  else()
    message(STATUS "Defaulting optimized CUDA settings for Ampere GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_AMPERE)
  endif()
endfunction()

function(set_target_hip_arch target)
  if(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET MATCHES "gfx906")
    message(STATUS "Using optimized HIP settings for MI50 GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_VEGA)
  elseif(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET MATCHES "gfx908")
    message(STATUS "Using optimized HIP settings for MI100 GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_MI2xx)
  elseif(HIP_AMDGPUTARGET AND HIP_AMDGPUTARGET MATCHES "gfx90a")
    message(STATUS "Using optimized HIP settings for MI210 GPU")
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_MI2xx)
  else()
    target_compile_definitions(${target} PUBLIC GPUCA_GPUTYPE_VEGA)
  endif()
endfunction()

# Detect and enable CUDA
STRING(REGEX REPLACE "\-std=[^ ]*" "" O2_GPU_CMAKE_CXX_FLAGS_NOSTD "${CMAKE_CXX_FLAGS}") # Need to strip c++17 imposed by alidist defaults

if(ENABLE_CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
  include(CheckLanguage)
  check_language(CUDA)
  if (NOT ENABLE_CUDA STREQUAL "AUTO")
    if (NOT CMAKE_CUDA_COMPILER)
      set(CMAKE_CUDA_COMPILER "nvcc")
    endif()
    set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler")
  endif()
  if(CMAKE_CUDA_COMPILER)
    if(GPUCA_CUDA_GCCBIN)
      message(STATUS "Using as CUDA GCC version: ${GPUCA_CUDA_GCCBIN}")
      set(CMAKE_CUDA_HOST_COMPILER "${GPUCA_CUDA_GCCBIN}")
    endif()
    if(CUDA_COMPUTETARGET)
      set(CMAKE_CUDA_ARCHITECTURES ${CUDA_COMPUTETARGET} CACHE STRING "" FORCE)
    else()
      set(CMAKE_CUDA_ARCHITECTURES 61-virtual CACHE STRING "" FORCE)
    endif()
    enable_language(CUDA)
    get_property(LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
    if (ENABLE_CUDA STREQUAL "AUTO")
      set(FAILURE_SEVERITY STATUS)
    else()
      set(FAILURE_SEVERITY FATAL_ERROR)
    endif()
    if(NOT CUDA IN_LIST LANGUAGES)
      message(${FAILURE_SEVERITY} "CUDA was found but cannot be enabled")
      set(CMAKE_CUDA_COMPILER OFF)
    endif()
    find_path(THRUST_INCLUDE_DIR thrust/version.h PATHS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} NO_DEFAULT_PATH)
    if(THRUST_INCLUDE_DIR STREQUAL "THRUST_INCLUDE_DIR-NOTFOUND")
      message(${FAILURE_SEVERITY} "CUDA found but thrust not available")
      set(CMAKE_CUDA_COMPILER OFF)
    endif()
    if (NOT CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "11.4")
      message(${FAILURE_SEVERITY} "CUDA Version too old: ${CMAKE_CUDA_COMPILER_VERSION}, 11.4 required")
      set(CMAKE_CUDA_COMPILER OFF)
    endif()
  endif()
  if(CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler \"${O2_GPU_CMAKE_CXX_FLAGS_NOSTD}\" --expt-relaxed-constexpr --extended-lambda --allow-unsupported-compiler -Xptxas -v")
    if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "12.3")
      string(APPEND CMAKE_CUDA_FLAGS " -Xcudafe --diag_suppress=20257")
    endif()
    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -lineinfo -Xcompiler \"${CMAKE_CXX_FLAGS_DEBUG}\" -Xptxas -O0 -Xcompiler -O0")
    if(NOT CMAKE_BUILD_TYPE STREQUAL "DEBUG")
      set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE}} -Xcompiler \"${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}\" -Xptxas -O4 -Xcompiler -O4")
    endif()
    if(DEFINED GPUCA_NO_FAST_MATH AND "${GPUCA_NO_FAST_MATH}")
      set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE}} --ftz=false --prec-div=true --prec-sqrt=true --fmad false")
    elseif(NOT CMAKE_BUILD_TYPE STREQUAL "DEBUG")
      set(CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_CUDA_FLAGS_${CMAKE_BUILD_TYPE}} -use_fast_math --ftz=true")#
    endif()
    if(CMAKE_CXX_FLAGS MATCHES "(^| )-Werror( |$)")
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Werror=cross-execution-space-call")
    endif()

    set(CUDA_ENABLED ON)
    message(STATUS "CUDA found (Version ${CMAKE_CUDA_COMPILER_VERSION})")
  elseif(NOT ENABLE_CUDA STREQUAL "AUTO")
    message(FATAL_ERROR "CUDA not found (Compiler: ${CMAKE_CUDA_COMPILER})")
  else()
    set(CUDA_ENABLED OFF)
  endif()
endif()

# Detect and enable OpenCL 1.2 from AMD
if(ENABLE_OPENCL1 OR ENABLE_OPENCL2)
  find_package(OpenCL)
  if((ENABLE_OPENCL1 AND NOT ENABLE_OPENCL1 STREQUAL "AUTO")
     OR (ENABLE_OPENCL2 AND NOT ENABLE_OPENCL2 STREQUAL "AUTO"))
    set_package_properties(OpenCL PROPERTIES TYPE REQUIRED)
  else()
    set_package_properties(OpenCL PROPERTIES TYPE OPTIONAL)
  endif()
endif()
if(ENABLE_OPENCL1)
  if(NOT AMDAPPSDKROOT)
    set(AMDAPPSDKROOT "$ENV{AMDAPPSDKROOT}")
  endif()

  if(OpenCL_FOUND
     AND OpenCL_VERSION_STRING VERSION_GREATER_EQUAL 1.2
     AND AMDAPPSDKROOT
     AND EXISTS "${AMDAPPSDKROOT}")
    set(OPENCL1_ENABLED ON)
    message(STATUS "Found AMD OpenCL 1.2")
  elseif(NOT ENABLE_OPENCL1 STREQUAL "AUTO")
    message(FATAL_ERROR "AMD OpenCL 1.2 not available")
  else()
    set(OPENCL1_ENABLED OFF)
  endif()
endif()

# Detect and enable OpenCL 2.x
if(ENABLE_OPENCL2)
  find_package(OpenCL)
  find_package(LLVM)
  if(LLVM_FOUND)
    find_package(Clang)
  endif()
  if(DEFINED ENV{ROCM_PATH})
    get_filename_component(ROCM_ROOT "$ENV{ROCM_PATH}" ABSOLUTE)
  else()
    set(ROCM_ROOT "/opt/rocm")
  endif()
  find_program(ROCM_AGENT_ENUMERATOR rocm_agent_enumerator PATHS "${ROCM_ROOT}/bin")
  if (GPUCA_OPENCL_CLANGBIN)
    set(LLVM_CLANG ${GPUCA_OPENCL_CLANGBIN})
    execute_process(COMMAND "which" "/usr/lib/llvm/15/bin/clang-15" OUTPUT_VARIABLE TMP_LLVM_SPIRV_PATH COMMAND_ERROR_IS_FATAL ANY)
    cmake_path(GET TMP_LLVM_SPIRV_PATH PARENT_PATH TMP_LLVM_SPIRV_PATH)
    find_program(LLVM_SPIRV llvm-spirv HINTS "${TMP_LLVM_SPIRV_PATH}")
  else()
    find_program(LLVM_CLANG clang HINTS "${Clang_DIR}/../../../bin-safe")
    find_program(LLVM_SPIRV llvm-spirv HINTS "${Clang_DIR}/../../../bin-safe")
  endif()
  if(Clang_FOUND
     AND LLVM_FOUND
     AND NOT LLVM_CLANG STREQUAL "LLVM_CLANG-NOTFOUND"
     AND LLVM_PACKAGE_VERSION VERSION_GREATER_EQUAL 13.0)
    set(OPENCL2_COMPATIBLE_CLANG_FOUND ON)
  endif()
  if(OpenCL_VERSION_STRING VERSION_GREATER_EQUAL 2.2
     AND NOT LLVM_SPIRV STREQUAL "LLVM_SPIRV-NOTFOUND"
     AND OPENCL2_COMPATIBLE_CLANG_FOUND)
    set(OPENCL2_ENABLED_SPIRV ON)
    message(STATUS "Using CLANG ${LLVM_CLANG} and ${LLVM_SPIRV} for SPIR-V compilation")
  endif ()
  if(OPENCL2_COMPATIBLE_CLANG_FOUND AND
     (OpenCL_VERSION_STRING VERSION_GREATER_EQUAL 2.2
     OR OPENCL2_ENABLED_SPIRV))
    set(OPENCL2_ENABLED ON)
    message(STATUS "Found OpenCL 2 (${OpenCL_VERSION_STRING} SPIR-V ${OPENCL2_ENABLED_SPIRV} with CLANG ${LLVM_PACKAGE_VERSION})")
  elseif(NOT ENABLE_OPENCL2 STREQUAL "AUTO")
    message(FATAL_ERROR "OpenCL 2.x not available")
  else()
    set(OPENCL2_ENABLED OFF)
  endif()
endif()

# Detect and enable HIP
if(ENABLE_HIP)
  set(CMAKE_HIP_STANDARD 17)
  set(CMAKE_HIP_STANDARD_REQUIRED TRUE)
  if(HIP_AMDGPUTARGET)
    set(AMDGPU_TARGETS "${HIP_AMDGPUTARGET}" CACHE STRING "AMD GPU targets to compile for" FORCE)
    set(GPU_TARGETS "${HIP_AMDGPUTARGET}" CACHE STRING "AMD GPU targets to compile for" FORCE)
    set(CMAKE_HIP_ARCHITECTURES "${HIP_AMDGPUTARGET}" CACHE STRING "AMD GPU targets to compile for" FORCE)
  endif()
  if(EXISTS "/opt/rocm/lib/cmake" AND (NOT DEFINED CMAKE_PREFIX_PATH OR NOT ${CMAKE_PREFIX_PATH} MATCHES "rocm") AND (NOT ENV{CMAKE_PREFIX_PATH} MATCHES "rocm"))
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/opt/rocm/lib/cmake")
    if (NOT DEFINED CMAKE_HIP_COMPILER)
      set(CMAKE_HIP_COMPILER "/opt/rocm/llvm/bin/clang++")
    endif()
    if (NOT DEFINED HIP_CLANG_PATH)
      set(HIP_CLANG_PATH "/opt/rocm/llvm/bin")
    endif()
  endif()
  include(CheckLanguage)
  check_language(HIP)
  find_package(hip)
  find_package(hipcub)
  find_package(rocprim)
  find_package(rocthrust)
  find_program(hip_HIPIFY_PERL_EXECUTABLE "hipify-perl")
  if(NOT hip_HIPIFY_PERL_EXECUTABLE)
    find_program(hip_HIPIFY_PERL_EXECUTABLE "hipify-perl" HINTS "/opt/rocm/bin")
  endif()
  if(ENABLE_HIP STREQUAL "AUTO")
    set_package_properties(hip PROPERTIES TYPE OPTIONAL)
    set_package_properties(hipcub PROPERTIES TYPE OPTIONAL)
    set_package_properties(rocprim PROPERTIES TYPE OPTIONAL)
    set_package_properties(rocthrust PROPERTIES TYPE OPTIONAL)
  else()
    set_package_properties(hip PROPERTIES TYPE REQUIRED)
    set_package_properties(hipcub PROPERTIES TYPE REQUIRED)
    set_package_properties(rocprim PROPERTIES TYPE REQUIRED)
    set_package_properties(rocthrust PROPERTIES TYPE REQUIRED)
  endif()
  if (CMAKE_HIP_COMPILER)
    enable_language(HIP)
    message(STATUS "HIP language enabled: ${CMAKE_HIP_COMPILER}")
  endif()
  if(hip_FOUND AND hipcub_FOUND AND rocthrust_FOUND AND rocprim_FOUND AND hip_HIPCC_EXECUTABLE AND hip_HIPIFY_PERL_EXECUTABLE)
    set(HIP_ENABLED ON)
    set_target_properties(roc::rocthrust PROPERTIES IMPORTED_GLOBAL TRUE)
    message(STATUS "HIP Found (${hip_HIPCC_EXECUTABLE} version ${hip_VERSION})")
    set(O2_HIP_CMAKE_CXX_FLAGS "-fgpu-defer-diag -mllvm -amdgpu-enable-lower-module-lds=false -Wno-invalid-command-line-argument -Wno-unused-command-line-argument -Wno-invalid-constexpr -Wno-ignored-optimization-argument -Wno-unused-private-field -Wno-pass-failed")
    set(O2_HIP_CMAKE_LINK_FLAGS "-Wno-pass-failed")
    string(REGEX REPLACE "(gfx1[0-9]+;?)" "" CMAKE_HIP_ARCHITECTURES "${CMAKE_HIP_ARCHITECTURES}") # ROCm currently doesn’t support integrated graphics
    if(HIP_AMDGPUTARGET)
      foreach(HIP_ARCH ${HIP_AMDGPUTARGET})
        set(O2_HIP_CMAKE_CXX_FLAGS "${O2_HIP_CMAKE_CXX_FLAGS} --offload-arch=${HIP_ARCH}")
        set(O2_HIP_CMAKE_LINK_FLAGS "${O2_HIP_CMAKE_LINK_FLAGS} --offload-arch=${HIP_ARCH}")
      endforeach()
      set(CMAKE_HIP_ARCHITECTURES "${HIP_AMDGPUTARGET}") # If GPU build is enforced we override autodetection
    endif()
    if(NOT DEFINED GPUCA_NO_FAST_MATH OR NOT ${GPUCA_NO_FAST_MATH})
      set(O2_HIP_CMAKE_CXX_FLAGS "${O2_HIP_CMAKE_CXX_FLAGS} -fgpu-flush-denormals-to-zero") # -ffast-math disabled, since apparently it leads to miscompilation and crashes in FollowLooper kernel
    endif()
    if (CMAKE_CXX_COMPILER MATCHES "bin/c\\+\\+\$" AND NOT CMAKE_CXX_COMPILER MATCHES "^/usr/bin")
      string(REGEX REPLACE "(.*)bin/c\\+\\+\$" "\\1" HIP_GCC_TOOLCHAIN_PATH "${CMAKE_CXX_COMPILER}")
      set(O2_HIP_CMAKE_CXX_FLAGS "${O2_HIP_CMAKE_CXX_FLAGS} --gcc-toolchain=${HIP_GCC_TOOLCHAIN_PATH}") # -ffast-math disabled, since apparently it leads to miscompilation and crashes in FollowLooper kernel
    endif()
    set(CMAKE_HIP_FLAGS "${O2_GPU_CMAKE_CXX_FLAGS_NOSTD} ${CMAKE_HIP_FLAGS} ${O2_HIP_CMAKE_CXX_FLAGS}")
    set(CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE} "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${CMAKE_HIP_FLAGS_${CMAKE_BUILD_TYPE}}")
  else()
    set(HIP_ENABLED OFF)
  endif()
  if(NOT HIP_ENABLED AND NOT ENABLE_HIP STREQUAL "AUTO")
    if (NOT hip_FOUND)
      message(WARNING "HIP not found")
    endif()
    if (NOT hipcub_FOUND)
      message(WARNING "hipcub not found")
    endif()
    if (NOT rocthrust_FOUND)
      message(WARNING "rocthrust not found")
    endif()
    if (NOT rocprim_FOUND)
      message(WARNING "rocprim not found")
    endif()
    if (NOT hip_HIPCC_EXECUTABLE)
      message(WARNING "hipcc executable not found")
    endif()
    if (NOT hip_HIPIFY_PERL_EXECUTABLE)
      message(WARNING "hipify-perl executable not found")
    endif()
    message(FATAL_ERROR "HIP requested with HIP_PATH=${HIP_PATH} but some of the above packages are not found")
  endif()

endif()

# if we end up here without a FATAL, it means we have found the "O2GPU" package
set(O2GPU_FOUND TRUE)
include("${CMAKE_CURRENT_LIST_DIR}/../GPU/GPUTracking/cmake/kernel_helpers.cmake")
