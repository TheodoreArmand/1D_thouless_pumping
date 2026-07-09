cmake_minimum_required(VERSION 3.16)
project(ecg1d_thouless_pumping LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core PATHS /usr/include/eigen3 /usr/local/include/eigen3)
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR "Eigen3 not found. Install Eigen or provide Eigen3_DIR/EIGEN3_INCLUDE_DIR.")
    endif()
endif()

find_package(OpenMP REQUIRED)

set(ECG_SRC
    src/basis_params.cpp
    src/csv_utils.cpp
    src/derivatives.cpp
    src/hamiltonian.cpp
    src/hamiltonian_gradient.cpp
    src/interaction_kernels.cpp
    src/observable_derivatives.cpp
    src/pair_cache.cpp
    src/permutation.cpp
    src/progress_output_n1.cpp
    src/progress_output_n2.cpp
    src/realtime_tdvp.cpp
    src/run_report.cpp
    src/snapshot_io.cpp
    src/tdvp_solver.cpp
    src/trace_io.cpp
)

set(PUMPCONFIG_SRC
    pumpconfig/pump_config.cpp
    pumpconfig/legacy_prb_3_3.cpp
    pumpconfig/lohse_10_5.cpp
    pumpconfig/lohse_n2_free.cpp
    pumpconfig/lohse_n2_gauss.cpp
)

set(PUMP2_SRC
    nointeraction_src/pump_common.cpp
)

add_executable(ecg1d_thouless_pumping
    main.cpp
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_thouless_pumping PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_thouless_pumping PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_thouless_pumping PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_thouless_pumping PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_thouless_pumping PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_thouless_pumping PRIVATE -O2 -march=x86-64-v3)

add_executable(ecg1d_legacy_prb_3_3_300step
    txt_task/ecg1d_legacy_prb_3_3_300step.cpp
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_legacy_prb_3_3_300step PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_legacy_prb_3_3_300step PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_legacy_prb_3_3_300step PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_legacy_prb_3_3_300step PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_legacy_prb_3_3_300step PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_legacy_prb_3_3_300step PRIVATE -O2 -march=x86-64-v3)

add_executable(ecg1d_2ninteraction
    main_2ninteraction.cpp
    ${PUMP2_SRC}
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_2ninteraction PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_2ninteraction PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_2ninteraction PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_2ninteraction PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_2ninteraction PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_2ninteraction PRIVATE -O2 -march=x86-64-v3)

add_executable(ecg1d_2gaussian
    main_2gaussian.cpp
    2gaussian_src/gaussian_terms.cpp
    ${PUMP2_SRC}
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_2gaussian PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_2gaussian PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_2gaussian PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_2gaussian PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_2gaussian PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_2gaussian PRIVATE -O2 -march=x86-64-v3)

add_executable(ecg1d_vs3_n2_dt001_full
    txt_task/ecg1d_vs3_n2_dt001_full.cpp
    2gaussian_src/gaussian_terms.cpp
    ${PUMP2_SRC}
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_vs3_n2_dt001_full PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_vs3_n2_dt001_full PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_vs3_n2_dt001_full PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_vs3_n2_dt001_full PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_vs3_n2_dt001_full PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_vs3_n2_dt001_full PRIVATE -O2 -march=x86-64-v3)

add_executable(ecg1d_vs3_n2_dt001_fine5_full
    txt_task/ecg1d_vs3_n2_dt001_fine5_full.cpp
    2gaussian_src/gaussian_terms.cpp
    ${PUMP2_SRC}
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_vs3_n2_dt001_fine5_full PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_vs3_n2_dt001_fine5_full PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_vs3_n2_dt001_fine5_full PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_vs3_n2_dt001_fine5_full PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_vs3_n2_dt001_fine5_full PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_vs3_n2_dt001_fine5_full PRIVATE -O2 -march=x86-64-v3)

add_executable(ecg1d_vs3_n2_dt001_fine20_full
    txt_task/ecg1d_vs3_n2_dt001_fine20_full.cpp
    2gaussian_src/gaussian_terms.cpp
    ${PUMP2_SRC}
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_vs3_n2_dt001_fine20_full PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_vs3_n2_dt001_fine20_full PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_vs3_n2_dt001_fine20_full PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_vs3_n2_dt001_fine20_full PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_vs3_n2_dt001_fine20_full PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_vs3_n2_dt001_fine20_full PRIVATE -O2 -march=x86-64-v3)

add_executable(ecg1d_gauss_gradient_check
    txt_task/ecg1d_gauss_gradient_check.cpp
    2gaussian_src/gaussian_terms.cpp
    ${PUMP2_SRC}
    ${PUMPCONFIG_SRC}
    ${ECG_SRC}
)

target_include_directories(ecg1d_gauss_gradient_check PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src
)

if(Eigen3_FOUND)
    target_link_libraries(ecg1d_gauss_gradient_check PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
else()
    target_include_directories(ecg1d_gauss_gradient_check PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(ecg1d_gauss_gradient_check PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(ecg1d_gauss_gradient_check PRIVATE EIGEN_DONT_PARALLELIZE)
target_compile_options(ecg1d_gauss_gradient_check PRIVATE -O2 -march=x86-64-v3)
