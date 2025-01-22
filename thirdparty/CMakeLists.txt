# ---- Add dependencies via CPM ----

# PackageProject
CPMAddPackage(
  NAME PackageProject.cmake
  GITHUB_REPOSITORY TheLartians/PackageProject.cmake
  VERSION 1.11.2
)

# Eigen
CPMAddPackage(
  NAME Eigen
  GIT_TAG 3.4
  GITLAB_REPOSITORY libeigen/eigen
  OPTIONS "EIGEN_BUILD_DOC OFF" "EIGEN_BUILD_BTL OFF" "EIGEN_BUILD_PKGCONFIG OFF"
)

# autodiff
CPMAddPackage(
  NAME autodiff
  GIT_TAG v1.1.2
  GITHUB_REPOSITORY autodiff/autodiff
  OPTIONS "AUTODIFF_BUILD_TESTS OFF" "AUTODIFF_BUILD_EXAMPLES OFF" "AUTODIFF_BUILD_DOCS OFF"
          "AUTODIFF_BUILD_PYTHON OFF"
)

# cspice
CPMAddPackage(
  NAME cspice
  GITHUB_REPOSITORY gcasadesus/cspice
  GIT_TAG master
  OPTIONS "CMAKE_BUILD_TYPE Release"
)

# matplotplusplus
CPMAddPackage(
  NAME matplotplusplus
  GIT_TAG v1.2.1
  GITHUB_REPOSITORY alandefreitas/matplotplusplus
  OPTIONS "MATPLOTPP_BUILD_INSTALLER ON" "MATPLOTPP_BUILD_PACKAGE ON" "CMAKE_BUILD_TYPE Release"
)
install(
  TARGETS cimg
  EXPORT Matplot++Targets
  INCLUDES
  DESTINATION include
)
install(
  TARGETS nodesoup
  EXPORT Matplot++Targets
  INCLUDES
  DESTINATION include
)

# HighFive
CPMAddPackage(
  NAME HighFive
  GIT_TAG v2.10.0
  GITHUB_REPOSITORY BlueBrain/HighFive
  OPTIONS "HIGHFIVE_UNIT_TESTS OFF" "HIGHFIVE_EXAMPLES OFF" "HIGHFIVE_BUILD_DOCS OFF"
)

# EigenRand
CPMAddPackage(
  NAME EigenRand
  GIT_TAG v0.5.1
  GITHUB_REPOSITORY bab2min/EigenRand
  OPTIONS "EIGENRAND_BUILD_TEST OFF" "EIGENRAND_BUILD_BENCHMARK OFF"
)

# magic_enum
CPMAddPackage(
  NAME magic_enum
  GIT_TAG v0.9.7
  GITHUB_REPOSITORY Neargye/magic_enum
  OPTIONS "MAGIC_ENUM_OPT_BUILD_EXAMPLES OFF" "MAGIC_ENUM_OPT_BUILD_TESTS OFF"
          "MAGIC_ENUM_OPT_INSTALL ON"
)
install(
  TARGETS magic_enum
  EXPORT magic_enumTargets
  INCLUDES
  DESTINATION include
)

# Yaml-cpp
CPMAddPackage(
  NAME yaml-cpp
  GITHUB_REPOSITORY jbeder/yaml-cpp
  GIT_TAG 0.8.0 # or the version you want
  OPTIONS "YAML_CPP_INSTALL ON" "YAML_CPP_BUILD_TESTS OFF" "YAML_CPP_BUILD_CONTRIB OFF"
          "YAML_CPP_BUILD_TOOLS OFF" "YAML_CPP_FORMAT_SOURCE OFF"
)
install(
  TARGETS yaml-cpp
  EXPORT yaml-cppTargets
  INCLUDES
  DESTINATION include
)

set_target_properties(cspice PROPERTIES FOLDER thirdparty)
set_target_properties(yaml-cpp PROPERTIES FOLDER thirdparty)
set_target_properties(matplot PROPERTIES FOLDER thirdparty)
set_target_properties(nodesoup PROPERTIES FOLDER thirdparty)
set_target_properties(demos PROPERTIES FOLDER thirdparty)
set_target_properties(uninstall PROPERTIES FOLDER thirdparty)

# LuPNT dependencies used for installing
set(LUPNT_DEPENDENCIES
    "Eigen3 3.4;autodiff 1.1;cspice;Matplot++ 1.2;OpenMP;HighFive 2.10;Boost;yaml-cpp 0.8.0;magic_enum 0.9.7;EigenRand 0.5.1;OpenCV"
    PARENT_SCOPE
)
