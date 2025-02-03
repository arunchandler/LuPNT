/**
 * @file file.h
 * @author Stanford NAV LAB
 * @brief File access utils
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  std::filesystem::path GetDataPath();
  std::filesystem::path GetOutputPath(std::string output_dir);
  H5Easy::File GetCacheFile(std::filesystem::path cache_path, bool recompute);

  template <int N, int M, typename scalar, typename Func, typename... Args>
  Matrix<scalar, N, M> LoadOrRecompute(const std::string& name, H5Easy::File& cache_file,
                                       bool recompute, Func func, Args&&... args) {
// #pragma omp critical
//     if (!recompute && cache_file.exist(name)) {
//       return H5Easy::load<Matrix<double, N, M>>(cache_file, name);
//     }
//     Matrix<scalar, N, M> result = func(std::forward<Args>(args)...);
// #pragma omp critical
//     H5Easy::dump(cache_file, name, result.template cast<double>(), H5Easy::DumpMode::Overwrite);
//     return result;
//   }

    Matrix<scalar, N, M> result;
    bool found = false;
    Matrix<double, N, M> cachedResult;

// 1) Use a critical block to safely check cache
    #pragma omp critical
    {
        if (!recompute && cache_file.exist(name)) {
            cachedResult = H5Easy::load<Matrix<double, N, M>>(cache_file, name);
            found = true;
        }
    }

// 2) If we found the matrix in cache, return it
    if (found) {
        return cachedResult.template cast<scalar>();
    }

    // 3) Otherwise, compute new result outside of the critical region
    result = func(std::forward<Args>(args)...);

    // 4) Another critical block to safely save to cache
    #pragma omp critical
    {
        H5Easy::dump(
            cache_file,
            name,
            result.template cast<double>(),
            H5Easy::DumpMode::Overwrite
        );
    }

// 5) Finally, return the computed result
  return result;
}

std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& base_path,
                                                    std::string_view filename);
std::filesystem::path GetFilePath(std::string_view filename);
std::filesystem::path GetCspiceKernelDir();
std::filesystem::path GetAsciiKernelDir();
std::filesystem::path GetGroundStationDataDir();
std::filesystem::path GetConfigFileDir();

std::chrono::time_point<std::chrono::high_resolution_clock> GetSystemTime();
std::string PrintDuration(const std::chrono::duration<double>& duration);

template <typename T> T OpenFile(const std::filesystem::path& filepath);
size_t CountLines(const std::filesystem::path& filepath);

}  // namespace lupnt
