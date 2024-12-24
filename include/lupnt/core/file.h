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
#include <string>
#include <vector>

#include "lupnt/core/definitions.h"
#include "lupnt/numerics/math_utils.h"

namespace lupnt {

  std::filesystem::path GetDataPath();
  std::filesystem::path GetOutputPath(std::string output_dir);
  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& base_path,
                                                     std::string_view filename);
  std::filesystem::path GetFilePath(std::string_view filename);
  std::filesystem::path GetCspiceKernelDir();
  std::filesystem::path GetAsciiKernelDir();
  std::filesystem::path GetGroundStationDataDir();

  std::chrono::time_point<std::chrono::high_resolution_clock> GetSystemTime();
  std::string PrintDuration(const std::chrono::duration<double>& duration);

  template <typename T> T OpenFile(const std::filesystem::path& filepath);
  size_t CountLines(const std::filesystem::path& filepath);

}  // namespace lupnt
