#include "lupnt/core/file.h"

namespace lupnt {
  std::filesystem::path GetDataPath() {
    const char* data_path_env = std::getenv("LUPNT_DATA_PATH");
    if (data_path_env == nullptr) {
      std::string msg = "Environment variable LUPNT_DATA_PATH is not set.";
      throw std::runtime_error(msg);
    }
    return std::filesystem::path(data_path_env);
  }

  std::filesystem::path GetOutputPath(std::string output_dir) {
    const char* output_path_env = std::getenv("LUPNT_OUTPUT_PATH");
    std::filesystem::path output_path;
    if (output_path_env == nullptr) {
      output_path = GetDataPath() / "output";
    } else {
      output_path = std::filesystem::path(output_path_env);
    }
    if (!output_dir.empty()) {
      output_path /= output_dir;
    }

    if (!std::filesystem::exists(output_path)) {
      std::filesystem::create_directories(output_path);
    }
    return output_path;
  }

  H5Easy::File GetCacheFile(std::filesystem::path cache_path, bool recompute) {
    if (recompute && std::filesystem::exists(cache_path)) std::filesystem::remove(cache_path);
    auto open_mode = (recompute || !std::filesystem::exists(cache_path))
                         ? H5Easy::File::OpenOrCreate
                         : H5Easy::File::ReadOnly;
    H5Easy::File cache_file(cache_path, open_mode);
    return cache_file;
  }

  std::optional<std::filesystem::path> FindFileInDir(const std::filesystem::path& base_path,
                                                     std::string_view filename) {
    for (const auto& entry : std::filesystem::recursive_directory_iterator(base_path)) {
      if (entry.is_directory()) continue;
      if (entry.path().filename().string() == filename) return entry.path();
      if (entry.path().stem().string() == filename) return entry.path();
    }
    return std::nullopt;
  };

  std::filesystem::path GetFilePath(std::string_view filename) {
    auto filepath = FindFileInDir(GetDataPath(), filename);
    if (!filepath.has_value()) {
      std::string msg = "File not found: " + std::string(filename);
      throw std::runtime_error(msg);
    }
    return filepath.value();
  }

  std::filesystem::path GetCspiceKernelDir() { return GetDataPath() / "ephemeris"; }
  std::filesystem::path GetAsciiKernelDir() { return GetDataPath() / "ephemeris" / "ascii"; }
  std::filesystem::path GetGroundStationDataDir() { return GetDataPath() / "ground_station"; }
  std::filesystem::path GetConfigFileDir() { return GetDataPath() / ".." / "config"; }

  std::chrono::time_point<std::chrono::high_resolution_clock> GetSystemTime() {
    return std::chrono::high_resolution_clock::now();
  }
  // print duration between start_time and end_time
  std::string PrintDuration(const std::chrono::duration<double>& duration) {
    auto total_seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    auto hours = total_seconds / 3600;
    auto minutes = (total_seconds % 3600) / 60;
    auto seconds = total_seconds % 60;

    std::string str;
    if (hours > 0) str += std::to_string(hours) + "h ";
    if (minutes > 0) str += std::to_string(minutes) + "m ";
    str += std::to_string(seconds) + "s ";
    return str;
  }

  template <typename T> T OpenFile(const std::filesystem::path& filepath) {
    T file(filepath);
    if (!file.is_open()) throw std::runtime_error("Unable to open file: " + filepath.string());
    return file;
  }
  template std::ifstream OpenFile<std::ifstream>(const std::filesystem::path& filepath);
  template std::ofstream OpenFile<std::ofstream>(const std::filesystem::path& filepath);

  size_t CountLines(const std::filesystem::path& filepath) {
    std::ifstream file = OpenFile<std::ifstream>(filepath);
    size_t line_count = 0;
    std::string line;
    while (std::getline(file, line)) ++line_count;
    file.close();
    return line_count;
  }

}  // namespace lupnt
