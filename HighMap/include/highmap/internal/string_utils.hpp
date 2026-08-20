/**
 * @file string_utils.hpp
 * @copyright Copyright (c) 2023 Otto Link. Distributed under the terms of the
 * GNU General Public License. The full license is in the file LICENSE,
 * distributed with this software.
 */
#pragma once
#include <filesystem>
#include <string>

namespace hmap
{

/**
 * @brief Adds a suffix to the filename of a given file path.
 *
 * This function appends a given suffix to the stem (base name without
 * extension) of a file while preserving the original directory and file
 * extension.
 *
 * @param  file_path The original file path.
 * @param  suffix    The suffix to append to the filename.
 * @return           A new std::filesystem::path with the modified filename.
 *
 * @note If the input file has no extension, the suffix is added directly to the
 * filename.
 *
 * **Example**
 * @code{.cpp}
 * std::filesystem::path path = "example.txt";
 * std::filesystem::path new_path = add_filename_suffix(path, "_backup");
 * std::cout << new_path; // Outputs "example_backup.txt"
 * @endcode
 */
std::filesystem::path add_filename_suffix(
    const std::filesystem::path &file_path,
    const std::string           &suffix);

/**
 * @brief Create a unique temporary directory.
 *
 * Creates and returns a unique directory inside the system temporary directory
 * using the given prefix.
 *
 * @param  prefix Directory name prefix.
 * @return        Path to the created temporary directory.
 *
 * @note The caller is responsible for deleting the directory.
 */
std::filesystem::path make_unique_temp_dir(const std::string &prefix);

/**
 * @brief Pads the beginning of a string with zeros until it reaches a specified
 * length.
 *
 * This function prepends '0' characters to the input string `str` so that the
 * resulting string has at least `n_zero` characters. If `str` is already equal
 * to or longer than `n_zero`, it is returned unchanged.
 *
 * @param  str    The input string to pad.
 * @param  n_zero The minimum total length of the resulting string after
 *                padding.
 * @return        A new string padded with leading zeros to reach the specified
 *                length.
 *
 * **Example**
 * @code{.cpp}
 * zfill("42", 5); // returns "00042"
 * zfill("12345", 5); // returns "12345"
 * @endcode
 */
std::string zfill(const std::string &str, int n_zero);

} // namespace hmap