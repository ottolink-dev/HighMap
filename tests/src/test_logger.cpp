#include <sstream>

#include "highmap/logger.hpp"

#include <gtest/gtest.h>

TEST(LoggerTest, BasicLogging)
{
  // Test that all logging functions compile and run cleanly
  hmap::log::trace("Trace test: {}", "ok");
  hmap::log::info("Info test: {}", 123);
  hmap::log::warn("Warn test: {}", 456);
  hmap::log::error("Error test: {}", 789);

  std::string dynamic_msg = "dynamic message";
  hmap::log::trace(dynamic_msg);
  hmap::log::info(dynamic_msg);
  hmap::log::warn(dynamic_msg);
  hmap::log::error(dynamic_msg);

  SUCCEED();
}

TEST(LoggerTest, FormattedOutput)
{
  // Capture std::cout to verify info/trace format
  std::stringstream buffer;
  std::streambuf   *old_cout = std::cout.rdbuf(buffer.rdbuf());

  hmap::log::info("Hello {}", "World");

  std::cout.rdbuf(old_cout);

#if HIGHMAP_ENABLE_LOGS
  std::string output = buffer.str();
  EXPECT_NE(output.find("[info ]"), std::string::npos);
  EXPECT_NE(output.find("Hello World"), std::string::npos);
  EXPECT_NE(output.find("test_logger.cpp"), std::string::npos);
  EXPECT_EQ(output.find("/test_logger.cpp"), std::string::npos);
#endif
}

TEST(LoggerTest, ErrorOutput)
{
  // Capture std::cerr to verify warn/error format
  std::stringstream buffer;
  std::streambuf   *old_cerr = std::cerr.rdbuf(buffer.rdbuf());

  hmap::log::error("Critical failure: code {}", 404);

  std::cerr.rdbuf(old_cerr);

#if HIGHMAP_ENABLE_LOGS
  std::string output = buffer.str();
  EXPECT_NE(output.find("[error]"), std::string::npos);
  EXPECT_NE(output.find("Critical failure: code 404"), std::string::npos);
  EXPECT_NE(output.find("test_logger.cpp"), std::string::npos);
#endif
}

TEST(LoggerTest, ExplicitSourceLocation)
{
  std::stringstream buffer;
  std::streambuf   *old_cout = std::cout.rdbuf(buffer.rdbuf());

  auto loc = std::source_location::current();
  hmap::log::info(loc, "Explicit loc formatted {}", 42);
  hmap::log::trace(loc, "Explicit loc plain message");

  std::string dyn = "dynamic string";
  hmap::log::info(loc, dyn);

  std::cout.rdbuf(old_cout);

#if HIGHMAP_ENABLE_LOGS
  std::string output = buffer.str();
  EXPECT_NE(output.find("[info ]"), std::string::npos);
  EXPECT_NE(output.find("Explicit loc formatted 42"), std::string::npos);
  EXPECT_NE(output.find("[trace]"), std::string::npos);
  EXPECT_NE(output.find("Explicit loc plain message"), std::string::npos);
  EXPECT_NE(output.find("dynamic string"), std::string::npos);
#endif
}

#if HIGHMAP_ENABLE_LOGS
TEST(LoggerTest, FixedWidthFormattingAndTruncation)
{
  EXPECT_EQ(hmap::log::detail::format_fixed_width("hello", 5), "hello");
  EXPECT_EQ(hmap::log::detail::format_fixed_width("hi", 5), "hi   ");
  EXPECT_EQ(hmap::log::detail::format_fixed_width("very_long_string", 8),
            "very_...");
  EXPECT_EQ(hmap::log::detail::format_fixed_width("abc", 2), "ab");
  EXPECT_EQ(hmap::log::detail::format_fixed_width("abc", 0), "abc");
}
#endif
