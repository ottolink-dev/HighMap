/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file dbg.hpp
 * @author  Otto Link (otto.link.bv@gmail.com)
 * @brief Header file containing the implementation of the Recorder and Timer
 * classes for high-resolution timing.
 * @version 0.1
 * @date 2023-04-30
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <string>

namespace hmap
{

/**
 * @brief The Recorder class is responsible for recording timing information for
 * individual events.
 */
struct Recorder
{
public:
  /**
   * @brief Constructs a new Recorder object.
   *
   * @param name The name of the event to be recorded.
   */
  Recorder(std::string name);

  /**
   * @brief Outputs the timing data to the console.
   */
  void dump();

  /**
   * @brief Starts the timer for this Recorder instance.
   */
  void start();

  /**
   * @brief Stops the timer for this Recorder instance and updates the
   * total elapsed time.
   */
  void stop();

  std::string name; ///< The name of the event.
  int nb_calls = 0; ///< The number of times the event has been recorded.
  std::chrono::high_resolution_clock::time_point
        t0;          ///< The start time of the event.
  float total = 0.f; ///< The total time recorded for the event.
};

/**
 * @brief The Timer class is a singleton that manages multiple Recorders and
 * provides an interface for timing events.
 *
 * The Timer class allows you to measure the duration of multiple events by
 * using start and stop commands. Each event is identified by a unique name. The
 * Timer class maintains a collection of these events and their corresponding
 * durations. The class is designed as a singleton, meaning only one instance of
 * Timer will exist throughout the lifetime of the program.
 *
 * ### Usage Example:
 *
 * @code
 * #include "timer.hpp"
 *
 * int main() {
 *     // Start timing an event named "step 1"
 *     Timer::Start("step 1");
 *
 *     // Perform some work here...
 *
 *     // Stop timing the event named "step 1"
 *     Timer::Stop("step 1");
 *
 *     // Start and stop another event Timer::Start("step 2");
 *
 *     // Perform more work here...
 *
 *     Timer::Stop("step 2");
 *
 *     // Dump the timing results for all recorded events Timer::Dump();
 *
 *     return 0;
 * }
 * @endcode
 */
class Timer
{
public:
  /**
   * @brief Gets the singleton instance of the Timer class.
   *
   * @return Timer& Reference to the singleton instance.
   */
  static Timer &get_instance();

  /**
   * @brief Starts a timer for the specified event name.
   *
   * @param name The name of the event to start timing.
   */
  static void Start(const std::string &name);

  /**
   * @brief Stops the timer for the specified event name.
   *
   * @param name The name of the event to stop timing.
   */
  static void Stop(const std::string &name);

  static void Clear();

  /**
   * @brief Dumps the timing information for all recorded events to the console.
   */
  static void Dump();

  static std::map<std::string, float> DumpDurations();

  std::map<std::string, Recorder *> get_records() const
  {
    return this->records;
  }

private:
  /**
   * @brief Constructs a new Timer object. This constructor is private to
   * enforce the singleton pattern.
   *
   * @param sid An optional identifier for the Timer instance.
   */
  Timer(std::string sid = "");

  /**
   * @brief Destroys the Timer object and dumps all timing data.
   */
  ~Timer();

  /**
   * @brief Starts a timer for the specified event name.
   *
   * @param name The name of the event to start timing.
   */
  void start(const std::string &name);

  /**
   * @brief Stops the timer for the specified event name.
   *
   * @param name The name of the event to stop timing.
   */
  void stop(const std::string &name);

  /**
   * @brief Dumps the timing information for all recorded events to the console.
   */
  void dump();

  std::map<std::string, float> dump_durations() const;

  // Deleting the copy constructor and assignment operator to enforce singleton
  // pattern.
  Timer(const Timer &) = delete;
  Timer &operator=(const Timer &) = delete;

private:
  std::string sid; ///< An optional identifier for the Timer instance.
  std::map<std::string, Recorder *>
      records; ///< A map of event names to their corresponding Recorder
               ///< objects.
  std::list<Recorder>
      data; ///< A list of Recorder objects that store timing information.
  int current_level = 0; ///< Current nesting level (if applicable).
};

struct ScopedTimer
{
  ScopedTimer(const char *id) : id(id)
  {
    Timer::Start(id);
  }
  ~ScopedTimer()
  {
    Timer::Stop(id);
  }
  const char *id;
};

} // namespace hmap
