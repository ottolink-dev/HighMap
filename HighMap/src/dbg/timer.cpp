/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <chrono>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <utility>

#include "highmap/dbg/timer.hpp"

namespace hmap
{

Recorder::Recorder(std::string name) : name(name)
{
}

void Recorder::dump()
{
  std::cout << std::setw(20) << this->name;
  std::cout << std::setw(10) << this->nb_calls;
  std::cout << std::setw(20) << this->total / (float)this->nb_calls << " ms";
  std::cout << std::endl;
}

void Recorder::start()
{
  this->t0 = std::chrono::high_resolution_clock::now();
  this->nb_calls++;
}

void Recorder::stop()
{
  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();
  this->total += (float)std::chrono::duration_cast<std::chrono::nanoseconds>(
                     t1 - t0)
                     .count() *
                 1e-6f;
}

// Static method to get the singleton instance
Timer &Timer::get_instance()
{
  static Timer instance; // This is the singleton instance
  return instance;
}

// Static methods to start and stop timers
void Timer::Start(const std::string &name)
{
  get_instance().start(name);
}

void Timer::Stop(const std::string &name)
{
  get_instance().stop(name);
}

void Timer::Clear()
{
  get_instance().records.clear();
}

void Timer::Dump()
{
  get_instance().dump();
}

std::map<std::string, float> Timer::DumpDurations()
{
  return get_instance().dump_durations();
}

// Private constructor to prevent instantiation
Timer::Timer(std::string sid) : sid(sid)
{
}
Timer::~Timer()
{
  this->dump();
}

// The actual methods for starting and stopping
void Timer::start(const std::string &name)
{
  if (records.find(name) == records.end())
  {
    Recorder new_recorder(name);
    data.push_back(new_recorder);
    records[name] = &data.back();
  }
  records[name]->start();
}

void Timer::stop(const std::string &name)
{
  if (records.find(name) != records.end())
    records[name]->stop();
  else
    std::cout << "Warning! Trying to stop an unknown timer: " << name
              << std::endl;
}

void Timer::dump()
{
  if (records.empty()) return;

  std::cout << "Timer dump: " << this->sid << std::endl;
  for (auto &n : records)
    n.second->dump();
}

std::map<std::string, float> Timer::dump_durations() const
{
  std::map<std::string, float> out;

  for (auto &n : records)
    out[n.first] = n.second->total / (float)n.second->nb_calls;

  return out;
}

} // namespace hmap
