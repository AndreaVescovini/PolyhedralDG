#ifndef _WATCH_HPP_
#define _WATCH_HPP_

#include <chrono>
#include <iostream>

namespace Timings
{

class Watch
{
public:
  using MyClock = std::chrono::high_resolution_clock;
  using MyTimePoint = std::chrono::time_point<MyClock>;
  using Nanosec = std::chrono::nanoseconds;

  // Constructor. Initializes start time to now.
  Watch();

  // Default copy-constructor.
  Watch(const Watch&) = default;

  // Default copy-assignment operator.
  Watch& operator=(const Watch&) = default;

  // Default move constructor.
  Watch(Watch&&) = default;

  // Default move-assignment operator.
  Watch& operator=(Watch&&) = default;

  //! Outputs time from last start and stop
  friend std::ostream& operator <<(std::ostream & out, const Watch& w);

  //! Initializes start time to now.
  void start();

  //! Stops counting time
  void stop();

  // Reset to zero the measured time.
  void reset();

  // Get the measured time [microseconds]
  double getTime() const;

  // Makes a measure now and get [microseconds]
  double getTimeNow() const;

  // Default virtual destructor.
  virtual ~Watch() = default;

private:
  MyTimePoint startTime_;
  Nanosec measuredTime_;
};

} // namespace Timings

#endif // _WATCH_HPP_
