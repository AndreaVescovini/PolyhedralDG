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

  Watch();

  virtual ~Watch() = default;

  //! Outputs time from last start and stop
  friend std::ostream& operator <<(std::ostream & out, const Watch& w);

  //! Starts/reset  counting time
  void start();

  //! Stops counting time
  void stop();

  void reset();

  double getTime() const;

  double getTimeNow() const;

private:
  MyTimePoint startTime_;
  Nanosec measuredTime_;
};

} // namespace Timings

#endif // _WATCH_HPP_
