/*!
    @file   Watch.cpp
    @author Andrea Vescovini
    @brief  Implementation for the class Watch
*/

#include "Watch.hpp"

namespace Utilities
{

Watch::Watch()
  : startTime_{MyClock::now()}, measuredTime_{0} {}

void Watch::start()
{
  startTime_ = MyClock::now();
}

void Watch::stop()
{
  measuredTime_ += std::chrono::duration_cast<Nanosec>(MyClock::now() - startTime_);
}

void Watch::reset()
{
  measuredTime_ = MyClock::duration::zero();
}

double Watch::getTime() const
{
  return measuredTime_.count() / 1000.0;
}

double Watch::getTimeNow() const
{
  const Nanosec timeSpan = measuredTime_ + std::chrono::duration_cast<Nanosec>(MyClock::now() - startTime_);
  return timeSpan.count() / 1000.0;
}

std::ostream& operator<<(std::ostream& out, const Watch& w)
{
  double finalTime = w.getTime();
  out << "Elapsed Time = ";

  if(finalTime < 1e3)
    out << finalTime << " microsec";
  else if(finalTime < 1e6)
    out << finalTime * 1e-3 << " millisec";
  else if(finalTime < 1e9)
    out << finalTime * 1e-6 << " sec";

  return out;
}

} // namespace Utilities
