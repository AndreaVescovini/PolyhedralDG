#include "Watch.hpp"

// #include<iostream>

namespace Timings
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
  if(finalTime < 1e3)
    out << "Elapsed Time = " << finalTime << " microsec"<< std::endl;
  else if(finalTime < 1e6)
    out << "Elapsed Time = " << finalTime * 1e-3 << " millisec"<< std::endl;
  else if(finalTime < 1e9)
    out << "Elapsed Time = " << finalTime * 1e-6 << " sec"<< std::endl;

  return out;
}

} // namespace Timings
