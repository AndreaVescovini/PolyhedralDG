/*!
    @file   Watch.hpp
    @author Andrea Vescovini
    @brief  Class that used to measure time
*/

#ifndef _WATCH_HPP_
#define _WATCH_HPP_

#include <chrono>
#include <iostream>

namespace Timings
{

/*!
    @brief Class used to measure time

    This class defines a "Watch" useful for measuring time in programs. You can
    call the method tart() to start the measure, Watch::stop() to stop
    it and the getTime(). After that you can reset() the Watch to 0
    or start() a new measure to be added to the preavious one.
*/

class Watch
{
public:
  //! Alias for the system clock
  using MyClock = std::chrono::high_resolution_clock;

  //! Alias fot the instant of time
  using MyTimePoint = std::chrono::time_point<MyClock>;

  //! Alias for nanoseconds
  using Nanosec = std::chrono::nanoseconds;

  /*!
      @brief Constructor

      The constructor initializes start time to now.
  */
  Watch();

  //! Copy constructor
  Watch(const Watch&) = default;

  //! Copy-assignment operator
  Watch& operator=(const Watch&) = default;

  //! Move constructor
  Watch(Watch&&) = default;

  //! Move-assignment operator
  Watch& operator=(Watch&&) = default;

  /*!
      @brief Overload for the ostream operator

      It prints the measured time between the last start and stop of the Watch.
  */
  friend std::ostream& operator <<(std::ostream & out, const Watch& w);

  //! Initialize start time to now
  void start();

  //! Stop counting time
  void stop();

  //! Reset to 0 the measured time
  void reset();

  /*!
      @brief  Get the measured time [ms]
      @return The mesured time expressed in microseconds.
  */
  double getTime() const;

  /*!
      @brief  Makes a measure now and get it [microseconds]
      @return The measured time between the last start and now, expressed in
              microseconds.
  */
  double getTimeNow() const;

  //! Destructor
  virtual ~Watch() = default;

private:
  //! Stores the start time
  MyTimePoint startTime_;

  //! Stores the measured time
  Nanosec measuredTime_;
};

} // namespace Timings

#endif // _WATCH_HPP_
