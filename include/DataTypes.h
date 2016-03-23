#pragma once
#ifndef DATA_TYPES_H
#define DATA_TYPES_H


// libraries
#include <cstdlib>
#include <vector>
#include <map>
#include <iostream>                  // for std::cout
#include <iomanip>
#include <utility>                   // for std::pair
#include <set>

// numerics
#include <limits>

// namespaces
// using namespace std;


// typedefs
enum status_t { FILE_MISSING, FILE_OK, INFEAS,INFEAS_OR_UNBND,FEAS, UNDEF_ERROR, FRACTIONAL, INTEGRAL, NO_ROUND }; // trace status_t of methods
enum color { WHITE, GRAY, BLACK }; // usually used in graph algorithms
typedef unsigned int _index;
typedef unsigned int integer;
typedef float  sreal;
typedef double dreal;
// typedef  std::map< _index,std::set<_index> > dict;
typedef std::string filename;


// string parsing


// use graphs
#include "GraphX.h"

#define USE_HOLGER_TIMER

#ifdef USE_BOOST_TIMER
#include <boost/timer/timer.hpp>


class Timer
{
private:
	boost::timer::cpu_timer t;
	boost::timer::cpu_times time_measure;

public:
	// default constructor
	Timer(){t.start(); t.stop();};

	// reset timer
	void reset() { t.stop(); }
	
	// start timer
	void start() { 	t.start(); }
	

	// continue timer
	void cont(){t.resume();}

	// stop timer
	void stop() { t.stop(); time_measure = t.elapsed(); }

	// reading data
	float secs() const { return (float)(time_measure.wall / 1000000000.0); } /* in seconds */




};

#endif

#ifdef USE_HOLGER_TIMER
// timer
#if defined( _MSC_VER )
#include<Windows.h>
#include <stdint.h>

struct timezone {
	int tz_minuteswest;     /* minutes west of Greenwich */
	int tz_dsttime;         /* type of DST correction */
};

int gettimeofday(struct timeval * tp, struct timezone * tzp);


#else 
#include <sys/time.h>

#endif



//! A SIMPLE CLASS FOR TIME MEASUREMENTS.
//
class Timer
{
private:

	//! The timer value (initially zero)
	off_t _usecs;

	//! The timer value at the last mark set (initially zero)
	off_t _usecs_at_mark;

	//! Used by the gettimeofday command.
	struct timeval _tstart;

	//! Used by the gettimeofday command.
	struct timeval _tend;

	//! Used by the gettimeofday command.
	struct timezone _tz;

	//! Indicates whether a measurement is running.
	bool _running;

public:

	//! The default constructor.
	Timer() { reset(); }

	//! Resets the timer value to zero and stops the measurement.
	void reset() { _usecs = _usecs_at_mark = 0; _running = false; }

	//! Mark the current point in time (to be considered by next usecs_since_mark)
	void mark() { stop(); _usecs_at_mark = _usecs; cont(); }

	//! Resets the timer value to zero and starts the measurement.
	void start()
	{
		_usecs = _usecs_at_mark = 0;
		gettimeofday(&_tstart, &_tz);
		_running = true;
	}

	//! Continues the measurement without resetting the timer value (no effect it running)
	void cont()
	{
		if (_running == false)
		{
			gettimeofday(&_tstart, &_tz);
			_running = true;
		}
	}

	//! Stops the measurement (does *not* return the timer value anymore)
	void stop()
	{
		gettimeofday(&_tend, &_tz);
		_running = false;
		_usecs += (off_t)(1000000) * (off_t)(_tend.tv_sec - _tstart.tv_sec) + (off_t)(_tend.tv_usec - _tstart.tv_usec);
		//return _usecs;
	}

	//! Time at last stop (initially zero)
	off_t value() const { return _usecs; } /* in microseconds */
	off_t usecs() const { return _usecs; } /* in microseconds */
	off_t msecs() const { return _usecs / 1000; } /* in milliseconds */
	float msecs_float() const { return (float)(_usecs / 1000.0); }
	float secs() const { return (float)(_usecs / 1000000.0); } /* in seconds */

	//! Time from last mark to last stop (initially zero)
	off_t usecs_since_mark() const { return _usecs - _usecs_at_mark; }

};

#endif

#endif

#if defined(USE_HOLGER_TIMER) || defined(USE_BOOST_TIMER)
#define USE_TIMER
#endif
