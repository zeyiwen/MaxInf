#include <stdio.h>
#include <sys/timeb.h>

static struct _timeb start_time, finish_time;

// start the timer
void start_timer(void) {
	_ftime_s(&start_time);
}

// stop the timer and print the time elasped since last start() call

double stop_timer(int i) {
	_ftime_s(&finish_time);
	double duration = finish_time.time * 1000 + finish_time.millitm - (start_time.time  * 1000 + start_time.millitm);
	return duration;
}

