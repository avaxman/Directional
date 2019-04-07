#ifndef DIRECTIONAL_SIMPLETIMER_H
#define DIRECTIONAL_SIMPLETIMER_H
#include <chrono>
class SimpleTimer
{
public:
	using Clock = std::chrono::high_resolution_clock;
	using Time = std::chrono::high_resolution_clock::time_point;
private:
	Time t;
	double m_total = 0;
	bool m_running = false;
public:
	void start()
	{
		if (m_running) return;
		m_running = true;
		t = Clock::now();
	}
	SimpleTimer& reset()
	{
		m_running = false;
		m_total = 0;
		return *this;
	}
	void stop()
	{
		if (!m_running) return;
		m_running = false;
		std::chrono::duration<double> diff = (Clock::now() - t);
		m_total += std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
	}
	double elapsed() const
	{
		return m_total;
	}
};



#endif