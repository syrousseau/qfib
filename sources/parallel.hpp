/// \file
// --------------------------------------------------------------------------
// This file is part of the reference implementation for the paper
//    QFib: Fast and Efficient Brain Tractogram Compression
//    C. Mercier*, S. Rousseau*, P. Gori, I. Bloch and T. Boubekeur
//    NeuroInformatics 2020
//    DOI: 10.1007/s12021-020-09452-0
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.
// --------------------------------------------------------------------------
#pragma once

#include <thread>
#include <mutex>
#include <future>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

///
/// \brief The Parallel class handles parallelism of loops, allowing for a functionality close to OpenMP (with some limitations)
///

class Parallel
{

private:

	///
	/// \brief createThread creates a new thread that will run the lambda l, with arguments args
	/// \param l is the lambda to run inside a new thread
	/// \param args are the arguments of the lambda
	///

	template<typename Lambda, typename... Args>
	void createThread(Lambda &&l, Args&&... args)
	{
		bool isReady = false;
		bool wasReceived = false;
		std::mutex mutex;
		std::condition_variable condition;

		unsigned currentThreadID = m_threads.size();
		m_threads.emplace_back([=, &mutex, &condition, &wasReceived, &isReady]()
		{
			std::unique_lock<std::mutex> lock(mutex);
			condition.wait(lock, [&]{return isReady;});
			m_IDs.insert(std::pair<std::thread::id, unsigned>(std::this_thread::get_id(), currentThreadID));

			wasReceived =true;
			lock.unlock();
			condition.notify_one();

			l(args...);

		});

		{
			std::lock_guard<std::mutex> lock(mutex);
			isReady = true;
		}
		condition.notify_one();
		{
			std::unique_lock<std::mutex> lock(mutex);
			condition.wait(lock, [&]{return wasReceived;});
		}

	}

	std::vector<std::thread> m_threads;
	std::map<std::thread::id, unsigned> m_IDs;
	unsigned m_nbThreads;

public:

	///
	/// \brief Parallel Constructor for a parallel context
	/// \param nbThreads Number of threads to use, 0 (default) means all available threads are used
	///

	Parallel(unsigned nbThreads = 0):m_nbThreads(nbThreads)
	{
		//If no number of threads specified, get the maximum available
		if (m_nbThreads == 0)
			m_nbThreads = std::max(unsigned(1), std::thread::hardware_concurrency());

		m_threads.reserve(m_nbThreads);
	}

	///
	/// \brief getID Returns current thread ID, as a value between 0 and nbThreads - 1 (useful for avoiding race condition)
	///

	unsigned getID()
	{
		return m_IDs[std::this_thread::get_id()];
	}

	///
	/// \brief For Parallelized for loop - Warning, works only in increasing order, with step positive, and end > start
	/// The type of start, end, and step has to be the same, but can be anything comparable with "<"
	/// \param start Starting value for the loop
	/// \param end Ending value for the loop (not included)
	/// \param step Step for increasing the loop from start to end
	/// \param func Lambda function to launch inside the loop (can call getID() to get current thread number)
	///

	template<typename Index, typename Callable>
	void For(Index start, Index end, Index step, Callable func)
	{
		if (start > end || step < 0)
		{
			std::cerr << "Parallel for only handle increasing loop order" << std::endl;
			return;
		}

		//[Helper] Inner loop
		auto launchRange = [&step, &func] (int k1, int k2) {
			for (Index k = k1; k < k2; k+=step) {
				func(k);
			}
		};

		Index n = (end - start) / step;
		Index slice = (Index) std::round(n / static_cast<double>(m_nbThreads));
		slice = std::max(slice, Index(1)) * step;
		Index i1 = start;
		Index i2 = std::min(start + slice, end);


		//Launch threads on func
		for (unsigned i=0; i<m_nbThreads - 1 && i1 < end; ++i)
		{
			createThread(launchRange, i1, i2);
			i1 = i2;
			i2 = std::min(i2 + slice, end);
		}
		if (i1 < end)
		{
			createThread(launchRange, i1, end);
		}

		//Wait for all threads to finish
		for (std::thread &t: m_threads)
		{
			if (t.joinable())
				t.join();
		}
		m_threads.clear();
		m_IDs.clear();
	}

};

#ifndef USE_OPENMP
	Parallel parallel;
#endif
