#pragma once

#include <boost/thread/thread.hpp>
#include <boost/thread/condition.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "ARRAY.h"

class MULTITHREADING
{
public: // Type Define
	typedef boost::thread	THREAD, *pTHREAD;

public: // Essential Data 
	int						num_threads;

	pTHREAD*				thread_list;

	T*						sync_value;
	int*					sync_value_int;
	int*					sync_value_int_array;

	int*					lock_int;

	int						sync_index;
	int						sync_index2;

public: // For Synchronization
	T						sync_value_temp;

	boost::mutex			sync_mutex;
	boost::condition		cond;
	int						num_of_waiting_threads;

public: // For domain decomposition
	ARRAY<int>				start_ix_1D, end_ix_1D;

public: // For other operations
	boost::random::mt19937* random_number_generator_per_thread;

public: // Constructors and Destructor
	MULTITHREADING(void)
		: num_threads(0), sync_value(0), sync_value_int(0), sync_value_int_array(0), thread_list(0), num_of_waiting_threads(0), random_number_generator_per_thread(0), lock_int(0)
	{}

	~MULTITHREADING(void)
	{
		DELETE_ARRAY(sync_value);
		DELETE_ARRAY(sync_value_int);
		DELETE_ARRAY(sync_value_int_array);

		if(thread_list != 0)
		{
			for(int i = 0; i < num_threads; i++)
			{
				if(thread_list[i] != 0)
				{
					delete thread_list[i];
					thread_list[i] = 0;
				}
			}
			delete [] thread_list;
			thread_list = 0;
		}

		DELETE_ARRAY(random_number_generator_per_thread);
		DELETE_ARRAY(lock_int);
	}

public: // Initialization Functions
	void Initialize(const int& num_threads_input)
	{
		num_threads = num_threads_input;

		num_of_waiting_threads = 0;

		if(sync_value != 0) delete [] sync_value;
		sync_value = new T [num_threads];

		if(sync_value_int != 0) delete [] sync_value_int;
		sync_value_int = new int [num_threads];

		if(sync_value_int_array != 0) delete [] sync_value_int_array;
		sync_value_int_array = new int [num_threads*num_threads];

		if(thread_list != 0) delete [] thread_list;
		thread_list = new pTHREAD [num_threads];
		for(int i = 0; i < num_threads; i++) thread_list[i] = 0;

		random_number_generator_per_thread = new boost::random::mt19937 [num_threads];

		DELETE_ARRAY(lock_int);
		lock_int = new int [num_threads];
		for (int i = 0; i < num_threads; i++)
		{
			lock_int[i] = -1;
		}

		// For 1D domain decomposition (for example, SPH particle list)
		start_ix_1D.Initialize(num_threads);
		end_ix_1D.Initialize(num_threads);
	}

public: // Member Functions
	void Sync(const int& thread_id)
	{
		boost::mutex::scoped_lock lock(sync_mutex);

		num_of_waiting_threads++;

		if(num_of_waiting_threads == num_threads)
		{
			num_of_waiting_threads = 0;
			cond.notify_all();
		}
		else
			cond.wait(lock);
	}

	void Sync(const int& thread_id, const int& max_thread)
	{
		boost::mutex::scoped_lock lock(sync_mutex);

		num_of_waiting_threads++;

		if(num_of_waiting_threads == max_thread)
		{
			num_of_waiting_threads = 0;
			cond.notify_all();
		}
		else
			cond.wait(lock);
	}

	void InitializeSyncIndex(const int& thread_id)
	{
		if(thread_id == 0) sync_index = 0;

		Sync(thread_id);
	}

	void InitializeSyncIndex2(const int& thread_id)
	{
		if(thread_id == 0) sync_index2 = 0;

		Sync(thread_id);
	}

	void GetIterationIndex(int& index, const int& increase)
	{
		boost::mutex::scoped_lock lock(sync_mutex);

		index = sync_index;

		sync_index += increase;
	}

	void GetIterationIndex2(int& index, const int& increase)
	{
		boost::mutex::scoped_lock lock(sync_mutex);

		index = sync_index2;

		sync_index2 += increase;
	}

	void Sleep()
	{
		boost::this_thread::sleep(boost::posix_time::microseconds(1));
	}

	void SyncMax(const int& thread_id, T& value)
	{
		Sync(thread_id);

		sync_value[thread_id] = value;

		Sync(thread_id);

		for(int i = 0; i < num_threads; i++)
		{
			if(value < sync_value[i]) 
			{
				value = sync_value[i];
			}
		}

		Sync(thread_id);
	}

	void SyncMax(const int& thread_id, int& value)
	{
		Sync(thread_id);

		sync_value_int[thread_id] = value;

		Sync(thread_id);

		for (int i = 0; i < num_threads; i++)
		{
			if (value < sync_value_int[i])
			{
				value = sync_value_int[i];
			}
		}

		Sync(thread_id);
	}

	void SyncMin(const int& thread_id, T& value)
	{
		Sync(thread_id);

		sync_value[thread_id] = value;

		Sync(thread_id);
	
		for(int i = 0; i < num_threads; i++)
		{
			if(value > sync_value[i]) 
			{
				value = sync_value[i];
			}
		}

		Sync(thread_id);
	}

	void SyncSum(const int& thread_id, T& value) 
	{
		Sync(thread_id);

		sync_value[thread_id] = value;

		Sync(thread_id);

		value = sync_value[0];			// To remove one assignment
		for(int i = 1; i < num_threads; i++)
		{
			sync_value_temp += sync_value[i];
		}
				
		Sync(thread_id);
	}

	void SyncSum(const int& thread_id, T& value, const int& start_index, const int& end_index)
	{
		Sync(thread_id);

		sync_value[thread_id] = value;

		Sync(thread_id);

		assert(start_index >= 0);
		assert(end_index <= num_threads);
		assert(start_index <= end_index);

		value = sync_value[start_index];
		
		for(int i = start_index + 1 ; i <= end_index; i++)
		{
			value += sync_value[i];
		}

		Sync(thread_id);
	}

	void SyncSum(const int& thread_id, int& value)
	{
		Sync(thread_id);

		sync_value_int[thread_id] = value;

		Sync(thread_id);

		value = sync_value_int[0];			// To remove one assignment
		for (int i = 1; i < num_threads; i++)
		{
			value += sync_value_int[i];
		}

		Sync(thread_id);
	}

	int SyncSumFromArray(const int& thread_id, const int* value_array)
	{
		Sync(thread_id);

		for(int i = 0; i < num_threads; i++)
		{
			sync_value_int_array[i + num_threads*thread_id] = value_array[i];
		}

		Sync(thread_id);

		int sum(sync_value_int_array[thread_id]);

		for(int i = 1; i < num_threads; i++)
		{
			sum += sync_value_int_array[thread_id + num_threads*i];
		}

		Sync(thread_id);

		return sum;
	}

	void SyncDomainIndices1D(const int& thread_id, const int& res, int& i_start, int& i_end)
	{
		sync_value_int[thread_id] = res;

		Sync(thread_id);

		i_start = 0;

		for(int i = 0; i <= thread_id - 1; i++)
		{
			i_start += sync_value_int[i];
		}

		i_end = i_start + res - 1;

		Sync(thread_id);
	}
	
	void JoinAll()
	{
		for(int i = 0; i < num_threads; i++)
		{
			thread_list[i]->join();
		}

		DeleteAllThreads();
	}

	void Reset()
	{
		num_of_waiting_threads = 0;

		if(sync_value != 0)
		{
			delete [] sync_value;
			sync_value = 0;
		}

		if(thread_list != 0)
		{
			for(int i = 0; i < num_threads; i++)
			{
				if(thread_list[i] != 0)
				{
					delete thread_list[i];
					thread_list[i] = 0;
				}
			}
			delete [] thread_list;
			thread_list = 0;
		}
	}

	void DeleteAllThreads()
	{
		for(int i = 0; i < num_threads; i++)
		{
			if(thread_list[i] != 0)
			{
				delete thread_list[i];
				thread_list[i] = 0;
			}
		}

		num_of_waiting_threads = 0;
	}

	void SplitDomainIndex1D(const int& k_start, const int& k_res)
	{
		const int k_end = k_start + k_res - 1;
		const int quotient = k_res / num_threads;
		const int remainder = k_res % num_threads;

		int k_start_p = k_start;

		for(int i = 0; i < num_threads; i++)
		{
			int k_depth = i < remainder ? (quotient + 1) : quotient;

			start_ix_1D[i] = k_start_p;
			end_ix_1D[i] = k_start_p + k_depth - 1;

			k_start_p += k_depth;
		}
	}

	void SplitDomainIndex1D(const int& thread_id, const int& k_start, const int& k_res)
	{
		if(thread_id == 0) SplitDomainIndex1D(k_start, k_res);

		Sync(thread_id);
	}

	void SplitDomainIndex1D(const int& thread_id, const int& k_start, const int& k_res, int& ix_start, int& ix_end)
	{
		if(thread_id == 0) SplitDomainIndex1D(k_start, k_res);

		Sync(thread_id);

		ix_start = start_ix_1D[thread_id];
		ix_end = end_ix_1D[thread_id];
	}

	// Random functions are implemented as a part of MULTITHREADING so that they can use seed numbers (gen[thread_id]) per thread
	inline T RandomNumber(const int& thread_id)
	{
		boost::random::uniform_int_distribution<> dist(0, 1000000);
		return (T)dist(random_number_generator_per_thread[thread_id])/(T)1000000;
	}

	inline VT RandomUnitVector(const int& thread_id)
	{
		// Spherical Coordinate
		const T theta = RandomNumber(thread_id)*(T)3.1415;
		const T phi = RandomNumber(thread_id)*(T)6.2830;

		return VT(sin(theta)*cos(phi), sin(theta)*sin(phi));
	}

	inline VT RandomVectorInUnitSphere(const int& thread_id)
	{
		// Spherical Coordinate
		const T r = RandomNumber(thread_id);
		const T theta = RandomNumber(thread_id)*(T)3.1415;
		const T phi = RandomNumber(thread_id)*(T)6.2830;

		return VT(sin(theta)*cos(phi)*r, sin(theta)*sin(phi)*r);
	}

	inline VT RandomVector(const int& thread_id)
	{
		return VT(RandomNumber(thread_id) - (T)0.5, RandomNumber(thread_id) - (T)0.5);
	}

	template<class F, class A1>
	void RunThreads(F f, A1 a1)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id);
		}
		JoinAll();
	}

	template<class F, class A1, class A2>
	void RunThreads(F f, A1 a1, A2 a2)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id, a2);
		}
		JoinAll();
	}

	template<class F, class A1, class A2, class A3>
	void RunThreads(F f, A1 a1, A2 a2, A3 a3)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id, a2, a3);
		}
		JoinAll();
	}

	template<class F, class A1, class A2, class A3, class A4>
	void RunThreads(F f, A1 a1, A2 a2, A3 a3, A4 a4)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id, a2, a3, a4);
		}
		JoinAll();
	}

	template<class F, class A1, class A2, class A3, class A4, class A5>
	void RunThreads(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id, a2, a3, a4, a5);
		}
		JoinAll();
	}

	template<class F, class A1, class A2, class A3, class A4, class A5, class A6>
	void RunThreads(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id, a2, a3, a4, a5, a6);
		}
		JoinAll();
	}

	template<class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7>
	void RunThreads(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id, a2, a3, a4, a5, a6, a7);
		}
		JoinAll();
	}

	template<class F, class A1, class A2, class A3, class A4, class A5, class A6, class A7, class A8>
	void RunThreads(F f, A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7, A8 a8)
	{
		for(int thread_id = 0; thread_id < num_threads; thread_id++)
		{
			thread_list[thread_id] = new MULTITHREADING::THREAD(f, a1, thread_id, a2, a3, a4, a5, a6, a7, a8);
		}
		JoinAll();
	}
};


			








