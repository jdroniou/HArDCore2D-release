#ifndef PARALLEL_FOR
#define PARALLEL_FOR

#include <thread>

namespace HArDCore2D
{
  
  /// Parallel for
  static void parallel_for(unsigned nb_elements,
			   std::function<void(size_t start, size_t end)> functor,
			   bool use_threads = true)
  {
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);

    unsigned batch_size = nb_elements / nb_threads;
    unsigned batch_remainder = nb_elements % nb_threads;

    std::vector<std::thread> my_threads(nb_threads);

    if(use_threads) {
      // Multithread execution
      for(unsigned i = 0; i < nb_threads; ++i) {
      	int start = i * batch_size;
      	my_threads[i] = std::thread(functor, start, start + batch_size);
      }
    } else {
      // Single thread execution (for easy debugging)
      for(unsigned i = 0; i < nb_threads; ++i) {
      	int start = i * batch_size;
      	functor(start, start + batch_size);
      }
    }

    // Execute the elements left
    int start = nb_threads * batch_size;
    functor(start, start + batch_remainder);

    // Wait for the other thread to finish their task
    if (use_threads) {
      std::for_each(my_threads.begin(), my_threads.end(), std::mem_fn(&std::thread::join));
    }
  }
  
} // end of namespace HArDCore2D
#endif
