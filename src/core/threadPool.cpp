#include <threadPool.hpp>

ThreadPool::ThreadPool(size_t numThreads) : stop(false), activeTasks(0)
{
  for (size_t i = 0; i < numThreads; ++i) {
    workers.emplace_back([this] { this->workerThread(); });
  }
}

ThreadPool::~ThreadPool()
{
  {
    std::unique_lock<std::mutex> lock(queueMutex);
    stop = true;
  }
  condition.notify_all();
  for (std::thread &worker : workers) {
    worker.join();
  }
}

void ThreadPool::workerThread()
{
  while (true) {
    std::function<void()> task;
    {
      std::unique_lock<std::mutex> lock(queueMutex);
      condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
      if (this->stop && this->tasks.empty()) {
        return;
      }
      task = std::move(this->tasks.front());
      this->tasks.pop();
    }
    task();
  }
}

void ThreadPool::finish()
{
  std::unique_lock<std::mutex> lock(queueMutex);
  taskCondition.wait(lock, [this] { return tasks.empty() && activeTasks == 0; });
}
