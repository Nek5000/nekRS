#include <iostream>
#include <vector>
#include <thread>
#include <functional>
#include <future>
#include <queue>
#include <mutex>
#include <condition_variable>

class ThreadPool
{
public:
  ThreadPool(size_t numThreads);
  ~ThreadPool();

  template <class F, class... Args>
  auto enqueue(F &&f, Args &&...args) -> std::future<typename std::invoke_result<F, Args...>::type>;

  void finish();

  void safePrint(const std::string &message)
  {
    std::lock_guard<std::mutex> guard(coutMutex);
    std::cout << message << std::endl;
  }

private:
  std::vector<std::thread> workers;
  std::queue<std::function<void()>> tasks;
  std::mutex queueMutex;
  std::condition_variable condition;
  bool stop;

  std::mutex coutMutex;

  std::mutex taskMutex;
  std::condition_variable taskCondition;
  size_t activeTasks;

  void workerThread();
};

template <class F, class... Args>
auto ThreadPool::enqueue(F &&f, Args &&...args) -> std::future<typename std::invoke_result<F, Args...>::type>
{
  using returnType = typename std::invoke_result<F, Args...>::type;

  auto task = std::make_shared<std::packaged_task<returnType()>>(
      std::bind(std::forward<F>(f), std::forward<Args>(args)...));
  std::future<returnType> result = task->get_future();

  {
    std::unique_lock<std::mutex> lock(queueMutex);
    if (stop) {
      throw std::runtime_error("Enqueue on stopped ThreadPool");
    }
    tasks.emplace([task, this] {
      {
        std::unique_lock<std::mutex> taskLock(taskMutex);
        ++activeTasks;
      }
      (*task)();
      {
        std::unique_lock<std::mutex> taskLock(taskMutex);
        --activeTasks;
        if (activeTasks == 0) {
          taskCondition.notify_all();
        }
      }
    });
  }

  condition.notify_one();
  return result;
}
