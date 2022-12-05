#include <iostream>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <queue>
#include <vector>
#include <memory>
#include <random>
#include <functional>
#include <string>
#include <sstream>

using namespace std;

template <typename T>
class ParallelQueue {
public:
    void push(const T& data);       
    T pop();                       
private:
    queue<T> queue_;
    mutex queue_mutex_;
    condition_variable condition1;    
};


template <typename T>
void ParallelQueue<T>::push(const T& data) {

    unique_lock<mutex>lock(queue_mutex_);
    queue_.push(data);
    condition1.notify_one();   
}

template <typename T>
T ParallelQueue<T>::pop(){

    unique_lock<mutex>lock(queue_mutex_); 
    while (queue_.empty()) {                        
        condition1.wait(lock);                       
    }
    T result = queue_.front();
    queue_.pop();
    return result;
}


class Producer {
public:
    Producer(unsigned int id, ParallelQueue<string>* queue)
        : id_(id), queue_(queue) {}

    void operator()() {
        int data = 0;
        while (true) {
           stringstream stream;
            stream << "Producer: " << id_ << " Data: " << data++ << endl;
            queue_->push(stream.str());
            cout << stream.str() << endl;
        }
    }

private:
    unsigned int id_;
    ParallelQueue<string>* queue_;
};

class Consumer {
public:
    Consumer(unsigned int id, ParallelQueue<string>* queue)
        : id_(id), queue_(queue) {}

    void operator()() {
        while (true) {
            stringstream stream;
            stream << "Consumer: " << id_ << " Data: " << queue_->pop().c_str()
                << endl;

           cout << stream.str() << endl;
        }
    }

private:
    unsigned int id_;
    ParallelQueue<string>* queue_;
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        return 1;
    }
    int number_producers = stoi(argv[1]);
    int number_consumers = stoi(argv[2]);

    ParallelQueue<string> queue;

    vector<thread*> producers;
    for (unsigned int i = 0; i < number_producers; ++i) {
        thread* producer_thread = new thread(Producer(i, &queue));
        producers.push_back(producer_thread);
    }

    vector<thread*> consumers;
    for (unsigned int i = 0; i < number_consumers; ++i) {
        thread* consumer_thread = new thread(Consumer(i, &queue));
        consumers.push_back(consumer_thread);
    }

    int stop;
    cin >> stop;

    return 0;
}