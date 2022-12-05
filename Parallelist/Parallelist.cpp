#include <random>
#include <iterator>
#include <algorithm>
#include <vector>
#include <chrono>
#include <atomic>

using namespace std;

template <class T>
struct CNode {
    T valor;
    CNode<T>* next;
    CNode(T _v) { valor = _v; next = 0; }

    mutex node_lock;
    bool marked = false;
    void lock() {
        node_lock.lock();
    }
    void unlock() {
        node_lock.unlock();
    }
};

template <class T>
class CList {

public:
    CList();
    bool insert(T x);
    bool remove(T x);
    void print();
    bool validate(CNode<T>* pred, CNode<T>* curr);
private:
    CNode<T>* head;
    CNode<T>* tail;
};

template <class T>
CList<T>::CList() {
    tail = new CNode<T>(99999999);
    head = new CNode<T>(-99999999);
    head->next = tail;
}


template <class T>
bool CList<T>::insert(T x) {
    bool found = false;
    while (true) {
        CNode<T>* pred = head;
        CNode<T>* curr = head->next;
        while (curr->valor < x) {
            pred = curr;
            curr = curr->next;
        }
        pred->lock();
        curr->lock();

        if (validate(pred, curr)) {
            if (curr->valor == x) {
                found = false;
            }
            else {
                CNode<T>* new_node = new CNode<T>(x);
                new_node->next = curr;
                pred->next = new_node;
                found = true;
            }
            pred->unlock();
            curr->unlock();
            return found;
        }
        pred->unlock();
        curr->unlock();
    }
}

template <class T>
bool CList<T>::remove(T x) {
    bool found = false;
    while (true) {
        CNode<T>* pred = head;
        CNode<T>* curr = head->next;
        while (curr->valor < x) {
            pred = curr;
            curr = curr->next;
        }
        pred->lock();
        curr->lock();
        if (validate(pred, curr)) {
            if (curr->valor != x) {
                curr->unlock();
                found = false;
            }
            else {
                curr->marked = true;
                pred->next = curr->next;
                curr->unlock();
                found = true;
            }
            pred->unlock();
            return found;
        }

        pred->unlock();
        curr->unlock();
    }
}

template <class T>
bool CList<T>::validate(CNode<T>* pred, CNode<T>* curr) {

    return !pred->marked && !curr->marked && pred->next == curr;
}


template <class T>
void CList<T>::print() {
    cout << "H->";
    for (CNode<T>* a = head->next; a && a != tail; a = a->next)
        cout << a->valor << "->";


}

int op = 100;

int get_random(int low, int high) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distribution(low, high);
    return distribution(gen);
}


struct Insertar {

    int min_;
    int max_;
    CList<int>* ptr_;

    Insertar(int min, int max, CList<int>* ptr) :min_(min), max_(max), ptr_(ptr) {}
    void operator()(int operaciones) {

        for (int i = 0; i < operaciones; i++) {
            int x = get_random(min_, max_);
            cout << "Insertando el numero: " << x << endl;
            ptr_->insert(x);
        }
    }
};


struct Eliminar {

    int min_;
    int max_;
    CList<int>* ptr_;

    Eliminar(int min, int max, CList<int>* ptr) :min_(min), max_(max), ptr_(ptr) {}
    void operator()(int operaciones) {

        for (int i = 0; i < operaciones; i++) {
            int x = get_random(min_, max_);
            cout << "Eliminando el numero: " << x << endl;
            ptr_->remove(x);
        }
    }
};




int main()
{

    CList<int> list;
    vector<thread>insertar;
    vector<thread>eliminar;

    int num_threads = 8;
    for (int i = 0; i < num_threads; i++) {


        insertar.push_back(thread(Insertar(1, 100, &list), op));
        eliminar.push_back(thread(Eliminar(1, 100, &list), op));

    }

    for (int i = 0; i < num_threads; i++) {
        insertar[i].join();
        eliminar[i].join();

    }

    return 0;
}
