#include <iostream>
#include <thread>
#include <mutex>
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
        pred->lock_();  
        curr->lock_();  

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
                curr->unlock_();    
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
