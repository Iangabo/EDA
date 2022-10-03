// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"


template <size_t N, typename ElemType>
class KDTree {
 public:
  typedef std::pair<Point<N>, ElemType> value_type;

  typedef std::set<std::pair<double, ElemType>> BoundedPQueue;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value);

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

  ElemType knn_value(const Point<N> &key, size_t k) const;

  std::vector<ElemType> knn_query(const Point<N> &key, size_t k) const;

  std::vector<ElemType> ortogonal_range_query(const std::vector<std::pair<double, double>>& ranges) const;

private:
    struct KDTreeNode
    {
    public:
        value_type point_;
        KDTreeNode* left_;
        KDTreeNode* right_;

        KDTreeNode(const value_type& value) :point_(value), left_(nullptr), right_(nullptr) {}
    };
    mutable KDTreeNode* root_ = nullptr;
    size_t dimension_;
    size_t size_ = 0;
    ElemType temporal;
    bool find(const Point<N>& pt, KDTreeNode**& it) const;
    void find_neighbors_(KDTreeNode* current_node, BoundedPQueue &nearest_neighbors_candidates, int depth, const Point<N> key) const;
    KDTreeNode* copiar(KDTreeNode* root_);
    
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  // TODO(me): Fill this in.
    dimension_ = N;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  // TODO(me): Fill this in.
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  // TODO(me): Fill this in.
    this->dimension_ = rhs.dimension_;
    this->size_ = rhs.size_;
    root_ = copiar(rhs.root_);

}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  // TODO(me): Fill this in.

  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  // TODO(me): Fill this in.
  return dimension_;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  // TODO(me): Fill this in.
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    if (size_ == 0) {
        return true;
    }
    else {
        return false;
    }
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
    KDTreeNode** it;
    return find(pt, it);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N> &pt, const ElemType &value) {
  // TODO(me): Fill this in.
    value_type point(pt,value);
    KDTreeNode** it;

    if (find(pt, it)) {
        (*it)->point_.second = value;
    }
    else {
        (*it) = new KDTreeNode(point);
        size_++;
    }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::operator[](const Point<N> &pt) {
    
    value_type point(pt, NULL);
    KDTreeNode** it;
    
    if (find(pt, it)) {
        return (*it)->point_.second;
    }
    else {
        size_++;
        (*it) = new KDTreeNode(point);
        return (*it)->point_.second;
    }
}

template <size_t N, typename ElemType>
ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) {
  // TODO(me): Fill this in.
    KDTreeNode** it;
    if (find(pt, it)) {
        return (*it)->point_.second;
    }
    else {
        throw std::out_of_range("Fuera de rango");
    }
}

template <size_t N, typename ElemType>
const ElemType &KDTree<N, ElemType>::at(const Point<N> &pt) const {
  // TODO(me): Fill this in.
    KDTreeNode** it;
    if (find(pt, it)) {
        return (*it)->point_.second;
    }
    else {
        throw std::out_of_range("Fuera de rango");
    }
}

template <size_t N, typename ElemType>

ElemType KDTree<N, ElemType>::knn_value(const Point<N> &key, size_t k) const {
  // TODO(me): Fill this in.
    KDTreeNode* it = root_;
    BoundedPQueue nearest_neighbors_candidates;
    int depth = 0;
    
    find_neighbors_(it,nearest_neighbors_candidates,depth, key);

    std::vector<std::pair<double, ElemType>> aux(nearest_neighbors_candidates.begin(), nearest_neighbors_candidates.end());
    double val = (*nearest_neighbors_candidates.begin()).first;
    int cont = 0;

    for (int i = 0; i < nearest_neighbors_candidates.size(); i++) {
        if (aux[i].first == val) {
            cont++;
        }
        else {
            break;
        }
    }

    return aux[cont-1].second;
   
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N> &key,
                                                     size_t k) const {
  // TODO(me): Fill this in.
  std::vector<ElemType> values;
  KDTreeNode* it = root_;
  BoundedPQueue nearest_neighbors_candidates;
  int depth = 0;

  find_neighbors_(it, nearest_neighbors_candidates, depth, key);

  std::vector<std::pair<ElemType, ElemType>> aux(nearest_neighbors_candidates.begin(), nearest_neighbors_candidates.end());
  ElemType val = (*nearest_neighbors_candidates.begin()).first;
  int cont = 0;
  for (size_t i = 0; i < k, i++) {
      values[i] = aux[i].second;
  }

  return values;
}

template <size_t N, typename ElemType>
bool KDTree<N,ElemType>::find(const Point<N>& pt, KDTreeNode**& it) const {
    it = &root_;
    int aux = 0;

 
    while ((*it) && (*it)->point_.first != pt) {

        if ((*it)->point_.first[aux % dimension_] > pt[aux % dimension_]) {
            it = &((*it)->left_);
        }
        else {
            it = &((*it)->right_);
        }
        aux++;
    }

    return *it != NULL;

}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::find_neighbors_(KDTreeNode* current_node, BoundedPQueue &nearest_neighbors_candidates, int depth, const Point<N> key) const {
    if (!current_node) {
        return;
    }   
    nearest_neighbors_candidates.insert(std::make_pair(distance(current_node->point_.first, key), current_node->point_.second));
    
    int axis = depth % N;
    bool right = false;
    if (key[axis] < current_node->point_.first[axis]) {
        right = false;
        find_neighbors_(current_node->left_, nearest_neighbors_candidates, ++depth, key);
    }
    else {
        right = true;
        find_neighbors_(current_node->right_, nearest_neighbors_candidates, ++depth, key);
    }

    if (fabs(current_node->point_.first[axis] - key[axis]) < (*nearest_neighbors_candidates.begin()).first) {
        if (right) {
            find_neighbors_(current_node->left_, nearest_neighbors_candidates, ++depth, key);
        }
        else {
            find_neighbors_(current_node->right_, nearest_neighbors_candidates, ++depth, key);
        }
    }
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDTReeNode* copiar(typename KDTree<N, ElemType>::KDTreeNode* root) {
    if (root == NULL) return NULL;
    KDTreeNode* newroot = new KDTreeNode(*root);
    newroot->left_ = copiar(root->left_);
    newroot->right_ = copiar(root->right_);
    
    return newRoot;
}
// TODO(me): finish the implementation of the rest of the KDTree class

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::ortogonal_range_query(const std::vector<std::pair<double, double>>& ranges) const{
    //Vector Final
    std::vector<ElemType> f;
    
    //Recorrer los puntos mas cercanos a cada esquina del rectangulo
    for (int i = 0; i < reanges.size()) {
        std::vector<ElemType> knn = knn_query(ranges[i]);
        for (int i = 0; i < knn.size(); i++) {
            //Si esta dentro del rectangulo entonces añadelo a f
            if (in(knn, ranges)) {
                f.push_back(knn[i]);
            }
        }
    }

    //retorna f
    return f;

}

#endif  // SRC_KDTREE_HPP_
