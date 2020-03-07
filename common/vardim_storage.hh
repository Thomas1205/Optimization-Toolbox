/**** written by Thomas Schoenemann as a private person without employment, July 2011 ****/

#ifndef VAR_DIM_STORAGE_HH
#define VAR_DIM_STORAGE_HH

#include "vector.hh"

template<typename T>
class VarDimStorage : public StorageBase<T> {
public:

  using Base = StorageBase<T>;

  explicit VarDimStorage(const Math1D::Vector<size_t>& dim);

  explicit VarDimStorage(const Math1D::Vector<size_t>& dim, T fill);

  //copy constructor
  VarDimStorage(const VarDimStorage& toCopy);

  //move constructor
  VarDimStorage(VarDimStorage&& toTake);

  ~VarDimStorage();

  size_t dim(uint n) const;

  const Math1D::Vector<size_t>& dims();

  size_t nDims() const;

  const Math1D::Vector<size_t>& dim_vector() const;

  const T& operator()(Math1D::Vector<size_t>& pos) const;

  T& operator()(Math1D::Vector<size_t>& pos);

  T data(uint pos) const;

  VarDimStorage<T>& operator=(const VarDimStorage& toCopy);

  VarDimStorage<T>& operator=(VarDimStorage&& toTake) = default;

protected:

  Math1D::Vector<size_t> dim_;
};

template<typename T, typename ST>
bool operator==(const VarDimStorage<T>& v1, const VarDimStorage<T>& v2);

template<typename T, typename ST>
bool operator!=(const VarDimStorage<T>& v1, const VarDimStorage<T>& v2);


/*********** implementation *******/


template<typename T> 
VarDimStorage<T>::VarDimStorage(const Math1D::Vector<size_t>& dim) : StorageBase<T>(), dim_(dim)
{
  Base::size_ = (dim.size() == 0) ? 0 : 1;

  for (uint k=0; k < dim.size(); k++)
    Base::size_ *= dim[k];

  Base::data_ = new T[Base::size_];
}

template<typename T> 
VarDimStorage<T>::VarDimStorage(const Math1D::Vector<size_t>& dim, T fill) : StorageBase<T>(), dim_(dim)
{
  Base::size_ = (dim.size() == 0) ? 0 : 1;

  for (size_t k=0; k < dim.size(); k++)
    Base::size_ *= dim[k];

  Base::data_ = new T[Base::size_];

  std::fill_n(Base::data_,Base::data_+Base::size_,fill);
}

template<typename T> 
VarDimStorage<T>::VarDimStorage(const VarDimStorage& toCopy)
{
  Base::size_ = toCopy.size();
  dim_ = toCopy.dim_;

  Base::data_ = new T[Base::size_];
  for (uint k=0; k < Base::size_; k++) {
    Base::data_[k] = toCopy.data(k);
  }
}

//move constructor
template<typename T> 
VarDimStorage<T>::VarDimStorage(VarDimStorage&& toTake)
{
  Base::data_ = toTake.data_;
  Base::size_ = toTake.size_;
  dim_ = std::move(toTake.dim_);
  
  toTake.data_ = 0;
}

template<typename T> 
VarDimStorage<T>::~VarDimStorage()
{
}

template<typename T>
VarDimStorage<T>& VarDimStorage<T>::operator=(const VarDimStorage& toCopy)
{
  Base::size_ = toCopy.size();
  dim_ = toCopy.dim_vector();

  Base::data_ = new T[Base::size_];
  for (uint k=0; k < Base::size_; k++) {
    Base::data_[k] = toCopy.data(k);
  }
  
  return *this;
}

template<typename T>
size_t VarDimStorage<T>::dim(uint n) const
{
  return dim_[n];
}

template<typename T>
size_t VarDimStorage<T>::nDims() const
{
  return dim_.size();
}

template<typename T>
const Math1D::Vector<size_t>& VarDimStorage<T>::dims()
{
  return dim_;
}

template<typename T>
T VarDimStorage<T>::data(uint pos) const
{
  return Base::data_[pos];
}

template<typename T>
const Math1D::Vector<size_t>& VarDimStorage<T>::dim_vector() const
{
  return dim_;
}

template<typename T>
const T& VarDimStorage<T>::operator()(Math1D::Vector<size_t>& pos) const
{
  assert(pos.size() == dim_.size());

  uint data_pos = 0;

  for (size_t k=0; k < dim_.size(); k++) {

    assert(pos[k] < dim_[k]);

    if (k > 0)
      data_pos *= dim_[k];

    data_pos += pos[k];
  }

  return Base::data_[data_pos];
}

template<typename T>
T& VarDimStorage<T>::operator()(Math1D::Vector<size_t>& pos)
{
  assert(pos.size() == dim_.size());

  uint data_pos = 0;

  for (size_t k=0; k < dim_.size(); k++) {

    assert(pos[k] < dim_[k]);

    if (k > 0)
      data_pos *= dim_[k];

    data_pos += pos[k];
  }

  return Base::data_[data_pos];
}

template<typename T, typename ST>
bool operator==(const VarDimStorage<T>& v1, const VarDimStorage<T>& v2)
{
  if (v1.dims() != v2.dims())
    return false;

  for (size_t k = 0; k < v1.size(); k++) {
    if (v1.direct_access(k) != v2.direct_access(k))
      return false;
  }

  return true;
}

template<typename T, typename ST>
bool operator!=(const VarDimStorage<T>& v1, const VarDimStorage<T>& v2)
{
  return !(v1 == v2);
}




#endif
