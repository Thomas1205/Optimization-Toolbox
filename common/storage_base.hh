/*-*-c++-*-*/
/*** written by Thomas Schoenemann as a private person, November 2019 ***/

#ifndef STORAGE_BASE_HH
#define STORAGE_BASE_HH

#include "makros.hh"
#include <initializer_list>

//operators == and != need to be defined on T
template<typename T, typename ST=size_t>
class StorageBase {
public:

  //for compatibility with the STL:
  using value_type = T;

  StorageBase();

  StorageBase(ST size);

  StorageBase(ST size, const T default_value);

  //copy constructor
  StorageBase(const StorageBase<T,ST>& toCopy);

  //move constructor
  StorageBase(StorageBase<T,ST>&& toTake);

  ~StorageBase();

  virtual const std::string& name() const;

  inline ST size() const;

  inline T* attr_restrict direct_access();

  inline const T* attr_restrict direct_access() const;

  inline T& ref_attr_restrict direct_access(ST i);

  inline const T& ref_attr_restrict direct_access(ST i) const;

  inline void set_constant(const T constant);

  inline bool is_constant() const;

protected:

  StorageBase(const std::initializer_list<T>& init);
  
  StorageBase<T,ST>& operator=(const StorageBase<T,ST>& toCopy);

  StorageBase<T,ST>& operator=(StorageBase<T,ST>&& toTake);

  //pointer must go first so that following variables can be grouped for optimal alignment
  T* data_; //if `T is a basic type this is 16-byte aligned as returned by new
  ST size_;
  static const std::string stor_base_name_;
};

/*************************** implementation ****************************/

template<typename T, typename ST> 
StorageBase<T,ST>::StorageBase() : data_(0), size_(0)
{
}

template<typename T, typename ST> 
StorageBase<T,ST>::StorageBase(ST size) : size_(size)
{
  data_ = new T[size];
}

template<typename T, typename ST> 
StorageBase<T,ST>::StorageBase(ST size, const T default_value) : size_(size)
{
  data_ = new T[size];
  std::fill_n(data_,size_,default_value);
}

template<typename T, typename ST> 
StorageBase<T,ST>::StorageBase(const std::initializer_list<T>& init)
{
  size_ = init.size();
  std::copy(init.begin(),init.end(),data_);
}

//copy constructor
template<typename T, typename ST> 
StorageBase<T,ST>::StorageBase(const StorageBase<T,ST>& toCopy) 
{
  size_ = toCopy.size();
  data_ = new T[size_];
  
  Makros::unified_assign(data_, toCopy.direct_access(), size_);
}

//move constructor
template<typename T, typename ST> 
StorageBase<T,ST>::StorageBase(StorageBase<T,ST>&& toTake) : data_(toTake.data_), size_(toTake.size_) 
{
  toTake.data_ = 0;
}

template<typename T, typename ST> 
StorageBase<T,ST>::~StorageBase()
{
  delete[] data_;
}

template<typename T, typename ST> 
StorageBase<T,ST>& StorageBase<T,ST>::operator=(const StorageBase<T,ST>& toCopy)
{
  if (size_ != toCopy.size_) {

    delete[] data_;

    size_ = toCopy.size_;
    data_ = new T[size_];
  }

  Makros::unified_assign(data_, toCopy.direct_access(), size_);

  return *this;
}

template<typename T, typename ST> 
StorageBase<T,ST>& StorageBase<T,ST>::operator=(StorageBase<T,ST>&& toTake)
{
  delete[] data_;
  data_ = toTake.data_;
  size_ = toTake.size_;
  toTake.data_ = 0;

  return *this;
}

template<typename T, typename ST>
/*virtual*/ const std::string& StorageBase<T,ST>::name() const
{
  return stor_base_name_;
}

template<typename T, typename ST>
inline ST StorageBase<T,ST>::size() const
{
  return size_;
}

template<typename T, typename ST>
inline T* attr_restrict StorageBase<T,ST>::direct_access()
{
  return data_;
}

template<typename T, typename ST>
inline const T* attr_restrict StorageBase<T,ST>::direct_access() const
{
  return data_;
}

template<typename T, typename ST>
inline T& ref_attr_restrict StorageBase<T,ST>::direct_access(ST i)
{
  return data_[i];
}

template<typename T, typename ST>
inline const T& ref_attr_restrict StorageBase<T,ST>::direct_access(ST i) const
{
  return data_[i];
}

template<typename T, typename ST>
inline void StorageBase<T,ST>::set_constant(const T constant)
{
  std::fill_n(data_,size_,constant);
}

template<typename T, typename ST>
inline bool StorageBase<T,ST>::is_constant() const
{
  if (size_ <= 1)
    return true;
  
  const T val = data_[0];
  for (ST i = 1; i < size_; i++) {
    if (val != data_[i])
      return false;
  }
  return true;
}

template<typename T,typename ST>
/*static*/ const std::string StorageBase<T,ST>::stor_base_name_ = "unnamed base storage";


#endif