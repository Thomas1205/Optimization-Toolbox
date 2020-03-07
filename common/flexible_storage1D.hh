/********* Written by Thomas Schoenemann, created from moved code, February 2020 **********/

#ifndef FLEXIBLE_STORAGE1D_HH
#define FLEXIBLE_STORAGE1D_HH

#include "makros.hh"
#include "storage_base.hh"
#include <type_traits>

#ifndef INIT_CAPACITY
#define INIT_CAPACITY 4
#endif

//this class is meant to replace std::vector with its push_back() functionality.
// It has slightly less functionality, though. E.g. erase() is not available.
template<typename T, typename ST=size_t>
class FlexibleStorage1D : public StorageBase<T,ST> {
public:

  using Base = StorageBase<T,ST>;

  using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

  explicit FlexibleStorage1D();

  explicit FlexibleStorage1D(ST reserved_size);

  explicit FlexibleStorage1D(ST reserved_size, ST fill_size, T val);

  explicit FlexibleStorage1D(const std::initializer_list<T>& init, ST extra_reserve = 0);

  //copy constructor
  FlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy);
  
  //move constructor
  FlexibleStorage1D(FlexibleStorage1D<T,ST>&& toTake);

  ~FlexibleStorage1D();

  inline T& operator[](ST i) const;

  inline ST reserved_size() const;

  //for compatibility with std::vector, e.g. for use in templates
  inline ST capacity() const;

  void resize(ST size, bool exact_fit = false);

  void fit_exactly();

  void shrink(ST size);

  inline void shrink_by(ST reduction);
  
  inline void grow_by_dirty(ST increase);

  void reserve(ST size);

  //will not free memory
  void clear();

  FlexibleStorage1D<T,ST>& operator=(const FlexibleStorage1D<T,ST>& toCopy);

  inline T back() const;

  inline T& back();

  ST append(PassType val);
    
  ST move_append(T&& val);

  //shortcut when you are sure the allocated memory suffices
  inline void append_trusting(PassType val);

  //for compatibility with std::vector, e.g. for use in templates
  inline void push_back(PassType val);

  void append(Storage1D<T,ST>& toAppend);

  void append(FlexibleStorage1D<T,ST>& toAppend);

  void erase(ST pos);
  
  void erase_several(ST pos, ST nToErase);

  void insert(ST pos, T val);

  void insert_several(ST pos, T* vals, ST nData);

  void swap(FlexibleStorage1D<T,ST>& toSwap);

protected:

  ST reserved_size_;
  static const std::string flex_stor1D_name_;
};

template<typename T, typename ST>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T,ST>& v);

template<typename T, typename ST>
bool operator==(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2);

template<typename T, typename ST>
bool operator!=(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2);

template<typename T, typename ST=size_t>
class NamedFlexibleStorage1D : public FlexibleStorage1D<T,ST> {
public:

  NamedFlexibleStorage1D();

  NamedFlexibleStorage1D(const std::string& name);

  NamedFlexibleStorage1D(ST reserved_size, const std::string& name);

  //copy constructors
  NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T,ST>& toCopy);

  NamedFlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy);

  virtual const std::string& name() const;

  //operators
  void operator=(const NamedFlexibleStorage1D<T,ST>& toCopy);

  void operator=(const FlexibleStorage1D<T,ST>& toCopy);

protected:
  std::string name_;
};

/******* implementation of FlexibleStorage1D *********/

template<typename T, typename ST>
/*static*/ const std::string FlexibleStorage1D<T,ST>::flex_stor1D_name_ = "unnamed flexible 1Dstorage";

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D() : StorageBase<T,ST>(0), reserved_size_(INIT_CAPACITY)
{
#if INIT_CAPACITY > 0
  Base::data_ = new T[reserved_size_];
#endif  
}

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D(ST reserved_size) : StorageBase<T,ST>(0), reserved_size_(reserved_size)
{
  Base::data_ = new T[reserved_size_];
}

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D(ST reserved_size, ST fill_size, T val) : StorageBase<T,ST>(0), reserved_size_(reserved_size)
{
  assert(fill_size <= reserved_size_);
  
  Base::data_ = new T[reserved_size_];
  Base::size_ = fill_size;
  std::fill_n(Base::data_,fill_size,val);
}

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D(const std::initializer_list<T>& init, ST extra_reserve) : StorageBase<T,ST>(0)
{
  reserved_size_ = init.size() + extra_reserve;
  Base::size_ = init.size();
  Base::data_ = new T[reserved_size_];
  
  std::copy(init.begin(),init.end(),Base::data_);
}

//copy constructor
template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy)
{
  Base::size_ = toCopy.size();
  reserved_size_ = toCopy.reserved_size();

  Base::data_ = new T[reserved_size_];

  Makros::unified_assign(Base::data_, toCopy.direct_access(), Base::size_);
  //for (uint k=0; k < toCopy.size(); k++)
  //  data_[k] = toCopy[k];
}

//move constructor
template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D(FlexibleStorage1D<T,ST>&& toTake) : StorageBase<T,ST>(toTake)
{
  reserved_size_ = toTake.reserved_size_;
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::swap(FlexibleStorage1D<T,ST>& toSwap)
{
  std::swap(Base::data_,toSwap.data_);
  std::swap(Base::size_,toSwap.size_);
  std::swap(reserved_size_,toSwap.reserved_size_);
}

template<typename T, typename ST>
FlexibleStorage1D<T,ST>& FlexibleStorage1D<T,ST>::operator=(const FlexibleStorage1D<T,ST>& toCopy)
{
  uint new_res = toCopy.reserved_size();
  if (new_res != reserved_size_) {
    reserved_size_ = new_res;

    if (Base::data_ != 0)
      delete[] Base::data_;
    Base::data_ = new T[reserved_size_];
  }

  Base::size_ = toCopy.size();

  Makros::unified_assign(Base::data_, toCopy.direct_access(), Base::size_);
  //for (uint k=0; k < Base::size_; k++)
  //  Base::data_[k] = toCopy[k];

  return *this;
}

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::~FlexibleStorage1D()
{
}

template<typename T, typename ST>
inline ST FlexibleStorage1D<T,ST>::reserved_size() const
{
  return reserved_size_;
}

template<typename T, typename ST>
inline ST FlexibleStorage1D<T,ST>::capacity() const
{
  return reserved_size_;
}

template<typename T, typename ST>
inline T FlexibleStorage1D<T,ST>::back() const
{
  assert(Base::size_ > 0);
  assert(Base::data_ != 0);
  return Base::data_[Base::size_-1];
}

template<typename T, typename ST>
inline T& FlexibleStorage1D<T,ST>::back()
{
  assert(Base::size_ > 0);
  assert(Base::data_ != 0);
  return Base::data_[Base::size_-1];
}

template<typename T, typename ST>
ST FlexibleStorage1D<T,ST>::append(FlexibleStorage1D<T,ST>::PassType val)
{
  if (Base::size_ == reserved_size_) {

    reserved_size_ = size_t(1.2 * reserved_size_) + 4;

    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, Base::size_);
    //for (uint k=0; k < Base::size_; k++)
    //  new_data[k] = Base::data_[k];

    delete[] Base::data_;
    Base::data_ = new_data;
  }

  const ST k = Base::size_;
  Base::data_[k] = val;

  Base::size_++;

  return k;
}

template<typename T, typename ST>
ST FlexibleStorage1D<T,ST>::move_append(T&& val)
{
  if (Base::size_ == reserved_size_) {

    reserved_size_ = size_t(1.2 * reserved_size_) + 4;

    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, Base::size_);
    //for (uint k=0; k < Base::size_; k++)
    //  new_data[k] = Base::data_[k];

    delete[] Base::data_;
    Base::data_ = new_data;
  }

  const ST k = Base::size_;
  Base::data_[k] = val;

  Base::size_++;

  return k;  
}

//shortcut when you are sure the allocated memory suffices
template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::append_trusting(FlexibleStorage1D<T,ST>::PassType val)
{
  assert(Base::size_ < reserved_size_);
  Base::data_[Base::size_] = val;
  Base::size_++;
}

template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::push_back(FlexibleStorage1D<T,ST>::PassType val)
{
  append(val);
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::append(Storage1D<T,ST>& toAppend)
{
  if (reserved_size_ < Base::size_ + toAppend.size()) {

    reserved_size_ = Base::size_ + toAppend.size() + 2;

    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, Base::size_);
    //for (uint k=0; k < Base::size_; k++)
    //  new_data[k] = Base::data_[k];

    delete[] Base::data_;
    Base::data_ = new_data;
  }

  for (ST k=0; k < toAppend.size(); k++) {
    Base::data_[Base::size_] = toAppend[k];
    Base::size_++;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::append(FlexibleStorage1D<T,ST>& toAppend)
{
  if (reserved_size_ < Base::size_ + toAppend.size()) {

    reserved_size_ = Base::size_ + toAppend.size() + 2;

    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, Base::size_);
    //for (uint k=0; k < Base::size_; k++)
    //  new_data[k] = Base::data_[k];

    delete[] Base::data_;
    Base::data_ = new_data;
  }

  Makros::unified_assign(Base::data_, toAppend.direct_access(), toAppend.size());
  Base::size_ += toAppend.size();  
  //for (ST k=0; k < toAppend.size(); k++) {
  //  Base::data_[Base::size_] = toAppend[k];
  //  Base::size_++;
  //}
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::resize(ST size, bool exact_fit)
{
  if (size > reserved_size_ || size < (reserved_size_ / 3) ) {

    reserved_size_ = size;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, std::min(Base::size_,size));
    //for (uint k=0; k < std::min(size_,size); k++)
    //  new_data[k] = data_[k];

    delete[] Base::data_;
    Base::data_ = new_data;
  }

  Base::size_ = size;

  if (exact_fit && Base::size_ != reserved_size_) {

    reserved_size_ = Base::size_;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, Base::size_);
    //for (uint k=0; k < Base::size_; k++)
    //  new_data[k] = Base::data_[k];

    delete[] Base::data_;
    Base::data_ = new_data;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::fit_exactly()
{
  if (reserved_size_ != Base::size_) {
    
    T* new_data = new T[Base::size_];
    Makros::unified_assign(new_data, Base::data_, Base::size_);

    delete[] Base::data_;
    Base::data_ = new_data;
    reserved_size_ = Base::size_;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::shrink(ST size)
{
  assert(size <= Base::size_);
  Base::size_ = size;
}

template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::shrink_by(ST reduction)
{
  assert(reduction <= Base::size_);
  Base::size_ -= reduction;
}

template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::grow_by_dirty(ST increase)
{
  const ST size = Base::size_ + increase;
 
  if (size > reserved_size_ || size < (reserved_size_ / 3) ) {

    reserved_size_ = size;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, std::min(Base::size_,size));
    //for (uint k=0; k < std::min(Base::size_,size); k++)
    //  new_data[k] = Base::data_[k];

    delete[] Base::data_;
    Base::data_ = new_data;
  }

  Base::size_ = size; 
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::reserve(ST size)
{
  if (size > Base::size_ && size != reserved_size_) {

    reserved_size_ = size;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, Base::data_, Base::size_);

    delete[] Base::data_;
    Base::data_ = new_data;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::clear()
{
  Base::size_ = 0;
}

template<typename T, typename ST>
inline T& FlexibleStorage1D<T,ST>::operator[](ST i) const
{
#ifdef SAFE_MODE
  if (i >= Base::size_) {
    INTERNAL_ERROR << "    invalid access on element " << i
                   << " for FlexibleStorage1D " <<  "\"" << this->name() << "\" of type "
                   //<< Makros::Typename<T>()
                   << typeid(T).name()
                   << " with " << Base::size_ << " (valid) elements. exiting." << std::endl;
    print_trace();
    exit(1);
  }
#endif
  return Base::data_[i];
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::erase(ST pos)
{
  if (pos < Base::size_)
  {
    Routines::downshift_array(Base::data_, pos, Base::size_, 1);
    shrink_by(1);
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::erase_several(ST pos, ST nToErase)
{
  if (pos < Base::size_) 
  {
    nToErase = std::min(nToErase, Base::size_ - pos);
    Routines::downshift_array(Base::data_, pos, Base::size_, nToErase);
    shrink_by(nToErase);    
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::insert(ST pos, T val)
{
  if (pos <= Base::size_)
  {
    grow_by_dirty(1);
    Routines::upshift_array(Base::data_, pos, Base::size_, 1);
    Base::data_[pos] = val;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::insert_several(ST pos, T* vals, ST nData)
{
  if (pos <= Base::size_)
  {
    grow_by_dirty(nData);
    Routines::upshift_array(Base::data_, pos, Base::size_+nData-1, nData);
    Makros::unified_assign(Base::data_+pos, vals, nData);
  }
}

template<typename T, typename ST>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T,ST>& v)
{
  s << "[ ";
  for (int i=0; i < ((int) v.size()) - 1; i++)
    s << v[i] << ",";
  if (v.size() > 0)
    s << v[v.size()-1];
  s << " ]";

  return s;
}

template<typename T, typename ST>
bool operator==(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2)
{
  if (v1.size() != v2.size())
    return false;

  for (ST k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return false;
  }
  return true;
}

template<typename T, typename ST>
bool operator!=(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2)
{
  if (v1.size() != v2.size())
    return true;

  for (ST k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return true;
  }
  return false;
}

/***********************************/

template<typename T, typename ST> 
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D() : name_("unfs1d") {}

template<typename T, typename ST> 
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const std::string& name) : name_(name)
{
}

template<typename T, typename ST> 
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(ST reserved_size, const std::string& name) :
  FlexibleStorage1D<T,ST>(reserved_size), name_(name) {}

//Note: the name is NOT copied
template<typename T, typename ST> 
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T,ST>& toCopy) :
  FlexibleStorage1D<T,ST>(toCopy), name_("unfs1d")
{
}

template<typename T, typename ST> 
NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy) :
  FlexibleStorage1D<T,ST>(toCopy), name_("unfs1d")
{
}

template<typename T, typename ST>
/*virtual*/ const std::string& NamedFlexibleStorage1D<T,ST>::name() const
{
  return name_;
}

template<typename T, typename ST>
void NamedFlexibleStorage1D<T,ST>::operator=(const NamedFlexibleStorage1D<T,ST>& toCopy)
{
  FlexibleStorage1D<T,ST>::operator=(toCopy);
}

template<typename T, typename ST>
void NamedFlexibleStorage1D<T,ST>::operator=(const FlexibleStorage1D<T,ST>& toCopy)
{
  FlexibleStorage1D<T,ST>::operator=(toCopy);
}



#endif