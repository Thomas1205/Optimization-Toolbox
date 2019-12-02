/*-*-c++-*-*/
/*** first version written by Thomas Schoenemann as a private person without employment, September 2009 ***/
/*** much refined by Thomas Schoenemann  at Lund University, Sweden, the University of Pisa, Italy, ***
 *** and the University of Düsseldorf, Germany 2010 - 2012 **/
/*** if you desire the checked version, make sure your compiler defines the option SAFE_MODE on the command line ***/

#ifndef STORAGE_3D_HH
#define STORAGE_3D_HH

#include "storage1D.hh"

template<typename ST>
struct Dim3D {
  
  Dim3D(ST xDim, ST yDim, ST zDim) : xDim_(xDim), yDim_(yDim), zDim_(zDim) {}
  
  ST xDim_;
  ST yDim_;
  ST zDim_;
};

template<typename T, typename ST=size_t>
class Storage3D : public StorageBase<T,ST> {
public:

  Storage3D();

  //copy constructor
  Storage3D(const Storage3D& toCopy);

  Storage3D(ST xDim, ST yDim, ST zDim);
  
  Storage3D(const Dim3D<ST> dims);

  Storage3D(ST xDim, ST yDim, ST zDim, const T default_value);

  Storage3D(const Dim3D<ST> dims, const T default_value);

  ~Storage3D();

  inline const T& operator()(ST x, ST y, ST z) const;

  inline T& operator()(ST x, ST y, ST z);

  virtual const std::string& name() const;

  inline ST xDim() const;

  inline ST yDim() const;

  inline ST zDim() const;

  Dim3D<ST> dims() const;

  void set_x(ST y, ST z, const Storage1D<T,ST>& vec);

  void set_y(ST x, ST z, const Storage1D<T,ST>& vec);

  void set_z(ST x, ST y, const Storage1D<T,ST>& vec);

  void operator=(const Storage3D<T,ST>& toCopy);

#ifdef SAFE_MODE
  //for some reason g++ allows to assign an object of type T, but this does NOT produce the effect one would expect
  // => define this operator in safe mode, only to check that such an assignment is not made
  void operator=(const T& invalid_object);
#endif

  //existing positions are copied, new ones are uninitialized
  void resize(ST newxDim, ST newyDim, ST newzDim);

  inline void resize(const Dim3D<ST> dims)
  {
    resize(dims.xDim_, dims.yDim_, dims.zDim_);
  }

  //existing positions are copied, new ones are uninitialized
  void resize(ST newxDim, ST newyDim, ST newzDim, const T default_value);

  inline void resize(const Dim3D<ST> dims, T default_value)
  {
    resize(dims.xDim_, dims.yDim_, dims.zDim_, default_value);
  }

  //all elements are uninitialized after this operation
  void resize_dirty(ST newxDim, ST newyDim, ST newzDim);

  inline void resize_dirty(const Dim3D<ST> dims)
  {
    resize_dirty(dims.xDim_, dims.yDim_, dims.zDim_);
  }

protected:
  ST xDim_;
  ST yDim_;
  ST zDim_;

  static const std::string stor3D_name_;
};


template<typename T, typename ST=size_t>
class NamedStorage3D : public Storage3D<T,ST> {
public:

  NamedStorage3D();

  NamedStorage3D(std::string name);

  NamedStorage3D(ST xDim, ST yDim, ST zDim, std::string name);

  NamedStorage3D(ST xDim, ST yDim, ST zDim, T default_value, std::string name);

  virtual const std::string& name() const;

  inline void operator=(const Storage3D<T,ST>& toCopy);

  //NOTE: the name is NOT copied
  inline void operator=(const NamedStorage3D<T,ST>& toCopy);

protected:
  std::string name_;
};

template<typename T, typename ST>
bool operator==(const Storage3D<T,ST>& v1, const Storage3D<T,ST>& v2);

template<typename T, typename ST>
bool operator!=(const Storage3D<T,ST>& v1, const Storage3D<T,ST>& v2);


namespace Makros {

  template<typename T, typename ST>
  class Typename<Storage3D<T,ST> > {
  public:

    std::string name() const
    {
      return "Storage3D<" + Makros::Typename<T>() + "," + Makros::Typename<ST>() + "> ";
    }
  };

  template<typename T>
  class Typename<Storage3D<T> > {
  public:

    std::string name() const
    {
      return "Storage3D<" + Makros::Typename<T>() + "> ";
    }
  };


  template<typename T, typename ST>
  class Typename<NamedStorage3D<T,ST> > {
  public:

    std::string name() const
    {
      return "NamedStorage3D<" + Makros::Typename<T>() + "," + Makros::Typename<ST>() + "> ";
    }
  };

  template<typename T>
  class Typename<NamedStorage3D<T> > {
  public:

    std::string name() const
    {
      return "NamedStorage3D<" + Makros::Typename<T>() + "> ";
    }
  };

}


/******************************************** implementation **************************************************/
template<typename T, typename ST>
/*static*/ const std::string Storage3D<T,ST>::stor3D_name_ = "unnamed Storage3D";

template<typename T, typename ST> 
Storage3D<T,ST>::Storage3D() : StorageBase<T,ST>(), xDim_(0), yDim_(0), zDim_(0) {}

template<typename T, typename ST> 
Storage3D<T,ST>::Storage3D(const Storage3D<T,ST>& toCopy) : StorageBase<T,ST>(toCopy.xDim()*toCopy.yDim()*toCopy.zDim())
{
  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  zDim_ = toCopy.zDim();

  Makros::unified_assign(StorageBase<T,ST>::data_, toCopy.direct_access(), StorageBase<T,ST>::size_);

  //for (ST i=0; i < size_; i++)
  //  data_[i] = toCopy.direct_access(i);

  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
}

template<typename T, typename ST> 
Storage3D<T,ST>::Storage3D(ST xDim, ST yDim, ST zDim) : StorageBase<T,ST>(xDim*yDim*zDim), xDim_(xDim), yDim_(yDim), zDim_(zDim)
{
}

template<typename T, typename ST> 
Storage3D<T,ST>::Storage3D(const Dim3D<ST> dims) 
 : StorageBase<T,ST>(dims.xDim_*dims.yDim_*dims.zDim_), xDim_(dims.xDim), yDim_(dims.yDim_), zDim_(dims.zDim_) {}

template<typename T, typename ST> 
Storage3D<T,ST>::Storage3D(ST xDim, ST yDim, ST zDim, const T default_value) 
 : StorageBase<T,ST>(xDim*yDim*zDim,default_value), xDim_(xDim), yDim_(yDim), zDim_(zDim)
{
}

template<typename T, typename ST> 
Storage3D<T,ST>::Storage3D(const Dim3D<ST> dims, const T default_value) 
 : StorageBase<T,ST>(dims.xDim_*dims.yDim_*dims.zDim_,default_value), xDim_(dims.xDim), yDim_(dims.yDim_), zDim_(dims.zDim_) {}

template<typename T, typename ST> 
Storage3D<T,ST>::~Storage3D()
{
}

template<typename T, typename ST> 
void Storage3D<T,ST>::set_x(ST y, ST z, const Storage1D<T,ST>& vec)
{
  assert(y < yDim_);
  assert(z < zDim_);
  assert(vec.size() == xDim_);
  
  for (ST x = 0; x < xDim_; x++) 
    (*this)(x, y, z) = vec.direct_access(x);
}

template<typename T, typename ST> 
void Storage3D<T,ST>::set_y(ST x, ST z, const Storage1D<T,ST>& vec) 
{
  assert(x < xDim_);
  assert(z < zDim_);
  assert(vec.size() == yDim_);
  
  for (ST y = 0; y < yDim_; y++) 
    (*this)(x, y, z) = vec.direct_access(y);  
}

template<typename T, typename ST> 
void Storage3D<T,ST>::set_z(ST x, ST y, const Storage1D<T,ST>& vec)
{
  assert(x < xDim_);
  assert(y < yDim_);
  assert(vec.size() == zDim_);
  
  T* data =  StorageBase<T,ST>::data_ + (y*xDim_+x) * zDim_;

  Makros::unified_assign(data, vec.direct_access(), zDim_);
  
  //for (ST z = 0; z < zDim_; z++) 
  //  data[z] = vec.direct_access(z);
}

template<typename T, typename ST>
OPTINLINE const T& Storage3D<T,ST>::operator()(ST x, ST y, ST z) const
{
#ifdef SAFE_MODE
  if (x >= xDim_ || y >= yDim_ || z >= zDim_) {

    INTERNAL_ERROR << "     invalid const access on element (" << x << "," << y << "," << z << ") of 3D-storage \""
                   << this->name() << "\" of type "
                   << Makros::Typename<T>()
                   //<< Makros::get_typename(typeid(T).name())
                   << ":" << std::endl;
    std::cerr << "     dimensions " << xDim_ << "x" << yDim_ << "x" << zDim_ << " exceeded. Exiting..." << std::endl;
    print_trace();
    exit(1);
  }
#endif
  return StorageBase<T,ST>::data_[(y*xDim_+x)*zDim_+z];
}

template<typename T, typename ST>
OPTINLINE T& Storage3D<T,ST>::operator()(ST x, ST y, ST z)
{
#ifdef SAFE_MODE
  if (x >= xDim_ || y >= yDim_ || z >= zDim_) {

    INTERNAL_ERROR << "     invalid access on element (" << x << "," << y << "," << z << ") of 3D-storage \""
                   << this->name() << "\" of type "
                   //<< Makros::Typename<T>()
                   << typeid(T).name()
                   //<< Makros::get_typename(typeid(T).name())
                   << ":" << std::endl;
    std::cerr << "     dimensions " << xDim_ << "x" << yDim_ << "x" << zDim_ << " exceeded. Exiting..." << std::endl;
    print_trace();
    exit(1);
  }
#endif
  return StorageBase<T,ST>::data_[(y*xDim_+x)*zDim_+z];
}

template<typename T, typename ST>
/*virtual*/ const std::string& Storage3D<T,ST>::name() const
{
  return stor3D_name_;
}

template<typename T, typename ST>
inline ST Storage3D<T,ST>::xDim() const
{
  return xDim_;
}

template<typename T, typename ST>
inline ST Storage3D<T,ST>::yDim() const
{
  return yDim_;
}

template<typename T, typename ST>
inline ST Storage3D<T,ST>::zDim() const
{
  return zDim_;
}

template<typename T, typename ST>
Dim3D<ST> Storage3D<T,ST>::dims() const
{
  return Dim3D<ST>(xDim_,yDim_,zDim_);
}

template<typename T, typename ST>
void Storage3D<T,ST>::operator=(const Storage3D<T,ST>& toCopy)
{
  if (StorageBase<T,ST>::size_ != toCopy.size()) {
    if (StorageBase<T,ST>::data_ != 0) {
      delete[] StorageBase<T,ST>::data_;
    }

    StorageBase<T,ST>::size_ = toCopy.size();
    StorageBase<T,ST>::data_ = new T[StorageBase<T,ST>::size_];
  }

  xDim_ = toCopy.xDim();
  yDim_ = toCopy.yDim();
  zDim_ = toCopy.zDim();
  
  const size_t size = StorageBase<T,ST>::size_;
  assert(size == xDim_*yDim_*zDim_);

  Makros::unified_assign(StorageBase<T,ST>::data_, toCopy.direct_access(), size);

  // for (ST i=0; i < size_; i++)
    // data_[i] = toCopy.direct_access(i);

  //this is faster for basic types but it fails for complex types where e.g. arrays have to be copied
  //memcpy(data_,toCopy.direct_access(),size_*sizeof(T));
}

#ifdef SAFE_MODE
//for some reason g++ allows to assign an object of type T, but this does NOT produce the effect one would expect
// => define this operator in safe mode, only to check that such an assignment is not made
template<typename T,typename ST>
void Storage3D<T,ST>::operator=(const T& invalid_object)
{
  INTERNAL_ERROR << "assignment of an atomic entity to Storage1D \"" << this->name() << "\" of type "
                 << Makros::Typename<T>()
                 << " with " << StorageBase<T,ST>::size_ << " elements. exiting." << std::endl;
}
#endif


//existing positions are copied, new ones are uninitialized
template<typename T, typename ST>
void Storage3D<T,ST>::resize(ST newxDim, ST newyDim, ST newzDim)
{
  ST new_size = newxDim*newyDim*newzDim;

  if (newxDim != xDim_ || newyDim != yDim_ || newzDim != zDim_) {
    T* new_data = new T[new_size];

    if (StorageBase<T,ST>::data_ != 0) {

      //copy existing elements
      for (ST x=0; x < std::min(xDim_,newxDim); x++) {
        for (ST y=0; y < std::min(yDim_,newyDim); y++) {
          for (ST z=0; z < std::min(zDim_,newzDim); z++) {
            new_data[(y*newxDim+x)*newzDim+z] = StorageBase<T,ST>::data_[(y*xDim_+x)*zDim_+z];
          }
        }
      }

      delete[] StorageBase<T,ST>::data_;
    }
    StorageBase<T,ST>::data_ = new_data;
    StorageBase<T,ST>::size_ = new_size;
    xDim_ = newxDim;
    yDim_ = newyDim;
    zDim_ = newzDim;
  }
}

//existing positions are copied, new ones are uninitialized
template<typename T, typename ST>
void Storage3D<T,ST>::resize(ST newxDim, ST newyDim, ST newzDim, const T default_value)
{
  const ST new_size = newxDim*newyDim*newzDim;

  if (newxDim != xDim_ || newyDim != yDim_ || newzDim != zDim_) {

    T* new_data = new T[new_size];
    std::fill_n(new_data,new_size,default_value);
    //for (ST i=0; i < new_size; i++)
    //  new_data[i] = default_value;

    if (StorageBase<T,ST>::data_ != 0) {

      //copy existing elements
      for (ST x=0; x < std::min(xDim_,newxDim); x++) {
        for (ST y=0; y < std::min(yDim_,newyDim); y++) {
          for (ST z=0; z < std::min(zDim_,newzDim); z++) {
            new_data[(y*newxDim+x)*newzDim+z] = StorageBase<T,ST>::data_[(y*xDim_+x)*zDim_+z];
          }
        }
      }

      delete[] StorageBase<T,ST>::data_;
    }
    StorageBase<T,ST>::data_ = new_data;
    StorageBase<T,ST>::size_ = new_size;
    xDim_ = newxDim;
    yDim_ = newyDim;
    zDim_ = newzDim;
  }
}

//all elements are uninitialized after this operation
template<typename T, typename ST>
void Storage3D<T,ST>::resize_dirty(ST newxDim, ST newyDim, ST newzDim)
{
  if (newxDim != xDim_ || newyDim != yDim_ || newzDim != zDim_) {

    if (StorageBase<T,ST>::data_ != 0)
      delete[] StorageBase<T,ST>::data_;

    xDim_ = newxDim;
    yDim_ = newyDim;
    zDim_ = newzDim;
    StorageBase<T,ST>::size_ = xDim_*yDim_*zDim_;

    StorageBase<T,ST>::data_ = new T[StorageBase<T,ST>::size_];
  }
}

/***********************/

template<typename T, typename ST> NamedStorage3D<T,ST>::NamedStorage3D() : Storage3D<T,ST>(), name_("yyy") {}

template<typename T, typename ST> NamedStorage3D<T,ST>::NamedStorage3D(std::string name) : Storage3D<T,ST>(), name_(name) {}

template<typename T, typename ST> NamedStorage3D<T,ST>::NamedStorage3D(ST xDim, ST yDim, ST zDim, std::string name) :
  Storage3D<T,ST>(xDim,yDim,zDim), name_(name) {}

template<typename T, typename ST> NamedStorage3D<T,ST>::NamedStorage3D(ST xDim, ST yDim, ST zDim, T default_value, std::string name)
  : Storage3D<T,ST>(xDim,yDim,zDim,default_value), name_(name) {}

template<typename T, typename ST>
/*virtual*/ const std::string& NamedStorage3D<T,ST>::name() const
{
  return name_;
}

template<typename T, typename ST>
inline void NamedStorage3D<T,ST>::operator=(const Storage3D<T,ST>& toCopy)
{
  Storage3D<T,ST>::operator=(toCopy);
}

//NOTE: the name is NOT copied
template<typename T, typename ST>
inline void NamedStorage3D<T,ST>::operator=(const NamedStorage3D<T,ST>& toCopy)
{
  Storage3D<T,ST>::operator=(static_cast<Storage3D<T,ST> >(toCopy));
}


template<typename T, typename ST>
bool operator==(const Storage3D<T,ST>& v1, const Storage3D<T,ST>& v2)
{
  if (v1.xDim() != v2.xDim() || v1.yDim() != v2.yDim() || v1.zDim() != v2.zDim())
    return false;

  for (ST i=0; i < v1.size(); i++) {
    if (v1.direct_access(i) != v2.direct_access(i))
      return false;
  }

  return true;
}

template<typename T, typename ST>
bool operator!=(const Storage3D<T,ST>& v1, const Storage3D<T,ST>& v2)
{
  return !operator==(v1,v2);
}

#endif
