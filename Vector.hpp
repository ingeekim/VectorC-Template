///////////////////////////////////////////////////////////////////////////////
// Filename: Vector.hpp
// Mathematical Vector class definition.
//
// The class modifies that the array index starts from 1
// Performs the basic operations:
// 1. equality: == : returns true/false
// 2. inequality: != : returns true/false
// 3. assignment: = : A = B : assign B contents to A
// 4. component: [i] : access to the i the element either as lvalue or rvalue
// 5. complex conjugate: returns the complex conjugate of the elements,
//      returns a Vector
// 6. scalar product: A*B : performs the scalar product of A and B,
//	returns a double
// 7. scalar multiplication: A*a : 
//      performs the scalar multiplication with a scalar (double) a, 
//	returns a Vector
// 8. vector product: A^B : performs the vector product of A and B,
//      returns a Vector
// 9. vector addition: A+B : performs the vector addition of A and B, 
//      returns a Vector
// 10. vector subtraction: A-B : performs the vector subtraction of A and B,
//      returns a Vector
//
// In Gee Kim, July, 2011 
///////////////////////////////////////////////////////////////////////////////

#ifndef VECTOR_HPP
#define VECTOR_HPP

// system headers
#include <iostream>
using std::cerr;
using std::endl;

#include <cmath>
using std::sqrt;

#include <complex>
using std::conj;

#include <string>
using std::string;

// user-defined headers
#include "Array.hpp"

// user-defined libraries
void quit(string &);

template<typename T>
class Vector : public Array<T>
{ 
public:
  // constructors
  Vector() : Array<T>(){};
  Vector(long i) : Array<T>(i) {};
  Vector(long i, const T &a) : Array<T>(i, a){};
  Vector(long i, const T *a) : Array<T>(i, a){};
  Vector(const Vector<T> &a);

  // destructors
  ~Vector(){};

  // get norm
  double getNorm() const;

  // assignment operators
  Vector<T> &operator=(const Vector<T> &);
  Vector<T> &operator=(const T &);

  // equality operator
  bool operator==(const Vector &) const;

  // inequality operator
  bool operator!=(const Vector &right) const
  {
    return !(*this == right); // invokes Vector::operator== 
  }
  
  // scalar product
  double operator*(const Vector &) const;

  // complex conjugate
  Vector<T> conjugate(const Vector &) const;

  // scalar multiplication
  Vector<T> operator*(const T) const;

  // scalar division
  Vector<T> operator/(const T) const;

  // vector product
  Vector<T> operator^(const Vector &) const;

  // vector addition
  Vector<T> operator+(const Vector &) const;

  // vector subtraction
  Vector<T> operator-(const Vector &) const;

}; // end class Vector

// copy constructor
template<typename T>
Vector<T>::Vector(const Vector<T> &rhs)
{
  iSize = rhs.iSize;
  iOffset = rhs.iOffset;
  ptrArray = new T[iSize];
  for(int i = 0; i < iSize;  i++)
    ptrArray[i] = rhs.ptrArray[i];
} // end of copy constructor

// get norm
template<typename T>
double Vector<T>::getNorm() const
{
  double square = 0.0;
  
  for(long i = 0; i < iSize; i++) {
    square += static_cast<double> (ptrArray[i] * ptrArray[i]);
  }

  return sqrt(square);
}

// operator overloading the lvalue assignment
template<typename T>
Vector<T> &Vector<T>::operator=(const Vector<T> &rhs)
{
  if(this != &rhs)
  {
    if(ptrArray != 0) delete [] (ptrArray);
    iSize = rhs.iSize;
    iOffset = rhs.iOffset;
    ptrArray = new T[iSize];
  }
  for(long i = 0; i < iSize; i++)
    ptrArray[i] = rhs.ptrArray[i];

  return *this;
} // end assignment operator

// operator overloading the rvalue assignment
template<typename T>
Vector<T> &Vector<T>::operator=(const T &a)
{
  for(long i = 0; i < iSize; i++)
    ptrArray[i] = a;

  return *this;
} // end assignment operator

// equality operator
template<typename T>
bool Vector<T>::operator==(const Vector &right) const
{
  if(iSize != right.iSize)
    return false;

  for(long i = 0; i < iSize; i++)
    if(ptrArray[i] != right.ptrArray[i])
      return false; 
  
  return true;
} // end function bool Vector::operator==

// complex conjugate
template<typename T>
Vector<T> Vector<T>::conjugate(const Vector &vec) const
{
  iSize = vec.size();
  T *vecPtr = new T[iSize];
  for(long i = 0; i < iSize; i++)
    vecPtr[i] = conj(vec.ptrArray[i]);

  return Vector<T>(iSize, vecPtr);
}

// scalar product
template<typename T>
double Vector<T>::operator*(const Vector &right) const
{
  double scalar = 0.0;
  long rightDim = right.size();
  long minDim = ( iSize <= rightDim ? iSize : rightDim);

  // implementation subscript should starts from 0
  for(long i = 0; i < minDim; i++)
      scalar += ptrArray[i] * right.ptrArray[i];

  return scalar;
} // end function double Vector::operator*

// scalar multiplication
template<typename T>
Vector<T> Vector<T>::operator*(const T alpha) const
{
  T *vecPtr = new double[iSize]; 

  for(long i = 0; i < iSize; i++)
    vecPtr[i] = ptrArray[i] * alpha;
  
  return Vector<T>(iSize, vecPtr);
} // end function Vector Vector::operator*

template<typename T>
Vector<T> Vector<T>::operator/(const T alpha) const
{
  if( alpha == 0.0 ) {
    string message = "Error: Division by zero!";
    quit(message);
  }

  T *vecPtr = new T[iSize]; 

  for(long i = 0; i < iSize; i++)
    vecPtr[i] = ptrArray[i] / alpha;

  return Vector<T>(iSize, vecPtr);
} // end function Vector Vector::operator/

template<typename T>
Vector<T> Vector<T>::operator^(const Vector<T> &right) const
{
  // take the largest dimension
  long rDim = right.size();
  long maxDim = ( iSize >= rDim ? iSize : rDim);
  if(maxDim != 3) {
    string message =
      "Error: vector product is implemented for 3d vectors only";
    quit(message);
  }

  // build the temporary vectors
  T *vecPtr = new T[maxDim];

  // initializes the temporary vectors
  for(long i = 0; i < maxDim; i++)
    vecPtr[i] = 0.0;

  // perform cross product
  for(long i = 1; i <= 3; i++) {
    long j = i + 1;
    long k = j + 1;
    if( j > 3 ) j %= 3;
    if( k > 3 ) k %= 3;
    vecPtr[i-1] += ptrArray[j-1] * right.ptrArray[k-1];
    vecPtr[i-1] -= ptrArray[k-1] * right.ptrArray[j-1];
  }
 
  return Vector<T>(maxDim, vecPtr);
} // end function Vector Vector::operator^

template<typename T>
Vector<T> Vector<T>::operator+(const Vector &right) const
{
  // take the largest dimension
  long rightDim = right.size();
  long maxDim = ( iSize >= rightDim ? iSize : rightDim);

  // build the temporary vectors
  double *vecPtr = new double[maxDim];

  // initializes the temporary vectors
  for(long i = 0; i < maxDim; i++)
    vecPtr[i] = 0.0;

  // implementation subscript should starts from 0
  for(long i = 0; i < maxDim; i++)
    vecPtr[i] = ptrArray[i] + right.ptrArray[i];

  return Vector<T>(maxDim, vecPtr);
} // end function Vector Vector::operator+

template<typename T>
Vector<T> Vector<T>::operator-(const Vector &right) const
{
  // take the largest dimension
  long rightDim = right.size();
  long maxDim = ( iSize >= rightDim ? iSize : rightDim);

  // build the temporary vectors
  T *vecPtr = new T[maxDim];
  T *lPtr = new T[maxDim];
  T *rPtr = new T[maxDim];

  // initializes the temporary vectors
  for(long i = 0; i < maxDim; i++)
    vecPtr[i] = lPtr[i] = rPtr[i] = 0.0;

  // implementation subscript should starts from 0
  for(long i = 0; i < maxDim; i++)
    vecPtr[i] = ptrArray[i] - right.ptrArray[i];

  return Vector<T>(maxDim, vecPtr);
} // end function Vector Vector::operator-

#endif
