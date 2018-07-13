/*
 * SCC_vector.h
 *
 *  Created on: Jun 26, 2015
 *      Author: Chris Anderson
 *      Modified by: Zach Boyd
 *
 *
 *  A minimal 1d double vector class with move semantic implementation.
 *
 * Decisions:
 * Not making the class backward compatible with C++ versions
 *
 * Use of move semantics
 * Use of nullptr instead of 0 for null pointers (or NULL)
 * Use of std::copy to copy data instead of low level loop
 * Bounds checking completely turned off when _DEBUG is not set
 * Use of assert to facilitate index bounds error location
 *
 * Revised: Nov. 26, 2015
 *        : Jan. 3,  2016   If MS compiler, use std::memcpy instead of std:copy
 *                          to avoid MS compiler warnings about unchecked index ranges.
 *        : Jan. 13,  2016  Added iostream output to facilitate debugging
 */
/*
#############################################################################
#
# Copyright 2015-16-16-16-16-16 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/


#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <new>
using namespace std;

#ifdef  _DEBUG
#include <cstdio>
#else
#define _NDEBUG
#endif
#include <cassert>

#undef _VERBOSE_OPS_

#ifndef _SCC_vector_
#define _SCC_vector_

namespace SCC
{
	class vector
	{
		public:

			vector()
			{
				dataPtr    = nullptr;
				length = 0;
			}

			vector(long n)
			{
				if(n > 0)
				{
					dataPtr    = new double[n];
					length = n;
				}
				else {dataPtr    = nullptr; length = 0;}
			}

			vector(long n, double val)
			{
				if(n > 0)
				{
					dataPtr    = new double[n];
					length = n;

					if(dataPtr==0)
					{
						cout << "Failed allocation in SCC_vector.h" << endl;
						length = 0;
					}

					for(int i=0; i<n; i++)
					{
						dataPtr[i] = val;
					}
				}
				else {dataPtr    = nullptr; length = 0;}
			}

			vector(const vector& V)
			{
#ifdef _VERBOSE_OPS_
				cout << "Standard Copy " << endl;
#endif

				if(V.dataPtr == nullptr)
				{dataPtr = nullptr; length = 0; return;}

				dataPtr     = new double[V.length];
				length = V.length;

#ifdef _MSC_VER
				std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*length);
#else
				std::copy(V.dataPtr, V.dataPtr + length, dataPtr);
#endif

			}

			vector(vector&& V)
			{
#ifdef _VERBOSE_OPS_
				cout << "Move Copy " << endl;
#endif

				dataPtr      = V.dataPtr;
				length   = V.length;
				V.dataPtr    = nullptr;
				V.length = 0;;
			}

			virtual ~vector()
			{
				if(dataPtr != nullptr) delete [] dataPtr;
			}

			void initialize()
			{
				if(dataPtr != nullptr) delete [] dataPtr;
				dataPtr    = nullptr;
				length = 0;
			}

			void initialize(long n)
			{
				if(length != n)
				{
					if(dataPtr != nullptr) delete [] dataPtr;
					if(n > 0)
					{
						dataPtr    = new double[n];
						length = n;
					}
					else {dataPtr    = nullptr; length = 0;}
				}
			}

			void initialize(const vector& V)
			{
				if(V.dataPtr == nullptr)
				{
					if(dataPtr != nullptr) delete [] dataPtr;
					dataPtr = nullptr;
					length = 0;
					return;
				}

				if(length != V.length)
				{
					if(dataPtr != nullptr) delete [] dataPtr;
					dataPtr     = new double[V.length];
					length = V.length;
				}

#ifdef _MSC_VER
				std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*length);
#else
				std::copy(V.dataPtr, V.dataPtr + length, dataPtr);
#endif

			}

			void initialize(vector&& V)
			{
				if(V.dataPtr == nullptr)
				{
					if(dataPtr != nullptr) delete [] dataPtr;
					dataPtr = nullptr;
					length = 0;
					return;
				}

				if(dataPtr != nullptr) delete [] dataPtr;

				dataPtr      = V.dataPtr;
				length   = V.length;
				V.dataPtr    = 0;
				V.length = 0;
			}

			// Assignment operators : Being careful with nullptr instances

			vector& operator=(const vector& V)
			{
#ifdef _VERBOSE_OPS_
				cout << "Standard Assignment" << endl;
#endif

				if (this != &V)
				{
					if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
					{
						length  = V.length;
						dataPtr     = new double[length];
#ifdef _MSC_VER
						std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*length);
#else
						std::copy(V.dataPtr, V.dataPtr + length, dataPtr);
#endif
					}
					else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; length = 0; return *this;}
					else
					{
						assert(sizeCheck(this->length,V.length));
#ifdef _MSC_VER
						std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*length);
#else
						std::copy(V.dataPtr, V.dataPtr + length, dataPtr);
#endif
					}
				}
				return *this;
			}

			vector& operator=(vector&& V)
			{
#ifdef _VERBOSE_OPS_
				cout << "Move Assignment" << endl;
#endif

				if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
				{
					dataPtr      = V.dataPtr;
					length   = V.length;
					V.dataPtr    = nullptr;
					V.length = 0;
				}
				else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; length = 0; return *this;}
				else
				{
					assert(sizeCheck(this->length,V.length));

					// Remove existing data

					delete [] dataPtr;

					dataPtr      = V.dataPtr;
					length   = V.length;
					V.dataPtr    = nullptr;
					V.length = 0;
				}
				return *this;
			}

			inline void transformValues(std::function<double(double)> F)
			{
				for(long k =0; k < length; k++)
				{dataPtr[k] = F(dataPtr[k]);}
			}

			vector applyFunction(std::function<double(double)> F)
			{
#ifdef _VERBOSE_OPS_
				cout  << "F(*this)" << endl;
#endif

				vector R(*this);
				R.transformValues(F);
				return std::move(R);
			}


			inline void operator+=(const  vector& D)
			{
				assert(sizeCheck(this->length,D.length));
				for(long i = 0; i < this->length; i++)
				{
					dataPtr[i] += D.dataPtr[i];
				}
			}


			friend vector operator+(const vector& A, const vector& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A + &B" << endl;
#endif

				assert(A.sizeCheck(A.length,B.length));
				vector R(A);
				R += B;
				return std::move(R);
			}

			friend vector operator+(const vector& A, vector&& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A + &&B " << endl;
#endif

				assert(A.sizeCheck(A.length,B.length));
				B += A;
				return std::move(B);
			}

			friend vector operator+(vector&& A, const vector& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A +  &B" << endl;
#endif

				assert(B.sizeCheck(A.length,B.length));
				A += B;
				return std::move(A);
			}

			friend vector operator+(vector&& A, vector&& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A +  &&B" << endl;
#endif

				assert(B.sizeCheck(A.length,B.length));
				A += B;
				return std::move(A);
			}

			inline void operator-=(const  vector& D)
			{
				assert(sizeCheck(this->length,D.length));
				for(long i = 0; i < this->length; i++)
				{
					dataPtr[i] -= D.dataPtr[i];
				}
			}

			friend vector operator-(const vector& A, const vector& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A - &B" << endl;
#endif

				assert(A.sizeCheck(A.length,B.length));
				vector R(A);
				R -= B;
				return std::move(R);
			}

			friend vector operator-(const vector& A, vector&& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A - &&B " << endl;
#endif

				assert(A.sizeCheck(A.length,B.length));
				for(long i = 0; i < B.length; i++)
				{
					B.dataPtr[i] = A.dataPtr[i] - B.dataPtr[i];
				}
				return std::move(B);
			}

			friend vector operator-(vector&& A, const vector& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A -  &B" << endl;
#endif

				assert(B.sizeCheck(A.length,B.length));
				A -= B;
				return std::move(A);
			}

			friend vector operator-(vector&& A, vector&& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A -  &&B" << endl;
#endif

				assert(B.sizeCheck(A.length,B.length));
				A -= B;
				return std::move(A);
			}

			friend vector operator-(const vector& A)
			{
#ifdef _VERBOSE_OPS_
				cout  << "-&A" << endl;
#endif

				vector R(A);
				R *= -1.0;
				return std::move(R);
			}

			friend vector operator-(vector&& A)
			{
#ifdef _VERBOSE_OPS_
				cout  << "- &&A" << endl;
#endif

				A *= -1.0;
				return std::move(A);
			}

			friend vector operator+(const vector& A)
			{
#ifdef _VERBOSE_OPS_
				cout  << "+&A" << endl;
#endif

				vector R(A);
				return std::move(R);
			}

			friend vector operator+(vector&& A)
			{
#ifdef _VERBOSE_OPS_
				cout  << "+&&A" << endl;
#endif

				return std::move(A);
			}

			inline void operator*=(const vector& B)
			{
				assert(sizeCheck(length,B.length));
				for(long i = 0; i < this->length; i++)
				{
					dataPtr[i] *= B.dataPtr[i];
				}
			}

			friend vector operator*(const vector& A, const vector& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A * &B" << endl;
#endif

				assert(A.sizeCheck(A.length,B.length));
				vector R(A);
				R *= B;
				return std::move(R);
			}

			friend vector operator*(const vector& A, vector&& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A * &&B " << endl;
#endif

				assert(A.sizeCheck(A.length,B.length));
				B *= A;
				return std::move(B);
			}

			friend vector operator*(vector&& A, const vector& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A *  &B" << endl;
#endif

				assert(B.sizeCheck(A.length,B.length));
				A *= B;
				return std::move(A);
			}

			friend vector operator*(vector&& A, vector&& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A *  &&B" << endl;
#endif

				assert(B.sizeCheck(A.length,B.length));
				A *= B;
				return std::move(A);
			}

			inline void operator*=(const double alpha)
			{
				for(long i = 0; i < this->length; i++)
				{
					dataPtr[i] *= alpha;
				}
			}

			friend vector operator*(const double alpha, const vector& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "alpha* &B" << endl;
#endif

				vector R(B);
				R *= alpha;
				return std::move(R);
			}

			friend vector operator*(const double alpha, vector&& B)
			{
#ifdef _VERBOSE_OPS_
				cout  << "alpha*+ &&B " << endl;
#endif

				B *= alpha;
				return std::move(B);
			}

			friend vector operator*(vector& A, const double alpha)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A*alpha" << endl;
#endif

				vector R(A);
				R *= alpha;
				return std::move(R);
			}

			friend vector operator*(vector&& A,const double alpha)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A * alpha" << endl;
#endif

				A *= alpha;
				return std::move(A);
			}


			inline void operator/=(const vector& B)
			{
				assert(sizeCheck(length,B.length));
				for(long i = 0; i < this->length; i++)
				{
					dataPtr[i] /= B.dataPtr[i];
				}
			}

			inline void operator/=(const double alpha)
			{
				for(long i = 0; i < this->length; i++)
				{
					dataPtr[i] /= alpha;
				}
			}

			friend vector operator/(vector& A, const double alpha)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&A/alpha" << endl;
#endif

				vector R(A);
				R /= alpha;
				return std::move(R);
			}

			friend vector operator/(vector&& A,const double alpha)
			{
#ifdef _VERBOSE_OPS_
				cout  << "&&A/alpha" << endl;
#endif

				A /= alpha;
				return std::move(A);
			}
			/*!  Sets all values of the vector to d. */

			void setToValue(double d)
			{
				for(long i = 0; i < length; i++)
				{dataPtr[i] =  d;}
			}

			void addValue(double d)
			{
				for(long i = 0; i < length; i++)
				{dataPtr[i] += d;}
			}


			/*!  Standard vector dot product.  */

			virtual double dot(const vector& v) const
			{
				double dotVal = 0.0;
				for(long i = 0; i < length; i++)
				{
					dotVal += (dataPtr[i]*v.dataPtr[i]);
				}
				return dotVal;
			}

			/*!  Maximal absolute value of the elements of the vector. */

			virtual double normInf() const
			{
				double valMax = 0.0;
				for(long i = 0; i < length; i++)
				{
					valMax = (valMax > abs(dataPtr[i])) ? valMax : abs(dataPtr[i]);
				}
				return valMax;
			}

			/*!  The Euclidean norm of the vector. */

			virtual double norm2() const
			{
				double val = 0.0;
				for(long i = 0; i < length; i++)
				{
					val += (dataPtr[i]*dataPtr[i]);
				}
				return sqrt(abs(val));
			}

			// Selected BLAS interface

			/*! BLAS Euclidean norm of the vector */

			virtual double nrm2()
			{
				return norm2();
			}

			/*! BLAS scalar multiplication  */

			virtual void scal(double alpha)
			{
				(*this) *= alpha;
			}

			/*! BLAS copy : this  <-v   */

			void copy(const vector& v)
			{
				this->operator=(v);
			}

			void copy(vector&& v)
			{
				this->operator=((vector&&)v);
			}


			/*! BLAS axpby : this  <- alpha*v + beta*this  */

			virtual void axpby(double alpha, const vector& v, double beta)
			{
				assert(sizeCheck(this->length,v.length));
				for(long i = 0; i < length; i++)
				{dataPtr[i]  = alpha*v.dataPtr[i] + beta*dataPtr[i];}
			}

			/*! BLAS axpy : this  <- alpha*v + this  */

			virtual void axpy(double alpha, const vector& v)
			{
				assert(sizeCheck(this->length,v.length));
				for(long i = 0; i < length; i++)
				{dataPtr[i]  += alpha*v.dataPtr[i];}
			}

			/*!  Returns the dimension of the vector */

			virtual long getDimension()
			{
				return length;
			}

#ifndef _NDEBUG
			double&  operator()(long i1)
			{
				assert(boundsCheck(i1, 0, length-1,1));
				return *(dataPtr +  i1);
			};

			const double&  operator()(long i1) const
			{
				assert(boundsCheck(i1, 0, length-1,1));
				return *(dataPtr +  i1);
			};

#else
			inline double&  operator()(long i1)
			{
				return *(dataPtr + i1);
			};

			inline const double&  operator()(long i1) const
			{
				return *(dataPtr + i1);
			};
#endif


			/*!  Outputs vector values to a stream. */

			friend ostream& operator<<(ostream& outStream, const vector& V)
			{

				long i;
				for(i = 0; i <  V.length; i++)
				{
					outStream <<  setprecision(3) <<  std::right << setw(10) << V(i) << " ";
					outStream << endl;
				}
				return outStream;
			}

			long getSize()  const {return length;}

			long getIndex1Size()  const {return length;}

			double* getDataPointer(){return dataPtr;}
			const  double* getDataPointer()  const  {return dataPtr;}

			double*       dataPtr;
			long length;




			//###################################################################
			//                      Bounds Checking
			//###################################################################
			//
#ifdef _DEBUG
			bool boundsCheck(long i, long begin, long end, int coordinate) const
			{
				if((i < begin)||(i  > end))
				{
					cerr << "SCC::vector index " << coordinate << " out of bounds " << endl;
					cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << endl;
					return false;
				}
				return true;
			}
#else
			bool boundsCheck(long, long, long, int) const {return true;}
#endif

#ifdef _DEBUG
			bool sizeCheck(long size1, long size2)
			{
				if(size1 != size2)
				{
					cerr << "SCC::vector sizes are incompatible : " << size1 << " != " << size2;
					return false;
				}
				return true;
			}

			bool sizeCheck(long size1, long size2) const
			{
				if(size1 != size2)
				{
					cerr << "SCC::vector sizes are incompatible : " << size1 << " != " << size2;
					return false;
				}
				return true;
			}
#else
			bool sizeCheck(long, long) {return true;}
			bool sizeCheck(long, long) const{return true;}
#endif
			void print()
			{
				for(long i=0; i<length; i++)
				{
					cout << dataPtr[i] << endl;
				}
			}

			void allocation_check()
			{
			}
	};
}

#endif /* SCC_vector_ */
