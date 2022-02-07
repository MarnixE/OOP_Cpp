#include <cmath>
#include <functional>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>
#include <exception>
#include<vector>
#include<array>
typedef std::pair<int, int> pair;


template <typename T>
class Vector
{
    // Your implementation of the Vector class starts here
    public:  

        int length;
        T* data;      
        // Default Constructor
        Vector():length(0),data(nullptr){};

        // Copy constructor
        Vector(const Vector& other)
        : Vector(other.length)
        {
            for (auto i=0; i<other.length; i++)
                data[i] = other.data[i];        
        }

        // Move Constructor
        Vector(Vector&& other)
        : length(other.length),data(other.data)
        {
            other.length = 0;
            other.data   = nullptr;
        }

        //Constructor 1
        Vector(const int l):length(l),data(new T [l]){};

        //Constructor 2
        Vector(const std::initializer_list<T>& list)
        : Vector((int)list.size())
        {
            std::uninitialized_copy(list.begin(), list.end(), data);
        }

        // destructor
        ~Vector()
        {
            delete[] data;
            data = nullptr;
            length = 0;
        }

        // Operators

        // Copy assignment
        Vector& operator=(const Vector& other)
        {
            if (this != &other)
            {
                delete[] data;
                length = other.length;
                data   = new double[other.length];
                for (int i=0; i<other.length; i++)
                    data[i] = other.data[i];
            }
            return *this;
        }

        // Move assignment
        Vector& operator=(Vector&& other)
        {
            if (this != &other)
            {
                delete[] data;
                length = other.length;
                data   = other.data;
                other.length = 0;
                other.data   = nullptr;
            }
            return *this;
        } 

        T& operator[](int i)const
        {
            return data[i];
        }

        // Operator+

        // Add another vector
        template<typename U>
        auto operator+(const Vector<U>& other) const -> Vector<decltype(data[0]+other.data[0])>
        {
            // Throw exception if the vectors have different length
            if (length!=other.length) throw "Vectors have different size!";
            else
            {
            Vector<decltype(data[0]+other.data[0])> v3(other.length);
            for (auto i=0; i<length; i++)
                v3.data[i] = data[i] + other.data[i];
            return v3;
            }
        }

        // // Operator-

        // // Subtract another vector
        template<typename U>
        auto operator-(const Vector<U>& other) const -> Vector<decltype(data[0]-other.data[0])>
        {
            // Throw exception if the vectors have different length
            if (length!=other.length) throw "Vectors have different size!";
            else
            {
            Vector<decltype(data[0]-other.data[0])> v3(other.length);
            for (auto i=0; i<length; i++)
                v3.data[i] = data[i] - other.data[i];
            return v3;
            }
        }

        template<typename A>
        auto operator*(const A scalar) -> Vector<decltype(data[0]*scalar)>
        {
            Vector<decltype(data[0]*scalar)> v3(length);
            for (auto i=0; i<length; i++)
                v3.data[i] = data[i] * scalar;
            return v3;
        }
};

template<typename A,typename B>
        auto operator*(const A scalar,Vector<B> vec) -> Vector<decltype(vec.data[0]*scalar)>
        {
            Vector<decltype(vec.data[0]*scalar)> v3(vec.length);
            for (auto i=0; i<vec.length; i++)
                v3.data[i] = vec.data[i] * scalar;
            return v3;
        }

template<typename A>
        int len(Vector<A> vec)
        {
            return vec.length;
        }

template<typename T, typename U>
typename std::common_type<T,U>::type 
dot(const Vector<T>& lhs, 
    const Vector<U>& rhs)
{
    // Your implementation of the dot function starts here
    if (lhs.length!=rhs.length) throw "Vectors have different size!";
    else
    {
        auto sum = 0;
        for (auto i=0; i<lhs.length; i++)
            sum = sum + lhs.data[i]*rhs.data[i];
        return sum;
    }

}

template<typename T>
T norm(const Vector<T>& vec)
{
    // Your implementation of the norm function starts here
    auto norm = 0;
    for (auto i=0; i<vec.length; i++)
        norm = norm + pow(vec.data[i],2);
    return sqrt(norm);
}

template <typename T>
class Matrix
{
    // Start your implementation of the matrix class here
public:
    Matrix(const int r, const int c)
    {
        T M[r][c];
        std::map<pair,T> Map;
    };

    T& operator[] (const pair& ij)
    {   
        return Map[ij];
    }


     T operator()(const pair& ij) const
    {
        // Throw exception if the vectors have different length
        auto it = Map.find(ij);
        if (it == Map.end()) throw "Element is empty!";
        
        return it->second;
        
    }

    // // destructor
    // ~Matrix()
    // {
    //     delete[] M;
    //     M = nullptr;
    // }

// private:
    std::map<pair,T> Map;
    T M;

};

template<typename T>
int bicgstab(const Matrix<T>& A, 
             const Vector<T>& b, 
             Vector<T>&       x, 
             T                tol     = (T)1e-8, 
             int              maxiter = 100)
{
    // Your implementation of the bicgstab function starts here
    return 0;
}

template<typename T>
void heun(const Vector<std::function<T(Vector<T> const&, T)> >& f,
          Vector<T>&                                            y, 
          T                                                     h,
          T&                                                    t)
{
    // Your implementation of the heun function starts here
};

template<typename T>
class SimplestWalker
{
    // Your implementation of the simplest walker class starts here
};

int main(int argc, char* argv[])
{
    // Your testing of the simplest walker class starts here
    Matrix<double> M(10, 20); 
    Vector<int> v1 = {8,7,5};
    Vector<float> v2 = {1.0,2.0,3.0};
    // Vector v3 = v1-v2;
    // v3 = v2+v1;

    // M[{0,0}] = 1.0; // set value at row 0, column 0 to 1.0
    // M[{1,2}] = 2.0;

    // std::cout << M[{0,0}] << std::endl; // prints 1.0
    // std::cout << M({3,3}) << std::endl;
    std::cout << norm(v2) << std::endl;

    return 0;
}
