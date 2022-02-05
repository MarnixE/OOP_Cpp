#include <cmath>
#include <functional>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

template <typename T>
class Vector
{
    // Your implementation of the Vector class starts here
    public:        
        // Default Constructor
        Vector():length(0),data(nullptr);

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
        Vector(const int l):length(l),data(new T [l]);

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

        // Operator+

        Vector operator+(const Vector& other) const
        {
            Vector new_c(other.length);
            for (int i=0; i<other.length; i++)
            {
                new_c.data[i] = data[i]+other.data[i];
            }
            return new_c;
        }

    private:
        int length;
        T data;

};

template<typename T, typename U>
typename std::common_type<T,U>::type 
dot(const Vector<T>& lhs, 
    const Vector<U>& rhs)
{
    // Your implementation of the dot function starts here
}

template<typename T>
T norm(const Vector<T>& vec)
{
    // Your implementation of the norm function starts here
}

template <typename T>
class Matrix
{
    // Start your implementation of the matrix class here
public:
    Matrix(int r, int c)
    {
        const std::map<T,T> M;
    }

    Matrix& operator[] (const std::pair<int, int>& ij)
    {
        if (M[ij] != 0)
        {
           return M[ij]; 
        }
    }

    // destructor
    ~Matrix()
    {
        delete[] M;
        M = nullptr;
    }



};

template<typename T>
int bicgstab(const Matrix<T>& A, 
             const Vector<T>& b, 
             Vector<T>&       x, 
             T                tol     = (T)1e-8, 
             int              maxiter = 100)
{
    // Your implementation of the bicgstab function starts here
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

    return 0;
}
