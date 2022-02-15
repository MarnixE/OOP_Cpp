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

        // Default Constructor
        Vector()
        {
            length = 0;
            data = nullptr;
        }
        
        //Constructor 1
        Vector(int len)
        {
            length = len;
            data = new T[len];
            
        }

        //Constructor 2
        Vector(const std::initializer_list<T>& list)
        : Vector((int)list.size())
        {
            std::uninitialized_copy(list.begin(), list.end(), data);
            length = list.size();
        }
        
        // Copy constructor
        Vector(const Vector& other): Vector(other.length)
        {
            for (auto i=0; i<other.length; i++)
                data[i] = other.data[i];        
        }

        // Move Constructor
        Vector(Vector&& other):length(other.length),data(other.data)
        {
            other.length = 0;
            other.data   = nullptr;
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
                if(length!=0)
                    delete[] data;
                data   = new T[other.length];
                length = other.length;
                for (int i=0; i<length; i++)
                    data[i] = other.data[i];
            }
            return *this;
        }

        // Move assignment
        Vector& operator=(Vector&& other)
        {
            if (this != &other)
            {
                if (length!=0)
                    delete[] data;
                data   = other.data;
                other.data   = nullptr;
                length = other.length;
                other.length = 0;             
            }
            return *this;
        } 

        const T& operator[](int i)const
        {
            if(i>= length)throw std::invalid_argument("Index is out of bounds!");
            T& d = this->data[i];
            return d;
        
        }

    
        T& operator[](int i)
        {
            if(i>= length)throw std::invalid_argument("Index is out of bounds!");
            T& d = this->data[i];
            return d;
        }

        // Operator+

        // Add another vector
        template<typename U>
        Vector<typename std::common_type<U,T>::type> operator+(const Vector<U>& other) const
        {
            // Throw exception if the vectors have different length
            if (length!=other.len()) throw std::invalid_argument("Vectors have different size!");
            else
            {
                Vector v3 = Vector(length);
                for (auto i=0; i<length; i++)
                    v3.data[i] = data[i] + other[i];
                return v3;
            }
        }

        // Operator-

        // Subtract another vector
        template<typename U>
        Vector<typename std::common_type<U,T>::type> operator-(const Vector<U>& other) const
        {
            // Throw exception if the vectors have different length
            if (length!=other.len()) throw std::invalid_argument("Vectors have different size!");
            else
            {
                Vector v3 = Vector(length);
                for (auto i=0; i<length; i++)
                    v3.data[i] = data[i] - other[i];
                return v3;
            }
        }

        template<typename A>
        Vector<typename std::common_type<A,T>::type> operator*(A scalar) const
        {
            Vector<typename std::common_type<A,T>::type> v3(length);
            for (auto i=0; i<length; i++)
                v3.data[i] = typename std::common_type<A,T>::type(data[i] * scalar);
            return v3;
        }
        
        int len()const
        {
            return length;
        }
    // private:        
            int length;
            T* data;

};

template<typename A,typename B>
        Vector<typename std::common_type<A,B>::type>  operator*(const A scalar,const Vector<B>& vec)
        {
            Vector v3 = Vector<typename std::common_type<A,B>::type> (vec.len());
            for (auto i=0; i<vec.len(); i++)
                v3.data[i] = typename std::common_type<A,B>::type(vec[i] * scalar);
            return v3;
        }





template<typename T, typename U>
typename std::common_type<T,U>::type 
dot(const Vector<T>& lhs, 
    const Vector<U>& rhs)
{
    // Your implementation of the dot function starts here
    if (lhs.len()!=rhs.len()) throw std::invalid_argument("Vectors have different size!");
    else
    {
        typename std::common_type<T,U>::type  sum = 0;
        for (auto i=0; i<lhs.len(); i++)
            sum = sum + lhs[i]*rhs[i];
        return sum;
    }
}

template<typename T>
T 
norm(const Vector<T>& vec)
{
    // Your implementation of the norm function starts here
    T norm = 0;
    for (auto i=0; i<vec.len(); i++)
        norm = norm + pow(vec[i],2);
    return sqrt(norm);
}

template <typename T>
class Matrix
{
    // Start your implementation of the matrix class here
public:
    Matrix(int r, int c) : rows(r), cols(c)
    {
        std::map<pair,T> Map;
    };

    T& operator[] (const pair& ij)
    {
        if (std::get<0>(ij) >= rows || std::get<1>(ij) >= rows) throw std::invalid_argument("Index is out of bounds!");

        return Map[ij];
    }


     const T& operator() (const pair& ij)  const
    {
        // Throw exception if the vectors have different length
        return Map.at(ij);
    }
    
    // destructor
    ~Matrix()
    {
        rows = 0;
        cols = 0;
    }

// private:
    std::map<pair,T> Map;
    int rows;
    int cols;

};

template<typename T, typename U>
Vector<typename std::common_type<T,U>::type>
operator* (const Matrix<T>& lhs, const Vector<U>& rhs)
{

    if (lhs.cols != rhs.length) throw std::invalid_argument("The Matrix and Vector are not compatible.");

    Vector<typename std::common_type<T,U>::type> v(lhs.rows);

    // for (typename std::map<pair, T>::iterator it = lhs.Map.begin(); it != lhs.Map.end(); it++)
    for (auto const& it : lhs.Map)
    {
        int i   = it.first.first;
        int j   = it.first.second;
        T value = it.second;

        v[j] = v[j] + value * rhs[i];    
    }
    return v;
}

template<typename T>
int bicgstab(const Matrix<T>& A, 
             const Vector<T>& b, 
             Vector<T>&       x, 
             T                tol     = (T)1e-8, 
             int              maxiter = 100)
{
    
    Vector<T> p_prev(b.len());
    Vector<T> v_prev(b.len());
    
    
    Vector<T> q_0= b-A*x;  
    Vector<T> r_prev= b-A*x;  
    
    T alpha=1;
    T rho;
    T beta;
    T rho_prev=1;
    T omega_prev=1;
    T omega;
    Vector<T> h;
    Vector<T> s;
    Vector<T> t;

    for (auto k=1; k<maxiter; k++)
    {
        // Vector<T> rho_prev = rho;
        rho = dot(q_0, r_prev);
        beta = (rho/rho_prev)*(alpha/omega_prev);
        p_prev = r_prev + beta*(p_prev - omega_prev*v_prev);
        v_prev = (A*p_prev);
        alpha = rho/dot(q_0,v_prev);
        h = x + alpha*p_prev;

        if (norm(b-A*h)< tol)
        {
            x = h;
            return k;
        }
        s = r_prev - alpha * v_prev;
        t = A*s;
        omega= dot(t,s) / dot(t,t);
        x = h+omega*s;

        if (norm(b-A*x)<tol)
        {
            return k;
        }
        r_prev = s-omega*t;
        omega_prev = omega;
        rho_prev=rho;
        // std::cout << norm(b-A*h) << std::endl;
    }

    return -1;
} 

template<typename T>
void heun(const Vector<std::function<T(const Vector<T>&, T)> >& f,
          Vector<T>&                                            y, 
          T                                                     h,
          T&                                                    t)
{
    // Your implementation of the heun function starts here

    Vector<T> s = y.len();
    s = {f[0](y, t),f[1](y, t),f[2](y, t),f[3](y, t)};
    
    Vector<T> y_bar = y + h * s;

    t = t + h;

    // const Vector<double> s_bar = {f[0](y, t),f[1](y, t),f[2](y, t),f[3](y, t)};

    // y = y + h/2 * (s_bar - s);

    Vector<T> s_bar = {f[0](y_bar, t),f[1](y_bar, t),f[2](y_bar, t),f[3](y_bar, t)};

    y = y + h/2 * (s + s_bar);
}


template<typename T>
class SimplestWalker
{
    // Your implementation of the simplest walker class starts here
public:
    Vector<T> y;
    T t;
    T g;
    const Vector<std::function<double(const Vector<double>&, double)> > f =
            {
            [](Vector<double> const& y, double t) { return y[2]; },
            [](Vector<double> const& y, double t) { return y[3]; },
            [=](Vector<double> const& y, double t) { return sin(y[1] - g) + pow(y[3],2) * sin(y[0]) - cos(y[1] - g) * sin(y[0]); },
            [=](Vector<double> const& y, double t) { return sin(y[1] - g); } };
    SimplestWalker(const Vector<T>& y0, 
                    T               t0, 
                    T               gamma)
        {
            y = y0;
            t = t0;
            g = gamma;
        }

    Vector<T> derivative(const Vector<T>& y) const 
    {
        Vector<double> y_dot = {y[2], y[3], sin(y[1] - g) + 
            pow(y[3],2) * sin(y[0]) - cos(y[1] - g) * sin(y[0]), sin(y[1] - g)};
        
        return y_dot;
    }

    const Vector<T>& step(T h) 
    {
        heun(f,y,h,t);

        y = y - h/2 * derivative(y);

        heun(f,y,h/2,t);
        return y;
    }
};

int main(int argc, char* argv[])
{
    Vector<double> y0 = {0.4,0.2,0.0,-0.2};
    Vector<int> y1 = {1,2,3,4};
    // Vector<double> y(8);
    Vector<double> x = {0.1,0.2,0.3,0.4};
    Vector<double> b = {2,4,0.6,0.8};
    Matrix<double> A(4,4);
    A[{0,0}] = 1;
    A[{1,1}] = 1;
    A[{2,2}] = 1;
    A[{3,3}] = 1;
    A[{2,3}] = -1;
    
    std::cout << bicgstab(A,b,x) << std::endl;

    double t0 = 0;
    double gamma = 0.009;
    double h = 0.001;
    SimplestWalker<double> SW(y0, t0, gamma);

    Vector<double> y_dot = SW.derivative(y0);
    while (SW.t < 2)
    {
        SW.step(h);

        std::cout << std::endl;
        std::cout << SW.y[0] << std::endl;
        std::cout << SW.y[1] << std::endl;
        std::cout << SW.y[2] << std::endl;
        std::cout << SW.y[3] << std::endl;
    }
    return 0;
}