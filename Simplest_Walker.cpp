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
        T* data;
        int length;

        // Default Constructor
        Vector():data(nullptr),length(0){}
        
        //Constructor 1
        Vector(int l):data(new T[l]()),length(l) {}

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


        // Copy constructor
        Vector(const Vector& other)
        : Vector<T>(other.length)
        {
            for (auto i=0; i<other.length; i++)
                data[i] = other.data[i];        
        }

        // Move Constructor
        Vector(Vector&& other)
        :data(other.data),length(other.length)
        {
            other.data   = nullptr;
            other.length = 0;
        }

        // Operators

        // Copy assignment
        Vector& operator=(const Vector& other)
        {
            if (this != &other)
            {
                delete[] data;
                data   = new T[other.length]();
                length = other.length;
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
                if (length!=0)
                {
                    delete[] data;
                }
                length = other.length;
                data   = other.data;
                other.length = 0;
                other.data   = nullptr;               
            }
            return *this;
        } 

        T& operator[](const int i)const
        {
            return data[i];
        }

        // Operator+

        // Add another vector
        template<typename U>
        auto operator+(const Vector<U>& other) const-> Vector<decltype(data[0]+other.data[0])>
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

        // Operator-

        // Subtract another vector
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
        auto operator*(const A scalar) const-> Vector<decltype(data[0]*scalar)>
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
        int len(const Vector<A>& vec)
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
    Matrix(int r, int c) : rows(r), cols(c)
    {
        std::map<pair,T> Map;
    };

    T& operator[] (const pair& ij)
    {
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

    if (lhs.cols != rhs.length) throw "The Matrix and Vector are not compatible.";

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
    //auto q_0 = b-A*x_0;
    //auto v_0,p_0 = 0;
    //auto alpha, rho_0,omega_0 = 1;

    Vector<T> q;
    q.push_back(b-A*x[0]);
    Vector<T> r;
    r.push_back(b-A*x[0]);
    Vector<T> p;
    p.push_back(0);
    Vector<T> v;
    v.push_back(0);
    auto alpha=1;
    Vector<T> rho;
    rho.push_back(1);
    Vector<T>omega;
    omega.push_back(1);




    for (auto k=1; k<maxiter; k++)
    {
        auto rho_k = dot(q[0], r[k-1]);
        r.push_back(rho_k);
        auto beta = (rho[k]/rho[k-1])*(alpha/omega[k-1]);
        auto p_k = r[k-1] + beta*(p[k-1] - omega[k-1]*v[k-1]);
        p.push_back(p_k);
        v.push_back(A*p[k]);
        alpha = rho[k]/dot(q[0],v[k]);
        auto h = x[k-1] + alpha*p[k];

        if (norm(b-A*h)< tol)
        {
            x.push_back(h);
            return k;
        }
        auto s = r[k-1] - alpha * v[k];
        auto t = A*s;
        auto omega_k= dot(t,s) / dot(t,t);
        omega.push_back(omega_k);
        auto x_k= h+omega[k]*s;
        x.push_back(x_k);

        if (norm(b-A*x[k])<tol)
        {
            return k;
        }
        auto r_k = s-omega[k]*t;
        r.push_back(r_k);
    }

    return 0;
}
template<typename T>
void heun(const Vector<std::function<T(const Vector<T>&, T)> >& f,
          Vector<T>&                                            y, 
          T                                                     h,
          T&                                                    t) 
{
    // Your implementation of the heun function starts here
    Vector<double> s = {f[0](y, 0.009),f[1](y, 0.009),f[2](y, 0.009),f[3](y, 0.009)};

    // Vector<double> y_bar = y + h * s;

    // Vector<double> s_bar = {f[0](y_bar, 0.009),f[1](y_bar, 0.009),f[2](y_bar, 0.009),f[3](y_bar, 0.009)};
    // y = y + h/2 * (s + s_bar);
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
            [](Vector<double> const& y, double gamma) { return y[2]; },
            [](Vector<double> const& y, double gamma) { return y[3]; },
            [](Vector<double> const& y, double gamma) { return sin(y[1] - gamma) + pow(y[4],2) * sin(y[0]) - cos(y[1] - gamma) * sin(y[0]); },
            [](Vector<double> const& y, double gamma) { return sin(y[1] - gamma); } };

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
            pow(y[4],2) * sin(y[0]) - cos(y[1] - g) * sin(y[0]), sin(y[1] - g)};
        
        return y_dot;
    }

    const Vector<T>& step(T h) 
    {
        heun(f,y,h,t);
        return y;
    }
};

int main(int argc, char* argv[])
{
    Vector<double> y0 = {0.4,0.2,0.0,-0.2};
    Vector<double> y = {1.};
    Vector<double> x = {0.1,0.2,0.3,0.4};
    Vector<double> b = {0.2,0.4,0.-0.6,-0.9};
    // Matrix<double> A[4][4] = {0.0, 1.0 ,2.0 ,3.0 ,4.0 , 5.0 , 6.0 , 7.0 , 8.0 , 9.0 , 10.0 , 11.0  , 12.0, 13.0, 14.0, 15.0};

    // std::cout << bicgstab(A,b,x) << std::endl;

    // Vector<int> y1 = {10,10, 20 ,30};
    // Vector y3 = std::move(y1);
    double t0 = 0;
    double gamma = 0.009;
    double h = 0.001;

    y = y0;
    SimplestWalker<double> SW(y0, t0, gamma);

    Vector<double> y_dot = SW.derivative(y0);
    std::cout << SW.y[1] << std::endl;
    Vector<double> y_new = SW.step(h);
    std::cout << SW.y[1] << std::endl;
    return 0;
}