#include <iostream>
#include "boost_1_61_0/boost/math/tools/polynomial.hpp"

using std::string;

using namespace boost::math;
using namespace boost::math::tools;
using boost::lexical_cast;


template <typename T>
string sign_str(T const &x)
{
    return x < 0 ? "-" : "+";
}

template <typename T>
string inner_coefficient(T const &x)
{
    string result(" " + sign_str(x) + " ");
    if (abs(x) != T(1))
        result += lexical_cast<string>(abs(x));
    return result;
}

template <typename T>
string formula_format(polynomial<T> const &a)
{
    string result;
    if (a.size() == 0)
        result += lexical_cast<string>(T(0));
    else
    {
        // First one is a special case as it may need unary negate.
        unsigned i = a.size() - 1;
        if (a[i] < 0)
            result += "-";
        if (abs(a[i]) != T(1))
            result += lexical_cast<string>(abs(a[i]));

        if (i > 0)
        {
            result += "x";
            if (i > 1)
            {
                result += "^" + lexical_cast<string>(i);
                i--;
                for (; i != 1; i--)
                    if (a[i] != 0)
                        result += inner_coefficient(a[i]) + "x^" + lexical_cast<string>(i);
                if (a[i] != 0)
                    result += inner_coefficient(a[i]) + "x";
            }
            i--;

            result += " " + sign_str(a[i]) + " " + lexical_cast<string>(abs(a[i]));
        }
    }
    return result;
} // string formula_format(polynomial<T> const &a)

using polys = std::vector<polynomial<int>>;


polynomial<int> powered_x(int degree) {
    int size = degree + 1;
    int arr[size];

    for(int i = 0; i < size; ++i)
        arr[i] = 0;
    arr[degree] = 1;

    return polynomial<int>(arr, degree);
}

bool is_divisible(polynomial<int> a, polynomial<int> b) {
    polynomial<int> current(a);

    while (current.degree() > b.degree()) {
        polynomial<int> poly_to_sub = powered_x(current.degree() - b.degree()) * b;
        current -= poly_to_sub;

        for (auto& x : current.data()) {
            x = abs(x) % 2;
        }
    }

    return current == b;
}

void find_irreducible(polys current, polys& previous_irreducible) {
    for (auto& irreducible : previous_irreducible) {
        current.erase(
                std::remove_if(
                        current.begin(),
                        current.end(),
                        [&irreducible](polynomial<int> &p) {
                            p.normalize();
                            if (p.size() == 0|| p[0] == 0 || p.size() == 1 ) return true;

                            return is_divisible(p, irreducible);
                        }
                ),
                current.end()
        );
    }

    for (auto poly : current) {
        previous_irreducible.push_back(poly);
    }
}

void get_next_current(std::vector<std::vector<int>>& current_vectors, polys& current) {
    // cannot use current polynoms as they are normalised during initialisation
    std::vector<std::vector<int>> new_current_vectors;
    std::vector<polynomial<int>> new_current;

    for (auto& v : current_vectors) {
        std::vector<int> new_0(v);
        std::vector<int> new_1(v);

        new_0.push_back(0);
        new_1.push_back(1);

        new_current_vectors.push_back(new_0);
        new_current_vectors.push_back(new_1);

        new_current.push_back(polynomial<int>(new_0.begin(), new_0.end()));
        new_current.push_back(polynomial<int>(new_1.begin(), new_1.end()));
    }

    current_vectors = new_current_vectors;
    current = new_current;
}

int main() {
    int n;
    std::cout << "Enter the highest degree." << std::endl;
    std::cin >> n;
    std::cout << std::endl;
    std::cout << std::endl;

    std::vector<std::vector<int>> current_vectors({ {0, 0}, {0, 1}, {1, 0}, {1, 1} });
    std::vector<polynomial<int>> previous({ {{0, 1}}, {{1, 1}} });
    std::vector<polynomial<int>> current;
    get_next_current(current_vectors, current);

    std::cout << "---------- All irreducible of degree 1 ----------" << std::endl;
    for (auto& p : previous) {
        std::cout << formula_format(p) << std::endl;
    }

    for (int i = 2; i < n + 1; ++i) {
        find_irreducible(current, previous);

        std::cout << "---------- All irreducible of degree " << i << " ----------" << std::endl;
        int count = 0;
        for (auto& p : previous) {
            if (p.degree() == i) {
                std::cout << formula_format(p) << std::endl;
                count++;
            }
        }
        std::cout << "Count of polynomials: " << count << std::endl;

        get_next_current(current_vectors, current);
    }


    return 0;
}
// doubrov@bsu.by