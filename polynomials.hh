#include<cmath>
#include<iostream>
#include<vector>
#include<sstream>
#include<limits>

/*
	Important.
	----------
	A polynomial is represented by a list of its coefficients in the usual order.
	E.g. x^2 -5x + 6 is represented by { 1, -5, 6 }. More specifically, such list
	is of the form vector<double>, such as
	vector<double> p_x = { 1, -5, 6 };
		  
*/

typedef unsigned int uint;

/*
*	Generates a string representing the list of coefficients of the polynomial.
	Ex: vector<double> p = {3,2,1};
		cout << vector_to_string(p); 
		Output: [ 3 2 1 ]
*/
std::string vector_to_string(std::vector<double> polynomial) {
	if (polynomial.empty()) return "[ ]";
	std::ostringstream os;
	os << "[ ";
	for (double& x : polynomial)
		os << x << " ";
	os << "]";
	return os.str();
}

/*
*	Generates a string representing a list of ordered pairs given by a 
*   variable of type vector<pair<double, double>>.
	Ex: vector<pair<double, double>> list;
		list.push_back({1,2});
		list.push_back({2,3});
		cout << intervals_to_string(list);
		Output: [ (1, 2) (2, 3) ]
*/
std::string intervals_to_string(std::vector<std::pair<double, double>> intervals) {
	if (intervals.empty())
		return "[ ]";
	std::ostringstream os;
	os << "[ ";
	for (auto interval : intervals)
		os << "(" << interval.first << ", " << interval.second << ") ";
	os << "]";
	return os.str();
}

/*
*   Removes all leading 0 coefficients of a polynomial.
	Ex: vector<double> p = {0,0,3,0,0,1,0,0};
		trim(p);
		p --> {3,0,0,1,0,0}
*/
std::vector<double> trim(std::vector<double> polynomial) {
	std::vector<double> trimmed_polynomial = polynomial;
	if (polynomial.empty()) return trimmed_polynomial;
	auto i = trimmed_polynomial.begin();
	while (*i == 0) {
		trimmed_polynomial.erase(i);
		i = trimmed_polynomial.begin();
		if (trimmed_polynomial.empty())
			break;
	}
	return trimmed_polynomial;
}

/*
	Adds n zeros (coefficientes) to the right of a polynomial.
	Ex: vector<double> p = { 3, 2, 1 };
	p = shift_left(p, 2);
	p --> { 3, 2, 1, 0, 0 }
*/
std::vector<double> shift_left(std::vector<double> polynomial, uint n)
{
	std::vector<double> zeroes(n, 0);
	polynomial.insert(polynomial.end(), zeroes.begin(), zeroes.end());
	return polynomial;
}

/*
	Adds n zeros (coefficients) to the left of a polynomial.
	Ex: vector<double> p = { 3, 2, 1 };
	p = shift_right(p, 3);
	p --> { 0, 0, 0, 3, 2, 1 }
*/
std::vector<double> shift_right(std::vector<double> polynomial, uint n)
{
	std::vector<double> zeroes(n, 0);
	polynomial.insert(polynomial.begin(), zeroes.begin(), zeroes.end());
	return polynomial;
}

/*
	Gives the degree of a polynomial.
	Ex: cout << degree({ 1, 0, 0});
	    Output: 2
*/
uint degree(std::vector<double> polynomial)
{
	polynomial = trim(polynomial);
	return polynomial.size() - 1;
}

/*
	Returns a polynomial which is the derivative of the argument.
	Ex: derivative({3,2,1}) --> {6,2}
*/
std::vector<double> derivative(std::vector<double> polynomial)
{
	polynomial = trim(polynomial);
	uint degree_p = degree(polynomial);
	if (degree_p == 0)
		return { 0 };
	std::vector<double> derivative_p(degree_p, 0);
	for (uint i = 0; i < degree_p; i++)
		derivative_p[i] = (degree_p - i) * polynomial[i];
	return derivative_p;
}

/*
	Overloaded addition operator allowing direct sum of polynomials,
	i.e. sum between vector<double> and vector<double>.

	TODO: don't trim. Check degree and shift. Only sum from end of 
	smaller degree to beginning of smaller degree.
*/
std::vector<double> operator+(std::vector<double> p, std::vector<double> q)
{
	p = trim(p); q = trim(q);
	int degree_gap = p.size() - q.size();
	if (degree_gap > 0)
		q = shift_right(q, degree_gap);
	else if (degree_gap < 0)
		p = shift_right(p, -degree_gap);
	std::vector<double> sum(p);
	for (uint i = 0; i < p.size(); i++)
		sum[i] += q[i];
	return sum;
}

std::vector<double> operator+(std::vector<double> p, double a)
{
	p.back() += a;
	return p;
}

std::vector<double> operator+(double a, std::vector<double> p)
{
	return p + a;
}

/*
	Overloaded negative operator. Inverts the signal of every term of a polynomial p;
*/
std::vector<double> operator-(std::vector<double> p)
{
	for (auto& pii : p)
		pii = -pii;
	return p;
}

/*
	Overloaded subtraction operator allowing direct sum of polynomials,
	i.e. subtraction between vector<double> and vector<double>.
*/
std::vector<double> operator-(std::vector<double> p, std::vector<double> q)
{
	p = trim(p); q = trim(q);
	int degree_gap = p.size() - q.size();
	if (degree_gap > 0)
		q = shift_right(q, degree_gap);
	else if (degree_gap < 0)
		p = shift_right(p, -degree_gap);
	std::vector<double> sum(p);
	for (uint i = 0; i < p.size(); i++)
		sum[i] -= q[i];
	return sum;
}

std::vector<double> operator-(std::vector<double> p, double a)
{
	p.back() -= a;
	return p;
}

std::vector<double> operator-(double a, std::vector<double> p)
{
	return -p + a;
}

/*
	Overloaded scalar multiplication operator allowing direct product of a polynomials 
	by a scalar, i.e. product between vector<double> and double.
*/
std::vector<double> operator*(std::vector<double> v, double a)
{
	for (auto& x : v) x *= a;
	return v;
}

/*
	Overloaded scalar multiplication operator allowing direct product of a polynomials
	by a scalar, i.e. product between double and vector<double>.
*/
std::vector<double> operator*(double a, std::vector<double> v)
{
	return v * a;
}

std::vector<double> operator/(std::vector<double> v, double a)
{
	return v * (1. / a);
}

std::vector<double> monomial_multiplication(std::vector<double> v, double u, uint degree_u)
{
	uint degree_v = degree(v);
	uint new_degree = degree_v + degree_u;
	v = trim(v);
	std::vector<double> product(new_degree + 1);
	for (uint i = 0; i < degree_v + 1; i++)
		product[i] = v[i] * u;
	return product;
}

std::vector<double> operator*(std::vector<double> p, std::vector<double> q)
{
	uint new_degree = 0;
	uint degree_p = degree(p);
	uint degree_q = degree(q);
	p = trim(p);
	q = trim(q);
	for (uint k = 0; k < degree_q + 1; k++)
		new_degree += (q[k] != 0) ? k : 0;
	std::vector<double> product(degree_p + new_degree + 1);
	for (uint k = 0; k < degree_q + 1; k++)
	{
		if (q[k] != 0)
			product = product + monomial_multiplication(p, q[k], degree_q - k);
	}
	return product;
}

/*
	General algorithm for polynomial divison. It was found as a refined version
	on StackOverflow (exact source missing).
	Ex: vector<double> f = {3,-6,13,-9,11,-1};
		vector<double> g = {1,-2,3};
		pair<vector<double>, vector<double>> result = polynomial_division(f, g);
		result.first --> {3, 0, 4, -1}
		result.second --> {-3, 2}
*/
std::pair<std::vector<double>, std::vector<double>> polynomial_division(std::vector<double> num, std::vector<double> den)
{
	num = trim(num);
	den = trim(den);

	int shift_degree = num.size() - den.size();
	if (shift_degree >= 0)
	{
		// Shift den towards left so it's the same degree
		for (uint i = 0; i < shift_degree; i++)
			den.push_back(0.);
	}
	else
		return { {0}, num };

	std::vector<double> quotient;
	double divisor = den[0];
	for (uint i = 0; i < shift_degree + 1; i++)
	{
		// Get the next coefficient of the quotient
		double mult = num[0] / divisor;
		quotient.push_back(mult);

		/* Subtract mult * den from num, but don't bother if mult == 0 
		   Note that when i==0, mult!=0; so quotient is automatically trimmed*/
		if (mult != 0)
			num = num - mult * den;

		num.erase(num.begin());
		den.pop_back();
	}
	return { quotient, num };
}

/*
	Overloaded division operator allowing direct division of a polynomials 
	by each other, i.e. division between vector<double> and vector<double>.
	Ex: vector<double> f = {3,-6,13,-9,11,-1};
		vector<double> g = {1,-2,3};
		vector<double> h = f / g;
		cout << vector_to_string(h);
		Output: [ 3, 0, 4, -1 ]
*/
std::vector<double> operator/(std::vector<double> num, std::vector<double> den)
{
	return polynomial_division(num, den).first;
}

/*
	Overloaded modulo operator allowing modulo operation of a polynomials
	by each other, i.e. division between vector<double> and vector<double>.
	Ex: vector<double> f = {3,-6,13,-9,11,-1};
		vector<double> g = {1,-2,3};
		vector<double> h = f % g;
		cout << vector_to_string(h);
		Output: [ -3, 2 ]
*/
std::vector<double> operator%(std::vector<double> num, std::vector<double> den)
{
	return polynomial_division(num, den).second;
}

/*
	Utilizes Horner's method to determine the value of a polynomial P(x) at x.
*/
double evaluate(std::vector<double> v, double x)
{   // Horner's method is used
	uint n = v.size();
	double sum = 0;
	for (uint i = 0; i < n; i++)
		sum = sum * x + v[i]; 
	return sum;
}

/*
	Returns the number of real roots of a polynomial in the interval [a, b].
	We utilize Sturm's theorem in order to make this calculation.
	Ex: cout << number_of_roots({ 1, -5, 6 }, 0, 4);
		Output: 2
*/
uint number_of_roots(std::vector<double> polynomial, double a, double b)
{   // We generate the Sturm sequence of polynomials
	std::vector<std::vector<double>> g;
	g.push_back(polynomial);
	g.push_back(derivative(polynomial));
	for (uint i = 0, j = 1; g[j].size() > 0; i++, j++)
		g.push_back(-1 * g[i] % g[j]);
	g.pop_back();

	// And use it to determine v(a) and v(b)
	uint n = g.size();
	int v_a = 0, v_b = 0;
	int previous_a = (evaluate(g[0], a) > 0) ? 1 : -1;
	int	previous_b = (evaluate(g[0], b) > 0) ? 1 : -1;;
	for (uint i = 1; i < n; i++)
	{   // we watch for a sign change
		double g_at_a = evaluate(g[i], a);
		double g_at_b = evaluate(g[i], b);
		if (g_at_a * previous_a < 0)
			v_a++;
		if (g_at_b * previous_b < 0)
			v_b++;
		// we can't update the previous value if it's zero
		if (g_at_a != 0)
			previous_a = (g_at_a > 0) ? 1 : -1;
		if (g_at_b != 0)
			previous_b = (g_at_b > 0) ? 1 : -1;
	}
	return v_a - v_b;
}

/*
	Returns a list of ordered pairs [a, b] such that each pair represents
	an interval that contains an isolated roots of the polynomial. It's a 
	recursive function of the number of roots (n) between a and b, with two
	base cases:
	1. n = 0
		Returns an empty list.
	2. n = 1
		Returns a list containing a single pair. 
	3. n > 1
		Chops the interval in half and searches separately on each half,
		accumulating the results if any root is isolated.
	Note: The user can specify a maximum size of interval (tolerance), so that
	even if only a single is found within an interval [s, t], the algorithm 
	keeps searching for a smaller bound on the isolated root.
	Ex: vector<double> p = {1,-5,6,4,-8};
		vector<pair<double, double>> intervals = isolate_roots(p, -1.9, 3, 100);
		cout << intervals_to_string(intervals);
		Output: [ (-1.41, -0.92) (1.53, 2.02) ]
*/
std::vector<std::pair<double, double>> isolate_roots(std::vector<double> polynomial, double a, double b, double tolerance = 1)
{
	std::vector<std::pair<double, double>> list;
	uint n_roots = number_of_roots(polynomial, a, b);
	if (n_roots == 0) // base case
		return list;
	if (n_roots == 1) // base case
	{
		if (b - a > tolerance) // interval's too large, let's take 50% of it
		{
			double step = (b - a) / 2;
			for (uint i = 0; i < 2; i++)
			{
				// Determining subinterval limits
				double a2 = a + i * step, b2 = a2 + step;
				// Making a recursive call over the subinterval
				std::vector<std::pair<double, double>> subinterval = isolate_roots(polynomial, a2, b2, tolerance);
				// Picking the subinterval containing the root
				if (!subinterval.empty())
				{
					list.insert(list.end(), subinterval.begin(), subinterval.end());
					break;
				}
			}
		}
		else // Subinterval has appropriate size
			list.push_back({ a, b });

		return list;
	}
	
	double step = (b - a) / 2;
	for (uint i = 0; i < 2; i++)
	{
		// Determining subinterval limits
		double a2 = a + i * step, b2 = a2 + step;
		// Making a recursive call over the subinterval
		std::vector<std::pair<double, double>> subinterval = isolate_roots(polynomial, a2, b2, tolerance);
		// If the call returns non-empty, accumulate the result
		if (!subinterval.empty())
			list.insert(list.end(), subinterval.begin(), subinterval.end());
		// If we have found all roots, finish
		if (list.size() == n_roots) break;
	}

	return list;
}

/*
*	Implements the bisection method to shorten the interval (x1, x2) around the
*	root. It's identical to the algorithm on the slide seen in class.
*/
double bisection(std::vector<double> p, double x1, double x2, double epsilon, uint n_max) {
	double Px1 = evaluate(p, x1);
	double Px2 = evaluate(p, x2);
	if (Px1 * Px2 > 0) {
		std::cout << "Error: function doesn't change sign between x1 and x2." << std::endl;
		return std::numeric_limits<double>::infinity();
	}
	double interval = abs(x2 - x1);
	double midpoint = x2;
	for (uint i = 0; i < n_max && interval * 5 > epsilon; i++) {
		midpoint = x1 + (x2 - x1) / 2;
		double Pmid = evaluate(p, midpoint);
		if (Px1 * Pmid > 0) {
			x1 = midpoint;
			Px1 = Pmid;
		}
		else {
			x2 = midpoint;
			Px2 = Pmid;
		}
		interval /= 2;
	}
	return midpoint;
}

/*
*	Searchs for a convergent root within an interval using the Newton Raphson method
*   for polynomials. It's identical to the algorithm on the slide seen in class.
*/
double newton_raphson(std::vector<double> polynomial , double x0, double epsilon1, double epsilon2, uint n_max) {
	double x = x0;
	double dx;
	int i;
	for (i = 0; i < n_max; i++) {
		double b = polynomial[0];
		double c = b;
		for (int j = 1; j < polynomial.size() - 1; j++) {
			b = b * x + polynomial[j];
			c = c * x + b;
		}
		b = b * x + polynomial.back();
		dx = b / c;
		x = x - dx;

		// para melhorar a precisão, ambas as condições de parada são avaliadas ao mesmo tempo.
		if (abs(b) < epsilon1 && abs(dx) < epsilon2) 
			return x;
	}
	std::cout << "Algorithm didn't converge after " << n_max << " iterations." << std::endl;
	std::cout << std::endl;
	return x;
}

/*
*	Returns the radius of a circular region that contains the roots of 
*	a polynomial. To do that, we apply the Second Circle Theorem.
	Ex: bound_on_roots({1,-1,1,-1}) --> 2
*/
double bound_on_roots(std::vector<double> polynomial) {
	polynomial = trim(polynomial);
	double radius = -std::numeric_limits<double>::infinity();
	double p_0 = polynomial.at(0);
	for (int i = 1; i < polynomial.size(); i++) { 
		double ratio = abs(polynomial.at(i)) / abs(polynomial.at(0));
		if (ratio > radius) {
			radius = ratio;
		}
	}
	return radius + 1;
}

/*
	Finds all the real convergent roots of a polynomial on the interval given by I.
	The interval argument is of the type vector<double> and must contain 2 elements,
	the lower and upper bound, respectively. If no such interval is provided,
	the efficient bound_on_roots function is called to find a suitable small
	interval of search. 
*/
std::vector<double> roots(std::vector<double> polynomial, double epsilon = 1e-6, std::vector<double> I = {}, uint n_max = 100, double tolerance = 1e-2)
{
	uint degree_p = degree(polynomial);
	double a = 0, b = 0, c = 0;
	std::vector<double> list_of_roots;
	switch (degree_p)
	{
	case 0:
		return { polynomial[0] };
	case 1:
		return { -polynomial[1] / polynomial[0] };
	case 2:
		a = polynomial[0];
		b = polynomial[1];
		c = polynomial[2];
		return { (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a),
				 (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a) };
	default:
		double radius = (I.empty()) ? bound_on_roots(polynomial) : I[1] - I[0];
		std::vector<std::pair<double, double>> intervals = isolate_roots(polynomial, -radius, radius, tolerance);
		for (auto& interval : intervals)
			list_of_roots.push_back(bisection(polynomial, interval.first, interval.second, epsilon, 5));
		for (auto& root : list_of_roots)
			root = newton_raphson(polynomial, root, epsilon, epsilon, n_max);
		return list_of_roots;
	}
}
/**/