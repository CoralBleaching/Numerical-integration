#include<iostream>
#include<fstream>
#include "polynomials.hh"

typedef unsigned int uint;

// NEWTON-COTES QUADRATURES

/* Each function implements a quadrature rule, up to five point interpolation (closed or open).
   Each takes as argument a function pointer and two doubles accounting for the limits of 
   integration. */

template <typename Function>
double closed_two_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	return delta_x / 2 * (f(x_i) + f(x_f));
}

template <typename Function>
double closed_three_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	double x_m = x_i + delta_x / 2;
	return delta_x / 6 * (f(x_i) + 4 * f(x_m) + f(x_f));
}

template <typename Function>
double closed_four_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	double h = delta_x / 3;
	return delta_x / 8 * (f(x_i) + 3 * f(x_i + h) + 3 * f(x_i + 2 * h) + f(x_f));
}

template <typename Function>
double closed_five_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	double h = delta_x / 4;
	return 2. / 45 * h * (7 * f(x_i) + 32 * f(x_i + h) + 12 * f(x_i + 2 * h) + 32 * f(x_i + 3 * h) + 7 * f(x_f));
}

template <typename Function>
double open_two_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	double h = delta_x / 3;
	return delta_x / 2 * (f(x_i + h) + f(x_i + 2 * h));
}

template <typename Function>
double open_three_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	double h = delta_x / 4;
	return delta_x / 3 * (2 * f(x_i + h) - f(x_i + 2 * h) + 2 * f(x_i + 3 * h));
}

template <typename Function>
double open_four_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	double h = delta_x / 5;
	return delta_x / 24 * (11 * f(x_i + h) + f(x_i + 2 * h) + f(x_i + 3 * h) + 11* f(x_i + 4 * h));
}

template <typename Function>
double open_five_point_rule(Function f, double x_i, double x_f)
{
	double delta_x = x_f - x_i;
	double h = delta_x / 6;
	return delta_x / 20 * (11 * f(x_i + h) - 14 * f(x_i + 2 * h) + 26 * f(x_i + 3 * h) - 14 * f(x_i + 4 * h) + 11 * f(x_i + 5 * h));
}


/* This is the basic integration routine as seen in the lecture.
   Parameters:
	   - f           a pointer to a function implementing the "f(x)" to be numerically integrated.
	   - x_i         lower limit of integration.
	   - x_f         upper limit of integration.
	   - quadrature  the number of interpolating points.
	   - closed      if true, a closed formula is used. Otherwise an open formula is used.
	   - precision   the desired precision of the computation.
	   - steps       optional pointer to an unsigned int that will store the total number of iterations.
	Returns:
	   - double; the result of integration. */
template <typename Function>
double integrate_newton_cotes(Function f, double x_i, double x_f, uint quadrature = 3, bool closed = true, double precision = 1e-6, uint* steps = nullptr, uint max_steps = 1e6)
{
	uint n = 1;
	if (steps != nullptr) (*steps) = 0;
	double I, I_0 = 0;
	double error = precision + 1;
	while (error > precision)
	{   // If the algorithm doesn't converge, let's stop it.
		if (steps != nullptr && *steps > max_steps)
		{
			std::cout << "Warning: integrate_newton_cotes - maximum number of steps reached.\n";
			break;
		}
		I = 0;
		double h = (x_f - x_i) / n;
		for (uint i = 0; i < n; i++)
		{
			if (steps != nullptr) (*steps)++;
			double x_0 = x_i + i * h;
			double x_1 = x_0 + h;
			switch (quadrature)
			{
			case 2:
				I += (closed) ? closed_two_point_rule(f, x_0, x_1) : open_two_point_rule(f, x_0, x_1);
				break;
			case 3:
				I += (closed) ? closed_three_point_rule(f, x_0, x_1) : open_three_point_rule(f, x_0, x_1);
				break;
			case 4:
				I += (closed) ? closed_four_point_rule(f, x_0, x_1) : open_four_point_rule(f, x_0, x_1);
				break;
			case 5:
				I += (closed) ? closed_five_point_rule(f, x_0, x_1) : open_five_point_rule(f, x_0, x_1);
				break;
			}
		}
		error = std::fabs(1 - I_0 / I);
		I_0 = I;
		n *= 2;
	}
	return I;
}

// GAUSS QUADRATURES

/*
	Utilizes a recurrence relation to calculate the coefficent of a Legendre polynomial (P_n)
	of arbitrary degree n. See https://proofwiki.org/wiki/Definition:Legendre_Polynomial.
	Returns a vector with the coefficients of P_n in the usual order. 
	Example: legendre_polynomial(4) --> { 4.375, 0, -3.75, 0, 0.375 }
	It stands for 4.375x^4 -3.75x^2 + 0.375, or (1/8) (35x^4 - 30x + 3).
*/
std::vector<double> legendre_polynomial(int n, std::vector<std::vector<double>>* p_sequence = nullptr)
{
	// First elements of the recurrence: P_0 = 1 and P_1 = x 
	std::vector<std::vector<double>> p = { {1}, {1,0} };
	for (int i = 1; i < n; i++)
	{
		std::vector<double> x = { 1, 0 };
		std::vector<double> p_next = (2. * i + 1) / (i + 1) * x * p[i] - (double) i / (i + 1) * p[i - 1];
		p.push_back(p_next);
	}

	if (p_sequence != nullptr)
		*p_sequence = p;
	return p[n];
}

/*
	Reverses the order of the coefficients in a polynomial representation.
	Example: vector<double> p = { 1 , -5, 6 }; // x^2 - 5x + 6
			 vector<double> q = reverse_polynomial(p);
			 cout << vector_to_string(q);
			 Output: [ 6 -5 1 ]
*/
std::vector<double> reverse_polynomial(std::vector<double> polynomial)
{
	uint n = polynomial.size();
	std::vector<double> reverse(n);
	for (uint i = 0; i < n; i++)
		reverse[i] = polynomial[n - 1 - i];
	return reverse;
}

/*
	Generates Hermite polynomials up to index m. It does so by using the recurrence
	relation seen in https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation
	to generate each coefficient. The user can provide a pointer to a vector of 
	polynomials (*vector<vector<double>>) in order to store the entire sequence, from
	H_0(x) up to H_m(x). Otherwise, the function only returns H_m(x).
*/
std::vector<double> hermite_polynomial(uint m, std::vector<std::vector<double>>* a_sequence = nullptr)
{
	std::vector<std::vector<double>> a = { {1}, {0, 2} };
	for (uint n = 2; n <= m; n++)
	{
		std::vector<double> p;
		for (uint k = 0; k < a.back().size(); k++)
		{
			double a_previous = (k < a.back().size() - 1) ? (k + 1) * a.back()[k + 1] : 0;
			double a_next = (k > 0) ? 2 * a.back()[k - 1] : 0;
			p.push_back(a_next - a_previous);
		}
		p.push_back(2 * a.back().back());
		a.push_back(p);
	}

	// The recurrence relation described in the source gives the polynomial in left to right order.
	for (auto& p : a)
		p = reverse_polynomial(p);

	if (a_sequence != nullptr)
		*a_sequence = a;
	return a[m];
}

/*
	Generates Laguerre polynomials up to intex n. It does so by using the recurrence
	relation seen in https://en.wikipedia.org/wiki/Laguerre_polynomials#Recursive_definition,_closed_form,_and_generating_function
	to generate each successive L_{k+1} in the sequence. The user can provide a pointer 
	to a vector of polynomials (*vector<vector<double>>) in order to store the entire 
	sequence, from L_0(x) up to L_m(x). Otherwise, the function only returns L_m(x).
*/
std::vector<double> laguerre_polynomial(uint n, std::vector<std::vector<double>>* p_sequence = nullptr)
{
	std::vector<std::vector<double>> p = { {1}, {-1,1} };
	
	for (uint k = 1; k < n; k++)
	{
		std::vector<double> x = { 1, 0 };
		std::vector<double> p_next = ((2 * k + 1 - x) * p[k] - k * p[k - 1]) / (k + 1);
		p.push_back(p_next);
	}
	if (p_sequence != nullptr)
		*p_sequence = p;
	return p[n];
}

/*
	We directly implement the definition of x as a function of alpha,
	which represents the change of variable from x in [x_i, x_f] to
	alpha in [-1, 1].
*/
double x_alpha(double alpha, double x_i, double x_f)
{
	return x_i / 2 + x_f / 2 + (x_f - x_i) / 2 * alpha;
}

/*
*	Returns the value of a L_i^(n)(alpha), an ingredient of the Lagrange 
	family of interpolating polynomials, corresponding to index i. The
	degree is implied from the number of interpolating points given
	in the vector nodes.
*/
double L_i(double i, std::vector<double> nodes, double alpha)
{
	uint n = nodes.size();
	double result = 1;
	for (uint k = 0; k < n; k++)
	{
		if (k == i) 
			continue;
		result *= (alpha - nodes[k]) / (nodes[i] - nodes[k]);
	}
	return result;
}

double factorial(int n)
{	
	for (int i = n - 1; i > 0; i--)
		n *= i;
	return (double)n;
}

/*
	Returns the value of the weight function for index i. It is a 
	straightfoward implementation of the mathematical formula for 
	w_i^(n)(alpha). The number of interpolating points used is 
	implied by the size of the argument 'roots'.
*/
std::vector<double> w_i(std::vector<double> roots, std::string quadrature = "Legendre")
{
	const double pi = 3.14159265358979323846;
	uint n = roots.size();
	std::vector<double> w(n);
	if (quadrature == "Hermite")
	{
		std::vector<double> H_nm1 = hermite_polynomial(n - 1);
		for (uint i = 0; i < n; i++)
			w[i] = std::pow(2, n - 1.) * factorial(n) * std::sqrt(pi) / ((double)n * n * evaluate(H_nm1, roots[i]) * evaluate(H_nm1, roots[i]));
	}
	else if (quadrature == "Laguerre")
	{
		std::vector<double> L_np1 = laguerre_polynomial(n + 1);
		for (uint i = 0; i < n; i++)
			w[i] = roots[i] / ((n + 1) * (n + 1) * evaluate(L_np1, roots[i]) * evaluate(L_np1, roots[i]));
	}
	else if (quadrature == "Chebyshev")
	{
		w = std::vector<double>(pi / n);
	}
	else // Default: Legendre quadraure
	{
		for (uint i = 0; i < n; i++)
		{
			auto L = [=](double x) {
				return L_i(i, roots, x);
			};
			w[i] = integrate_newton_cotes(L, -1, 1);
		}
	}

	return w;
}

/*
	Numerically integrates a function f over [x_i, x_f] using n interpolating
	points. It does so in a naive and inneficient way, by calculating all of
	the parameters of the quadrature before integration. The code is short and 
	neat, though.
*/
template<typename Function>
double naive_gauss_legendre(Function f, double x_i, double x_f, uint n)
{
	std::vector<double> p_n_roots = roots(legendre_polynomial(n));
	std::vector<double> w = w_i(roots);
	double sum = 0;
	for (uint i = 0; i < n; i++)
		sum += f(x_alpha(p_n_roots[i], x_i, x_f)) * w[i];
	return (x_f / 2 - x_i / 2) * sum;
}

/*
	Saves the parameters for a specific Gauss-Legendre quadrature of n interpolating
	points, in order to save computation. The user can save this quadrature to a file
	and import it later for further computations.
*/
class Quadrature
{
public:
	uint n;                     // The number of points in the quadrature.
	std::vector<double> alpha;  // The roots of the associated polynomial.
	std::vector<double> w;      // The weight function evaluated at the roots.

	// This constructor allows us to initialize a quadrature of size n with parameters equal to zero.
	Quadrature(uint n)
	{
		this->n = n;
		this->alpha = std::vector<double>(n);
		this->w = std::vector<double>(n);
	}

	// Facilites printing.
	std::string to_string()
	{
		std::ostringstream os;
		os << "alpha: " << vector_to_string(this->alpha) << std::endl;
		os << "w:     " << vector_to_string(this->w) << std::endl;
		return os.str();
	}
};

/*
	Generates the sufficient quadrature parameters for Gauss integration.
	The format it returns in is a self-explanatory class.
*/
Quadrature generate_quadrature(uint n, std::string name = "Legendre")
{
	const double pi = 3.14159265358979323846;
	Quadrature quadrature(n);
	std::vector<double> polynomial;
	if (name == "Legendre")
	{
		polynomial = legendre_polynomial(n);
		quadrature.alpha = roots(polynomial);
		quadrature.w = w_i(quadrature.alpha);
	}
	else if (name == "Laguerre")
	{
		polynomial = laguerre_polynomial(n);
		quadrature.alpha = roots(polynomial);
		quadrature.w = w_i(quadrature.alpha, "Laguerre");
	}
	else if (name == "Hermite")
	{
		polynomial = hermite_polynomial(n);
		quadrature.alpha = roots(polynomial);
		quadrature.w = w_i(quadrature.alpha, "Hermite");
	}
	else if (name == "Chebyshev")
	{
		for (uint i = 0; i < n; i++)
		{
			quadrature.alpha[n - i - 1] = std::cos((2. * (i + 1) - 1) / (2. * n) * pi);
			quadrature.w[i] = pi / n;
		}
	}
	return quadrature;
}

/* This is the basic integration routine as seen in the lecture.
   Parameters:
	   - f                 a pointer to a function implementing the "f(x)" to be numerically integrated.
	   - x_i               lower limit of integration.
	   - x_f               upper limit of integration.
	   - number_of_points  [optional] the number of interpolating points.
	   - precision         [optional] the desired precision of the computation.
	   - number_of_steps   [optional] optional pointer to an unsigned int that will store the total number of iterations.
	   - quadrature        [optional] if the user has pre-computed quadrature parameters, they can be passed via this argument.
	Returns:
	   - double; the result of the numerical integration. 
*/
template<typename Function>
double integrate_gauss_legendre(Function f, double x_i, double x_f, uint number_of_points = 3, double precision = 1e-6, uint* number_of_steps = nullptr, Quadrature quadrature = Quadrature(0), uint max_steps = 1e6)
{
	// We preload all the data for a quadrature of 'number_of_points' points.
	// That way we only do this calculation once. But we only do this if the
	// user hasn't already provided his own quadrature parameters.
	if (quadrature.n == 0)
		quadrature = generate_quadrature(number_of_points);

	if (number_of_steps != nullptr)
		*number_of_steps = 0;

	// Start the process of interval partitioning and integration
	uint n = 1;
	double I, I_0 = 0;
	double sum = 0;
	double error = precision + 1;
	while (error > precision)
	{   // If the algorithm doesn't converge, let's stop it.
		if (number_of_steps != nullptr && *number_of_steps > max_steps)
		{
			std::cout << "Warning: integrate_gauss_legendre - maximum number of steps reached.\n";
			break;
		}
		I = 0;
		double h = (x_f - x_i) / n;
		for (uint i = 0; i < n; i++)
		{
			if (number_of_steps != nullptr) 
				(*number_of_steps)++;
			double x_0 = x_i + i * h;
			double x_1 = x_0 + h;
			// Integrating on [x_0, x_1] 
			for (uint j = 0; j < number_of_points; j++)
				I += h / 2 * f(x_alpha(quadrature.alpha[j], x_0, x_1)) * quadrature.w[j];
		}
		error = std::fabs(1 - I_0 / I);
		I_0 = I;
		n *= 2;
	}
	return I;
}

/* This is the basic integration routine as seen in the lecture.
   Parameters:
	   - f                   a pointer to a function implementing the "f(x)" to be numerically integrated.
	   - type_of_quadrature  
	   - number_of_points    [optional] the number of interpolating points.
	   - precision           [optional] the desired precision of the computation.
	   - quadrature          [optional] if the user has pre-computed quadrature parameters, they can be passed via this argument.
	Returns:
	   - double; the result of the numerical integration. 
*/
template<typename Function>
double integrate_gaussian_quadrature(Function f, std::string type_of_quadrature, uint number_of_points = 3, double precision = 1e-6, Quadrature quadrature = Quadrature(0))
{
	// We preload all the data for a quadrature of 'number_of_points' points.
	// That way we only do this calculation once. But we only do this if the
	// user hasn't already provided his own quadrature parameters.
	if (quadrature.n == 0)
		quadrature = generate_quadrature(number_of_points, type_of_quadrature);
	// Integrating
	double I = 0;
	for (uint j = 0; j < number_of_points; j++)
		I += f(quadrature.alpha[j]) * quadrature.w[j];
	return I;
}