#include<algorithm>
#include<cmath>
#include<iostream>
#include<map>
#include<memory>
#include<stack>
#include<sstream>
#include<string>
#include<vector>

typedef unsigned int uint;

std::vector<std::pair<std::string, double>> tokenize(std::string);
std::string print_expression(std::vector<std::pair<std::string, double>>, bool);
std::vector<std::pair<std::string, double>> to_postfix(std::vector<std::pair<std::string, double>>);
double evaluate_postfix(std::vector<std::pair<std::string, double>>); 
bool is_operator(std::pair<std::string, double>);

/*
Allows for the conversion of a string representing a simple mathematical expression into
a container that can be numerically evaluated later. It allows for late substitution of
variables. It accepts only the most common mathematical functions and operators:
sin, cos, tan, exp, log (the natural logarithm), +, -, *, / and ^
It processes a string composed of the above words and symbols plus numbers.
OBS1: Function arguments must be separated from function names by " " or "(",
i.e. do not write "logx". Correct: "log(x)" or "log x".
OBS2: Multiplication is always explicit, i.e. always write out the "*".
Example: 
	Expression example = Expression("a*x^2 + b*x + log(c)");
	std::cout << example.subs({ {"a", 1}, {"b", 1}, {"c", 2.718}, {"x", 2} }).evaluate();
Output:
	7
*/
class Expression
{

private:

	/* The mathematical expression itself will be stored as a list of tokens. Each
	   token is a pair of a string and a double. The expression will be stored in 
	   two forms: one in postfix order, useful for evaluation, and another in infix 
	   order, useful for printing. Tokens representing functions, operators and 
	   variables will store their values in the first member of the pair and will 
	   have value 0. 
	   Example: the cosine function is stored as: ("cos", 0.00000000000)
	   Tokens representing numbers will store their values in the second member of the
	   pair and will have label (key) "n".
	   Example: 2021 is stored as: ("n", 2021.0)
	   */
	std::vector<std::pair<std::string, double>> expression; // will be in postfix order
	std::vector<std::pair<std::string, double>> infix_expression; // useful for printing

public:

	void set_expression(std::string input)
	{
		this->infix_expression = tokenize(input);
		this->expression = to_postfix(this->infix_expression);
	}

	Expression(std::string input)
	{
		this->set_expression(input);
	}

	/* Returns a string representing the expression (infix notation). It utilizes the 
	   function print_expression which converts raw expressions to strings. It can also 
	   print the returned string to the console. */
	std::string to_string(bool output_to_console = false)
	{
		return print_expression(this->infix_expression, output_to_console);
	}

	/* Returns a list (vector) of the variables that appear in the given expression. 
	   The variables are represented in string format so that they're easy to print. */ 
	std::vector<std::string> variables()
	{
		uint n = this->expression.size();
		std::vector<std::string> variables = std::vector<std::string>();
		for (auto token : this->expression)
		{   /* We check if the token is not an operator or a number. If it's not, then
			   it must be a variable. */
			if (std::isalpha(token.first[0]) && !is_operator(token) && token.first != "n")
				/* We must only add the variable to the vector once. We search to see if it's there. */
				if (!std::binary_search(variables.begin(), variables.end(), token.first))
				{   /* We need to sort our vector of variables in order to do a binary search. */
					variables.push_back(token.first);
					std::sort(variables.begin(), variables.end());
				}
		}
		return variables;
	}

	/* Substitutes the specified variables for numeric values, then returns a new 
	   Expression object as a result. 
	   Parameter: 
		   - values: A map in which its keys contain the variable name and its values
		             are the desired values for each respective variable. 
	   Returns: 
		   - Expression; A new Expression with the desired symbols substituted. */
	Expression subs(std::map<std::string, double> values)
	{
		Expression expression = *this;
		for (auto &pair : expression.infix_expression)
		{
			for (auto value : values)
			{ // for each token in the infix expression we check if it's the desired variable.
				if (pair.first == value.first)
				{ // if it is, we make the substitution, making it into the specified number.
					pair.first = "n";
					pair.second = value.second;
				}
			}
		}
		expression.expression = to_postfix(expression.infix_expression);
		return expression;
	}

	/* Numerically evaluates the mathematical expression. Note that it only works
	   if all variables and unknowns have been substituted with numerical values. */
	double evaluate()
	{
		return evaluate_postfix(this->expression);
	}

	/* Returns a lambda function that encapsulates the process of substitution of variables
	   and numerical evaluation. If all you want is a C++ function that will return the 
	   numerical value of an arbitrary expression entered as a string by the user, this is the
	   simplest tool.
	   Parameter: 
	       - variable_names: a list (vector) of strings representing the names of the variables 
		     in the order you want them to be substituted/entered as arguments to the lambda 
			 function to be generated.
	   Returns:
	       - a lambda function with the following specifications:
	         Lambda function(std::vector<double> x)
		     Parameter: 
			     - x: a vector of doubles representing the numerical values to be substituted
			          into the expression.
		     Return value: 
			     - double; the numerical value of the original expression evaluated at the 
			               desired point.
		Example:
			Expression expr("alpha^2 + sin beta");
			// Note that we inverted the order of the variables!
			auto func = expr.lambdify({ "beta", "alpha" }); 
			cout << func({ 3.14159/2, 2 });
			Output: 3 
	*/
	auto lambdify(std::vector<std::string> variable_names)
	{   /* We need to make a shared pointer in order for the lambda function to be able to 
		   capture "this" expression as it will go out of scope. */
		auto expression = std::make_shared<Expression>(this->to_string());
		uint n = expression->variables().size();
		return [expression, variable_names, n](std::vector<double> x) {
			/* We'll be using subs() and evaluate() in order to numerically evaluate our expression
			   on the spot given input "x". In order to use subs() we must pass to it as argument 
			   a map of the form {(variable, numerical value)}. */
			std::map<std::string, double> variable_values = std::map<std::string, double>();
			for (uint i = 0; i < n; i++)
				variable_values.insert({ variable_names[i], x[i] });
			return expression->subs(variable_values).evaluate();
		};
	}
};

/* Returns true if a char is not a digit, a letter or '.' */
bool isother(char c)
{
	return !(std::isdigit(c) || std::isalpha(c)) && c != '.';
}

/* Adds whitespace between different elements of a mathematical expression in a string.
   Useful for later splitting the string into tokens using whitespace as the delimiter. */
std::string format_spacing(std::string input)
{
	std::string output = "";
	char aux = ' ';
	for (auto character : input)
	{ // whenever the next char changes from number to letter, or to operator etc, add an whitespace.
		if (std::isalpha(aux) && std::isdigit(character)
			|| std::isalpha(character) && std::isdigit(aux)
			|| std::isalpha(aux) && isother(character)
			|| isother(aux) && std::isalpha(character)
			|| std::isdigit(aux) && isother(character)
			|| isother(aux) && std::isdigit(character)
			|| isother(aux) && isother(character))
			output += ' ';
		output += character;
		aux = character;
	}
	return output;
}

/* Creates a sequence of tokens of elements of a mathematical expression out of a string.
   Each token is implemented as a pair with first element string and second element double:
   first element stores functions, variables and operators and second element stores numbers.
   Returns a vector of such pairs in the order given by the original string. */
std::vector<std::pair<std::string, double>> tokenize(std::string input)
{
	std::string formatted_input = format_spacing(input);
	std::istringstream ss(formatted_input);
	std::string aux;
	std::vector<std::pair<std::string, double>> tokenized_expression = std::vector<std::pair<std::string, double>>();
	while (std::getline(ss, aux, ' '))
	{
		if (std::isdigit(aux[0]))
			tokenized_expression.push_back({ "n", std::stod(aux) });
		else if (aux == "+")
			tokenized_expression.push_back({ "+", 0 });
		else if (aux == "-")
			tokenized_expression.push_back({ "-", 0 });
		else if (aux == "*")
			tokenized_expression.push_back({ "*", 0 });
		else if (aux == "/")
			tokenized_expression.push_back({ "/", 0 });
		else if (aux == "^")
			tokenized_expression.push_back({ "^", 0 });
		else if (aux == "(")
			tokenized_expression.push_back({ "(", 0 });
		else if (aux == ")")
			tokenized_expression.push_back({ ")", 0 });
		else if (aux.size() > 0)
			tokenized_expression.push_back({ aux, 0 });
	}
	return tokenized_expression;
}

/* Takes a vector of tokens from a mathematical expression and rebuilds a string 
   that represents that expression. It also directly prints the result to the console
   if the optional argument "output_to_console" is not set to false. */
std::string print_expression(std::vector<std::pair<std::string, double>> expression, bool output_to_console = true)
{
	uint n = expression.size();
	std::string expression_string = "";
	for (uint i = 0; i < n; i++)
	{
		if (expression[i].first == "n")
			expression_string += std::to_string(expression[i].second);
		else
			expression_string += expression[i].first;
		expression_string += " ";
	}
	if (output_to_console == true)
		std::cout << expression_string << std::endl;
	return expression_string;
}

/* Gives the precedence of an operator or function. */
uint precedence(std::string op)
{
	if (op == "+" || op == "-") return 1;
	else if (op == "*" || op == "/") return 2;
	else if (op == "^") return 3;
	else if (op == "sin" || op == "cos" || op == "tan" || op == "log" || op == "exp") return 4;
	else return 0;
}

/* Converts a mathematical expression from infix notation to postfix notation. It uses a standard 
   algorithm found everywhere on the internet. Works on an expression made into tokens according
   to the above conventions: a vector of pairs { string, double }. It returns another vector
   of pairs { string, double }, but in postfix order. */
std::vector<std::pair<std::string, double>> to_postfix(std::vector<std::pair<std::string, double>> expression)
{
	uint n = expression.size();
	std::vector<std::pair<std::string, double>> postfix_expression = std::vector<std::pair<std::string, double>>();
	std::stack<std::pair<std::string, double>> operator_stack = {};
	operator_stack.push({ "#", 0 }); 
	for (auto token : expression)
	{
		if (std::isalpha(token.first[0]))
		{
			if (token.first.size() == 1)
				postfix_expression.push_back(token);
			else
				operator_stack.push(token);
		}
		if (std::isdigit(token.first[0]))
			postfix_expression.push_back(token);
		else if (token.first == "(")
			operator_stack.push(token);
		else if (token.first == "^")
			operator_stack.push(token);
		else if (token.first == ")")
		{
			while (operator_stack.top().first != "#" && operator_stack.top().first != "(")
			{
				postfix_expression.push_back(operator_stack.top());
				operator_stack.pop();
			}
			operator_stack.pop();
		}
		else if (token.first == "+" || token.first == "-" || token.first == "*" || token.first == "/")
		{
			if (precedence(token.first) > precedence(operator_stack.top().first))
				operator_stack.push(token);
			else
			{
				while (operator_stack.top().first != "#" && precedence(token.first) <= precedence(operator_stack.top().first))
				{
					postfix_expression.push_back(operator_stack.top());
					operator_stack.pop();
				}
				operator_stack.push(token);
			}
		}
	}
	while (operator_stack.top().first != "#")
	{
		postfix_expression.push_back(operator_stack.top());
		operator_stack.pop();
	}
	return postfix_expression;
}

/* Returns true if a token codifies a binary operator. */
bool is_binary_operator(std::pair<std::string, double> token)
{
	return token.first == "+" || token.first == "-" || token.first == "/"
		|| token.first == "*" || token.first == "^";
}

/* Returns true if a token codifies an unary operator. */
bool is_unary_operator(std::pair<std::string, double> token)
{
	return token.first == "sin" || token.first == "cos"
		|| token.first == "tan" || token.first == "log"
		|| token.first == "exp";
}


/* Returns true if a token codifies an operator. */
bool is_operator(std::pair<std::string, double> token)
{
	return is_binary_operator(token) || is_unary_operator(token);
}

/* Carries out all the arithmetic operations on a mathematical expression, returning
   a single numerical result from the expression. Works on expressions made into tokens 
   according to the conventions: a vector of pairs { string, double }. It uses a
   standard algorithm found everywhere on the internet. 
   WARNING: Will not make sense for expressions containing free variables or unknowns. */
double evaluate_postfix(std::vector<std::pair<std::string, double>> expression)
{
	uint n = expression.size();
	std::stack<double> _stack = {};
	for (uint i = 0; i < n; i++)
	{
		std::pair<std::string, double> token = expression[i];
		if (!is_operator(token))
		{
			_stack.push(token.second);
		}
		else
		{
			if (is_binary_operator(token))
			{
				double rhs = _stack.top();
				_stack.pop();
				double lhs = _stack.top();
				_stack.pop();

				if (token.first == "+") _stack.push(lhs + rhs);
				if (token.first == "-") _stack.push(lhs - rhs);
				if (token.first == "*") _stack.push(lhs * rhs);
				if (token.first == "/") _stack.push(lhs / rhs);
				if (token.first == "^") _stack.push(std::pow(lhs, rhs));
			}
			else // it's an unary operator
			{
				double operand = _stack.top();
				_stack.pop();

				if (token.first == "sin") _stack.push(std::sin(operand));
				if (token.first == "cos") _stack.push(std::cos(operand));
				if (token.first == "tan") _stack.push(std::tan(operand));
				if (token.first == "log") _stack.push(std::log(operand));
				if (token.first == "exp") _stack.push(std::exp(operand));
			}
		}
	}
	return _stack.top();
}