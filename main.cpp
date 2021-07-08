#include "math_parser.hh"
#include "integrate.hh"
#include<limits>

/* Neste programa, vamos lançar mão do nossos header autorais:
   - integrate.hh, que implementa funções de integração numérica desenvolvidas ao longo da 
     disciplina. Note que o nosso método de Gauss-Legendre possui capacidade genérica para
	 quadratura com um número arbitrário de pontos de interpolação.
   - polynomials.hh, que é o primeiro trabalho realizado na disciplina de Métodos Numéricos I
     com o prof. Bento. O trabalho consistia em encontrar raízes de polinômios, e utilizaremos
	 essa funcionalidade aqui para achar raízes de polinômios de Legendre.
   - math_parser.hh, que utilizaremos para capturar expressões digitadas pelo usuário no console 
     e as transformaremos em funções chamáveis (lambda) de C++. O arquivo header tem descrições 
     extensas sobre sua funcionalidade em seus comentários. 
*/

using std::string; using std::vector; using std::pair;
using std::cin; using std::cout; using std::pow;
using std::endl; using std::getline; using std::sin;
using std::exp; using std::sqrt;

typedef double (*function)(double);

void COMPLETAR_TABELA_AULA_12()
{
	for (uint i = 2; i <= 10; i++)
	{
		if (i == 5)
			cout << "\nAlguns passos extras:\n";

		cout << "\nn = " << i << endl;
		auto hermite_quadrature = generate_quadrature(i, "Hermite");
		cout << "Hermite:\n" << hermite_quadrature.to_string();

		auto laguerre_quadrature = generate_quadrature(i, "Laguerre");
		cout << "Laguerre:\n" << laguerre_quadrature.to_string();

		auto chebyshev_quadrature = generate_quadrature(i, "Chebyshev");
		cout << "Chebyshev:\n" << chebyshev_quadrature.to_string();
	}
}

void teste_pre_definido()
{	
	for (uint number_of_points = 2; number_of_points <= 5; number_of_points++)
	{
		double test_result;
		/**/

		auto f1 = [](double x) { return sin(x) * (1 + sin(x)); };

		test_result = integrate_gaussian_quadrature(f1, "Hermite", number_of_points, 1e-6);
		cout << "Funcao: " << "e ^ (-x ^ 2) (sin(x) ^ 2 + sin(x))" << " - " << number_of_points << " pontos\n";
		cout << "Wolfram: I = " << 0.5602022593661119221195044350470 << endl;
		cout << "Gauss-Hermite:" << endl;
		cout << "I = " << test_result << endl << endl;

		auto f2 = [](double x) { return x * x * x; };

		test_result = integrate_gaussian_quadrature(f2, "Laguerre", number_of_points, 1e-6);
		cout << "Funcao: " << "e ^ (-x ) x^3" << " - " << number_of_points << " pontos\n";
		cout << "Wolfram: I = " << 6 << endl;
		cout << "Gauss-Laguerre:" << endl;
		cout << "I = " << test_result << endl << endl;

		auto f3 = [](double x) { return exp(x) * sin(x); };

		test_result = integrate_gaussian_quadrature(f3, "Chebyshev", number_of_points, 1e-6);
		cout << "Funcao: " << "(e ^ x sin(x)) / sqrt(1 - x ^ 2)" << " - " << number_of_points << " pontos\n";
		cout << "Wolfram: I = " << 1.55989 << endl;
		cout << "Gauss-Chebyshev:" << endl;
		cout << "I = " << test_result << endl << endl;
		/**/
	}
	cout << endl;
}

int main()
{
	//COMPLETAR_TABELA_AULA_12();

	cout << "Programa de integracao numerica.\n"
		<< "Tarefa 3 - Disciplina de Metodos Numericos 2021.1.\n\n";

	while (true)
	{
		/* Vamos pedir que o usuário digite uma expressão matemática de acordo com algumas regras e
		   guardá-la na forma de string para processá-la em uma função utilizável de C++. */
		cout << "Digite a expressao que deseja integrar numericamente abaixo; ou\n"
			<< "digite 'teste' para executar um teste pre-definido.\n"
			<< "Observacoes:\n"
			<< "- Digite 'sair' para encerrar.\n"
			<< "- Utilize apenas UMA variavel. Faremos integrais unidimensionais.\n"
			<< "- Variaveis podem ter qualquer nome, exceto palavras reservadas.\n"
			<< "- Palavras reservadas: sin, cos, tan, exp, log.\n"
			<< "- Operadores disponiveis: + - * / ^\n"
			<< "- Separe nomes de funcoes de nomes de variaveis com '(' ou ' '.\n"
			<< "- Exemplo: sin(x^3 + 1)^(cos(3.1415 + 1))\n"
			<< "- Exemplo: 0.5 * (32 + x) + log(2.718) + sin(x)/cos(x)\n"
			<< "\n: ";
		string entrada;
		getline(cin, entrada);
		if (entrada == "sair") break;

		if (entrada == "teste")
		{
			teste_pre_definido();
			continue;
		}
		/* Cria uma expressão matemática a partir de uma string (entrada do usuário). */
		Expression expressao(entrada);
		/* Extraímos o nome da variável utilizada na expressão que o usuário defininiu (neste programa,
		   esperamos apenas uma variável).*/
		vector<string> variavel = expressao.variables();
		/* Obtemos uma função lambda que avalia numericamente a expressão entrada pelo usuário.
		   Para tanto, precisamos especificar a(s) variável(eis) livre(s) que serão tomadas
		   como argumento da função na ordem em que desejamos que elas sejam tomadas. Como só
		   esperamos uma variável, a ordem não é importante.*/
		auto f_ = expressao.lambdify(variavel);
		/* A nova função lambda "f_" aceita um argumento do tipo vector. Por praticidade, vamos
		   encapsular essa função dentro de uma nova função "f" que aceita como argumento apenas
		   uma double. */
		auto f = [f_](double x) { return f_({ x }); };

		cout << "\nEscolha o metodo de quadratura. Opcoes:\nNC - Newton-Cotes\nGL - Gauss-Legendre";
		cout << "\nGH - Gauss-Hermite\nGLA - Gauss-Laguerre\nGC - Gauss-Chebyshev\n: ";
		string quadratura;
		cin >> quadratura;

		double x_i;
		double x_f;
		if (quadratura == "NC" || quadratura == "GL")
		{
			cout << "\nDigite o limite inferior:\n: ";
			cin >> x_i;

			cout << "\nDigite o limite superior:\n: ";
			cin >> x_f;
		}

		cout << "\nQuantos pontos da funcao deseja interpolar\n"
			<< "por particao?\n";
		if (quadratura == "NC")
			cout << "Opcoes:\n- 2, 3, 4 ou 5\n";
		if (quadratura == "GL" || quadratura == "GH" || quadratura == "GLA" || quadratura == "GC")
			cout << "Opcoes:\n- Qualquer inteiro positivo\n";
		uint numero_de_pontos;
		cin >> numero_de_pontos;

		string formula;
		bool formula_fechada;
		if (quadratura == "NC")
		{
			formula;
			cout << "\nEscolha o tipo de formula que deseja utilizar. Opcoes:\n"
				<< "- F: fechada\n- Qualquer outra resposta: aberta\n: ";
			cin >> formula;
			formula_fechada = (formula == "F") ? true : false;
		}

		/* Vamos salvar o número de iterações. */
		uint numero_de_iteracoes = 1;

		/* Finalmente chamamos a nossa função de integração. */
		double I = NAN;
		if (quadratura == "NC")
			I = integrate_newton_cotes(f, x_i, x_f, numero_de_pontos, formula_fechada, 1e-6, &numero_de_iteracoes);
		else if (quadratura == "GL")
			I = integrate_gauss_legendre(f, x_i, x_f, numero_de_pontos, 1e-6, &numero_de_iteracoes);
		else if (quadratura == "GH")
			I = integrate_gaussian_quadrature(f, "Hermite", numero_de_pontos);
		else if	(quadratura == "GLA") 
			I = integrate_gaussian_quadrature(f, "Laguerre", numero_de_pontos);
		else if (quadratura == "GC")
			I = integrate_gaussian_quadrature(f, "Chebyshev", numero_de_pontos);
		cout << "Resultado:\n" << I << endl;
		cout << "Numero total de iteracoes: " << numero_de_iteracoes << endl << endl;

		/* Limpando o buffer do cin() */
		cin.clear();
		cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
};