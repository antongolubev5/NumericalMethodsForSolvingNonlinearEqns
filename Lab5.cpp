
//метод бисекции и ньютона, значения на концах знаков разные, 1 - локализация корня, с каким-то шагом разбивать отрезок
//выбирать кусочки, где корни и искать там, подумать, как проводить, с каким шагом, предложить рекомендации по разбиению
//2 - во  2 тесте должно получиться, что след х0 выйдет за пределы, предложить метод возврата на отрезок
//подумать, что делать в случае зацикливания, когда метод будет прыгать из точки в точку
//решаем для 2 уравнений, кодом обращать матрицу не нужно
//в зависимости от скорости сходимости показать цветом ???

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <iomanip>
#include <cmath>
using namespace std;

void printVector(double* Vector, int k)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		cout << Vector[i] /*<< " " << i+1 */ << endl;
	}

}

void printMatrix(double** extMatrix, int k, int m, bool extended)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < m; j++)
		{
			char simb = ' ';
			//if ((j == k - 1) && extended) { simb = '='; }
			//else { simb = ' '; }
			//cout <<  setprecision(2) << fixed<< extMatrix[i][j] << simb;
			cout << extMatrix[i][j] << simb;
		}
		cout << endl;
	}

}

double system1function1(double x, double y)
{
	return (x*x - y*y - 15);
}

double system1function2(double x, double y)
{
	return (x*y+4);
}

double system2function1(double x, double y)
{
	return (x*x + y*y +x+y-8);
}

double system2function2(double x, double y)
{
	return (x*x + y*y +x*y-7);
}

double system3function1(double x, double y)
{
	return (y-x*x);
}

double system3function2(double x, double y)
{
	return (x*x+y*y-4);
}

double system1function1derivedanalyticsx(double x, double y)
{
	return (2*x);
}

double system1function2derivedanalyticsx(double x, double y)
{
	return (y);
}

double system2function1derivedanalyticsx(double x, double y)
{
	return (2*x + 1);
}

double system2function2derivedanalyticsx(double x, double y)
{
	return (2*x + y);
}

double system3function1derivedanalyticsx(double x, double y)
{
	return (-2*x);
}

double system3function2derivedanalyticsx(double x, double y)
{
	return (2*x);
}

double system1function1derivedanalyticsy(double x, double y)
{
	return (- 2*y);
}

double system1function2derivedanalyticsy(double x, double y)
{
	return (x);
}

double system2function1derivedanalyticsy(double x, double y)
{
	return (2*y + 1);
}

double system2function2derivedanalyticsy(double x, double y)
{
	return (2*y + x);
}

double system3function1derivedanalyticsy(double x, double y)
{
	return (1);
}

double system3function2derivedanalyticsy(double x, double y)
{
	return (2*y);
}

double function1(double x)
{
	return (x - 0.1)*(x - 0.22)*(x - 0.55)*(x - 0.7)*(x - 0.75);
}

double function2(double x)
{
	return sqrt(x + 1) - 1;
}

double function3(double x)
{
	return 35 * x*x*x - 67 * x*x - 3 * x + 3;
}

double function4(double x)
{
	return (x - 1)*(x - 1);
}

double function5(double x)
{
	return x*x-1;
}

double function6(double x)
{
	return exp(x)-1;
}

double function1derivedanalytics(double x)
{
	return (x - 0.75)*(x - 0.7)*(x - 0.55)*(x - 0.22) + (x - 0.75)*(x - 0.7)*(x - 0.55)*(x - 0.1) +
		(x - 0.75)*(x - 0.7)*(x - 0.22)*(x - 0.1) + (x - 0.75)*(x - 0.55)*(x - 0.22)*(x - 0.1) + 
		(x - 0.7)*(x - 0.55)*(x - 0.22)*(x - 0.1);
}

double function2derivedanalytics(double x)
{
	return 1/(2*sqrt(x + 1));
}

double function3derivedanalytics(double x)
{
	return  105 * x*x - 134 * x -3;
}

double function4derivedanalytics(double x)
{
	return  2*(x-1);
}

double function5derivedanalytics(double x)
{
	return  2 *x;
}

double function6derivedanalytics(double x)
{
	return  exp(x);
}

double functionderivedcalculated(double x, double(*ff) (double), double eps)
{
	double result = (ff(x + eps) - ff(x))/eps;
	return result;
}

double Jacobian(double* vect, double(*derivative1) (double,double), double(*derivative2) (double,double),
	double(*derivative3) (double,double), double(*derivative4) (double,double))
{
	double result;
	double **MatrixJacobi;
	MatrixJacobi = new double*[2];

	for (int ii = 0; ii < 2; ii++)
	{
		MatrixJacobi[ii] = new double[2];
	}

	MatrixJacobi[0][0] = system1function1derivedanalyticsx(vect[0], vect[1]);
	MatrixJacobi[0][1] = system1function1derivedanalyticsy(vect[0], vect[1]);
	MatrixJacobi[1][0] = system1function2derivedanalyticsx(vect[0], vect[1]);
	MatrixJacobi[1][1] = system1function2derivedanalyticsy(vect[0], vect[1]);

	result = MatrixJacobi[0][0] * MatrixJacobi[1][1] - MatrixJacobi[0][1] * MatrixJacobi[1][0];

	return result;
}

double** JacobiMatrix(double* vect, double(*derivative1) (double, double), double(*derivative2) (double, double),
	double(*derivative3) (double, double), double(*derivative4) (double, double))
{
	double **MatrixJacobi;
	MatrixJacobi = new double*[2];

	for (int ii = 0; ii < 2; ii++)
	{
		MatrixJacobi[ii] = new double[2];
	}

	MatrixJacobi[0][0] = derivative1(vect[0], vect[1]);
	MatrixJacobi[0][1] = derivative2(vect[0], vect[1]);
	MatrixJacobi[1][0] = derivative3(vect[0], vect[1]);
	MatrixJacobi[1][1] = derivative4(vect[0], vect[1]);

	return MatrixJacobi;
}

double* multiplyMatrixVector(double** Matrix, double* Vector, int n)
{
	double *result;
	result = new double[n];
	double s = 0;

	for (int i = 0; i < n; i++)
	{
		s = 0;

		for (int j = 0; j < n; j++)
		{
			s += Matrix[i][j] * Vector[j];
		}
		result[i] = s;

	}
	return result;
}

double* substractVector(double* Vector1, double* Vector2, int n)
{
	double *result;
	result = new double[n];

	for (int i = 0; i < n; i++)
	{
		result[i] = Vector1[i] - Vector2[i];
	}

	return result;
}

double** RootLocalization(double(*ff) (double), double a, double b, double eps, int& npoints)
{
	//отталкиваемся от точности и выбираем шаг (точность eps = 1e-3)
	//double h = exp(log(eps) *0.5);
	double h = 0.01;
	double quantity = abs(a - b) / h;//количество промежутков на отрезке a,b
	int n = round(quantity + 1);
	//npoints = n;//для дальнейшего использования

	double **net;
	net = new double*[2];
	for (int ii = 0; ii<2; ii++)
		net[ii] = new double[n];
	double xi, yi;

	for (int i = 0; i < n; i++)
	{//создаем таблицу из аргументов и значений функции на отрезке с определенным выше шагом
		xi = a + h*i;
		yi = ff(xi);
		//cout << xi << ' ' << yi << endl;

		net[0][i] = xi;
		net[1][i] = yi;
	}
	int kol = 0;

	for (int i = 1; i < n; i++)
	{ 
		if ((net[1][i]>0 && net[1][i - 1]<0) 
			|| (net[1][i]<0 && net[1][i - 1]>0) 
			|| (net[1][i] == 0 && net[1][i - 1]>0 && net[0][i] == b)//правый конец-корень в точке в
			//|| (net[1][i] == 0 && net[1][i - 1]>0 && net[0][i - 1] == a)
			|| (net[1][i]>0 && net[1][i - 1]==0)//левый конец отрезка-корень
			//|| (net[1][i] == 0 && net[1][i - 1]<0 && net[0][i - 1] == a)
			|| (net[1][i] == 0 && net[1][i - 1]<0 && net[0][i] == b)//правый конец-корень в точке в
			|| (net[1][i]<0 && net[1][i - 1]==0)
			)
		{kol++;}
		else if (net[1][i] == 0 && net[1][i - 1] == 0)
		{
			kol++;
			kol++;
		}
		
	}

	npoints = kol;
	//создаем новую таблицу с концами отрезков, содержащих корни, размерность: 2 на kol
	double **rootstable;
	rootstable = new double*[2];
	for (int ii = 0; ii<2; ii++)
		rootstable[ii] = new double[kol];
	
	kol = 0;
	for (int i = 1; i < n; i++)
	{
		if ((net[1][i]>0 && net[1][i - 1]<0)
			|| (net[1][i]<0 && net[1][i - 1]>0)
			|| (net[1][i] == 0 && net[1][i - 1]>0 && net[0][i-1]==a)
			|| (net[1][i]>0 && net[1][i - 1] == 0)//ноль справа
			|| (net[1][i] == 0 && net[1][i - 1]<0 && net[0][i - 1] == a)
			|| (net[1][i]<0 && net[1][i - 1] == 0)//ноль справа
			)
		{
			rootstable[0][kol] = net[0][i - 1];
			rootstable[1][kol] = net[0][i];
			
			kol++;
		}

		else if (net[1][i] == 0 && net[1][i - 1] == 0)
		{
				rootstable[0][kol] = net[0][i - 1];
				rootstable[1][kol] = (net[0][i] + net[0][i + 1])*0.5;
				kol++;

				rootstable[0][kol] = (net[0][i] + net[0][i + 1])*0.5;
				rootstable[1][kol] = net[0][i];


				kol++;
		}

	}

	return rootstable;
}

double** RootLocalizationkrat(double(*ff) (double), double a, double b, double eps, int& npoints)
{
	//отталкиваемся от точности и выбираем шаг (точность eps = 1e-3)
	double h = exp(log(eps) *0.5);
	double quantity = abs(a - b) / h;//количество промежутков на отрезке a,b
	int n = round(quantity + 1);
	//npoints = n;//для дальнейшего использования

	double **net;
	net = new double*[2];
	for (int ii = 0; ii<2; ii++)
		net[ii] = new double[n];
	double xi, yi;

	for (int i = 0; i < n; i++)
	{//создаем таблицу из аргументов и значений функции на отрезке с определенным выше шагом
		xi = a + h * i;
		yi = ff(xi);
		//cout << xi << ' ' << yi << endl;

		net[0][i] = xi;
		net[1][i] = yi;
	}
	int kol = 0;

	for (int i = 1; i < n; i++)
	{
		if ((net[1][i]>1e-6 && net[1][i - 1]<1e-6) && (net[1][i-1]<1e-6 && net[1][i +1]>1e-6))
			//|| (net[1][i]<0 && net[1][i - 1]>0)
			//|| (net[1][i] == 0 && net[1][i - 1]>0 && net[0][i] == b)//правый конец-корень в точке в
			//														//|| (net[1][i] == 0 && net[1][i - 1]>0 && net[0][i - 1] == a)
			//|| (net[1][i]>0 && net[1][i - 1] == 0)//левый конец отрезка-корень
			//									  //|| (net[1][i] == 0 && net[1][i - 1]<0 && net[0][i - 1] == a)
			//|| (net[1][i] == 0 && net[1][i - 1]<0 && net[0][i] == b)//правый конец-корень в точке в
			//|| (net[1][i]<0 && net[1][i - 1] == 0)
			//)
		{
			kol++;
		}
		else if (net[1][i] == 0 && net[1][i - 1] == 0)
		{
			kol++;
			kol++;
		}

	}

	npoints = kol;
	//создаем новую таблицу с концами отрезков, содержащих корни, размерность: 2 на kol
	double **rootstable;
	rootstable = new double*[2];
	for (int ii = 0; ii<2; ii++)
		rootstable[ii] = new double[kol];

	kol = 0;
	for (int i = 1; i < n; i++)
	{
		if ((net[1][i]>1e-6 && net[1][i - 1]<1e-6) && (net[1][i - 1]<1e-6 && net[1][i + 1]>1e-6))

			//((net[1][i]>0 && net[1][i - 1]<0)
			//|| (net[1][i]<0 && net[1][i - 1]>0)
			//|| (net[1][i] == 0 && net[1][i - 1]>0 && net[0][i - 1] == a)
			//|| (net[1][i]>0 && net[1][i - 1] == 0)//ноль справа
			//|| (net[1][i] == 0 && net[1][i - 1]<0 && net[0][i - 1] == a)
			//|| (net[1][i]<0 && net[1][i - 1] == 0)//ноль справа
			//)
		{
			rootstable[0][kol] = net[0][i - 1];
			rootstable[1][kol] = net[0][i];

			kol++;
		}

		else if (net[1][i] == 0 && net[1][i - 1] == 0)
		{
			rootstable[0][kol] = net[0][i - 1];
			rootstable[1][kol] = (net[0][i] + net[0][i + 1])*0.5;
			kol++;

			rootstable[0][kol] = (net[0][i] + net[0][i + 1])*0.5;
			rootstable[1][kol] = net[0][i];


			kol++;
		}

	}

	return rootstable;
}

double BisectionMethod(double(*ff) (double), double xleft, double xright, double eps, bool ans)
{
	double curleft = xleft,
	curright = xright,
	curmiddle = (xright + xleft)*0.5;
	int cnt = 0;

	//double **net;
	//net = new double*[2];
	//for (int ii = 0; ii<2; ii++)
	//	net[ii] = new double[n];
	//double xi, yi;

	//for (int i = 0; i < n; i++)
	//{//создаем таблицу из аргументов и значений функции на отрезке с определенным выше шагом
	//	xi = a + h * i;
	//	yi = ff(xi);
	//	//cout << xi << ' ' << yi << endl;

	//	net[0][i] = xi;
	//	net[1][i] = yi;
	//}


	while (curright - curleft > 2*eps)
	{

		if ((curleft>0 && curmiddle<0)
			|| (curleft<0 && curmiddle>0)
			|| (curleft == 0 && curmiddle>0)
			|| (curleft>0 && curmiddle == 0)
			|| (curleft == 0 && curmiddle<0)
			|| (curleft<0 && curmiddle == 0)
			)
		{curright = curmiddle;}
		else if (curleft == 0 && curmiddle == 0) {curright = curmiddle;}
		
		else  {curleft = curmiddle;}

		curmiddle = (curleft + curright)*0.5;

		cnt++;
	}

	if (ans)
	{
		cout << "Количество итераций = 1. Значение корня попало в сетку" << endl;
		cout << endl;

	}
	else
	{
		cout << "Количество итераций = " << cnt << endl;
		cout << endl;
	}
	return curmiddle;
}

double NewtonMethod(double(*ff) (double), double xleft, double xright, double eps, double(*derivative) (double))
{
	double result;
	double resultprevious;
	double der;

	result = (xright + xleft)*0.5;//начальное приближение - середина отрезка
	if ((ff(xleft) - ff(xright))*1e+10 > (ff(xleft)*xright - ff(xright)*xleft))
	{
		result = (ff(xleft)*xright - ff(xright)*xleft) / (ff(xleft) - ff(xright));//начальное приближение по методу хорд 
	}

	result = 0;

	resultprevious = result + 1;

	int cnt = 0;
	while (abs(result - resultprevious) > eps)
	{
		resultprevious = result;

		if (derivative != nullptr)
		{
			result = resultprevious - ff(resultprevious) / derivative(resultprevious);
		}
		
		else 
		{
			der = (ff(resultprevious + eps) - ff(resultprevious)) / eps;
			result = resultprevious - ff(resultprevious) / der;
		
		}

		if (result<xleft || result>xright)
		{
			//возвращаем на отрезок по методу хорд
			result = (ff(resultprevious)*result - ff(result)*resultprevious) / (ff(resultprevious) - ff(result));
		}

		cnt++;
	}

	cout << endl;
	cout << "Количество итераций = " << cnt << endl;
	cout << endl;

	return result;
}

double NewtonMethodff2with8(double(*ff) (double), double xleft, double xright, double eps, double(*derivative) (double))
{
	double result;
	double resultprevious;
	double der;

	result = 8;
	resultprevious = result + 1;

	int cnt = 0;
	while (abs(result - resultprevious) > eps)
	{
		resultprevious = result;

		if (derivative != nullptr)
		{
			result = resultprevious - ff(resultprevious) / derivative(resultprevious);
		}

		else
		{
			der = (ff(resultprevious + eps) - ff(resultprevious)) / eps;
			result = resultprevious - ff(resultprevious) / der;

		}

		if (result<xleft || result>xright)
		{
			//возвращаем на отрезок по методу хорд
			//подкоренное выражение должно быть положительным
			 
			if (result >= -1)
			{
				result = (ff(resultprevious)*result - ff(result)*resultprevious) / (ff(resultprevious) - ff(result));
			}
			else
			{
				result = xleft;
				result = (ff(resultprevious)*result - ff(result)*resultprevious) / (ff(resultprevious) - ff(result));
			}
		}

		cnt++;
	}

	cout << endl;
	cout << "Количество итераций = " << cnt << endl;
	cout << endl;

	return result;
}

double* gaussmethod(double** extMatrix, int n)
{
	double *solution;
	solution = new double[n];
	double maxvalue = 0;
	int imax;

	for (int cnt = 0; cnt < n; cnt++)
	{
		solution[cnt] = 0;
	}
	/*double det = determinant(extMatrix, n);
	if (abs(det) < 1e-30)
	{
	cout << "Определитель равен 0. Не существует единственного решения." << endl;
	solution = nullptr;
	}
	else
	{*/

	for (int i = 0; i < n - 1; i++)//цикл по строкам, которые вычитаются из нижележащих
	{
		//выбор макс элемента из i-го столбца
		maxvalue = 0;
		for (int il = i; il < n; il++)
		{
			if (maxvalue < abs(extMatrix[il][i]))
			{
				maxvalue = abs(extMatrix[il][i]);
				imax = il;
			}
		}

		if (maxvalue < 1e-10)
		{
			cout << "Не существует единственного решения." << endl;
			return nullptr;
		}

		if (imax != i)
		{
			double* buf = extMatrix[imax];
			extMatrix[imax] = extMatrix[i];
			extMatrix[i] = buf;
		}

		//extMatrix[i][n] = extMatrix[i][n] / extMatrix[i][i];
		double aii = extMatrix[i][i];

		if (abs(aii) < 1e-10)
		{
			cout << "Не существует единственного решения. Последняя строка диагонализированной матрицы - нулевая" << endl;
			return nullptr;
		}

		for (int j = i; j <= n; j++)//цикл по элементам строками, которая вычитается из нижележащих  от i+1???
		{
			extMatrix[i][j] = extMatrix[i][j] / aii;
		}

		for (int ii = i + 1; ii < n; ii++)//вычитание из низлежащих строк i-ой строки
		{
			double a_ii_i = extMatrix[ii][i];
			for (int jj = i; jj <= n; jj++)
			{
				extMatrix[ii][jj] -= a_ii_i * extMatrix[i][jj];
			}
		}
	}
	//нормируем нижнюю строку
	double	 aii = extMatrix[n - 1][n - 1];
	if (abs(aii) < 1e-10)
	{
		cout << "Не существует единственного решения. Последняя строка диагонализированной матрицы - нулевая" << endl;
		return nullptr;
	}
	for (int j = n - 1; j <= n; j++)//цикл по элементам строками, которая вычитается из нижележащих  от i+1???
	{
		extMatrix[n - 1][j] = extMatrix[n - 1][j] / aii;
	}
	//printMatrix(extMatrix, n, n + 1, true);
	//обратный ход

	double sum = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (int j = i + 1; j < n; j++) //суммируем все более старшие переменные  взвешенные на коэффициенты текущей строки
		{
			sum += solution[j] * extMatrix[i][j];
		}
		solution[i] = extMatrix[i][n] - sum;//вычитаем из правой части 
	}

	//printMatrix(extMatrix, n);//печать диагонализированной (для проверки)
	return solution;
}

double** reverseMatrix(double** Matrix, int n)
{
	double **extMatrix;
	extMatrix = new double*[n];
	double **result;
	result = new double*[n];
	double *ordinary;
	ordinary = new double[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new double[n];
	for (int ii = 0; ii<n; ii++)
		extMatrix[ii] = new double[n + 1];

	for (int j = 0; j < n; j++)//заполняем правый столбец extMatrix нулями
	{
		extMatrix[j][n] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int ii = 0; ii < n; ii++)//переносим элементы из матрицы в расширенную матрицу
		{
			for (int jj = 0; jj < n; jj++)
			{
				extMatrix[ii][jj] = Matrix[ii][jj];
			}
		}
		for (int jj = 0; jj < n; jj++)//заполняем правый столбец extMatrix нулями
		{
			extMatrix[jj][n] = 0;
		}

		extMatrix[i][n] = 1;

		ordinary = gaussmethod(extMatrix, n);
		if (ordinary == nullptr)
		{
			cout << "Не существует обратной матрицы." << endl;
			return nullptr;
		}
		for (int j = 0; j < n; j++)
		{
			result[j][i] = ordinary[j];
		}

	}

	return result;
}

double* NewtonMethodSystem(double(*ff1) (double, double), double(*ff2) (double, double), double* x0, double eps,
	double(*derivative1) (double, double), double(*derivative2) (double, double),
	double(*derivative3) (double, double), double(*derivative4) (double, double),int* cntnumb=nullptr)
{
	double det = 0;
	double **Matrix;
	Matrix = new double*[2];
	double **MatrixRev;
	MatrixRev = new double*[2];
	for (int ii = 0; ii < 2; ii++)
	{
		MatrixRev[ii] = new double[2];
		Matrix[ii] = new double[2];
	}
	
	double *result;
	result = new double[2];
	double *resultprevious;
	resultprevious = new double[2];
	double *boof;
	boof = new double[2];
	double *boofnew;
	boofnew = new double[2];

	result[0] = x0[0];//присваиваем вектору результата начальное приближение
	result[1] = x0[1];

	/*result[0] = 4.0;
	result[1] = -1.0;*/

	int cnt = 0;

	while (abs(result[0] - resultprevious[0]) > eps && abs(result[1] - resultprevious[1]) > eps)
	{
		resultprevious[0] = result[0];
		resultprevious[1] = result[1];

		if (derivative1 != nullptr&& derivative2 != nullptr&& derivative3 != nullptr&& derivative4 != nullptr)
		{
			Matrix = JacobiMatrix(resultprevious, derivative1, derivative2,	derivative3, derivative4);
		}//формируем матрицу Якоби analit
		else
		{	
			Matrix[0][0] = (ff1(resultprevious[0]+eps, resultprevious[1])- ff1(resultprevious[0] , resultprevious[1]))/eps;
			Matrix[0][1] = (ff1(resultprevious[0] , resultprevious[1]+ eps) - ff1(resultprevious[0], resultprevious[1])) / eps;
			Matrix[1][0] = (ff2(resultprevious[0] + eps, resultprevious[1]) - ff2(resultprevious[0], resultprevious[1])) / eps;
			Matrix[1][1] = (ff2(resultprevious[0], resultprevious[1] + eps) - ff2(resultprevious[0], resultprevious[1])) / eps;

		}//формируем матрицу Якоби num


		/*cout << "matrix "<< endl;
		printMatrix(Matrix,2,2,true);
		cout << endl;*/

		//обращаем матрицу явно
		det = Matrix[1][1] * Matrix[0][0] - Matrix[0][1] * Matrix[1][0];
		MatrixRev[0][0] = Matrix[1][1] / det;
		MatrixRev[0][1] = -Matrix[0][1]/det;
		MatrixRev[1][0] = -Matrix[1][0]/det;
		MatrixRev[1][1] = Matrix[0][0] / det;

		/*cout << "matrixrev " << endl;
		printMatrix(MatrixRev, 2, 2, true);
		cout << endl;*/

		boofnew[0] = ff1(resultprevious[0], resultprevious[1]);
		boofnew[1] = ff2(resultprevious[0], resultprevious[1]);

		/*cout << endl;
		printVector(boofnew, 2);
		cout << endl;*/

		boof = multiplyMatrixVector(MatrixRev, boofnew, 2);

		/*cout << "boof " << endl;
		printVector(boof, 2);
		cout << endl;*/

		result = substractVector(resultprevious, boof, 2);

		cnt++;
	}

	
	if (cntnumb != nullptr) { *cntnumb = cnt; }
	else {
		cout << endl;
	cout << "Количество итераций = " << cnt << endl;
	cout << endl;
	}
	return result;
}

int main()
{
	setlocale(LC_ALL, "rus");
	double eps = 1e-3;
	//double a = -1;
	//double b = 10;
	//double a = 0;
	//double b = 1;

	double a = 0;
	double b = 1;
	int k;
	double result;

	//параметры тестов-----------------------------------
	double(*ftest) (double);
	ftest = function3;
	double(*fderivative) (double);
	fderivative = function3derivedanalytics;
	//------------------------------------------------

	cout << "Локализация корней: " << endl;
	double** test = RootLocalization(ftest, a, b, eps, k);
	printMatrix(test, 2, k, true);
	cout << endl;

	for (int i = 0; i < k; i++)
	{
		cout << "Точное решение x" << i + 1 << " лежит между " << test[0][i] << " и " << test[1][i] << endl;
		cout << endl;
	}

	cout << endl;
	cout << "Поиск точных корней уравнения методом бисекции: " << endl;
	cout << endl;

	for (int i = 0; i < k; i++)
	{
		result = BisectionMethod(ftest, test[0][i], test[1][i], eps, false);
		cout << "Точное решение x" << i + 1 << "= " << result << endl;
		cout << endl;
	}

	cout << endl;
	cout << "Поиск точных корней уравнения методом Ньютона: " << endl;
	cout << endl;

	for (int i = 0; i < k; i++)
	{
		result = NewtonMethod(ftest, test[0][i], test[1][i], eps, fderivative);//аналитическая производная
		//result = NewtonMethod(ftest, test[0][i], test[1][i], eps, nullptr);//численная производная
		//result = NewtonMethodff2with8(ftest, test[0][i], test[1][i], eps, fderivative);//2 пример, нач прибл = 8
		cout << "Точное решение x" << i + 1 << "= " << result << endl;
		cout << endl;
	}
//
//	/*cout << "Тест Якобиана: " << endl;
//	cout << endl;
//
//	double *vecttest;
//	vecttest = new double[2];
//	vecttest[0] = 1;
//	vecttest[1] = 2;
//	double** resulttestj = JacobiMatrix(vecttest, system1function1derivedanalyticsx, 
//		system1function1derivedanalyticsy, system1function2derivedanalyticsx, system1function1derivedanalyticsy);
//	cout << "J = " << endl;
//	printMatrix(resulttestj, 2, 2, true);
//	cout << endl;

	cout << "Решение системы нелинейных уравнений (1) методом Ньютона: " << endl;
	//решения (4,-1) (-4,1) 
	double *x0test;
	x0test = new double[2];
	x0test[0] = 1;
	x0test[1] = 3;
	cout << "Начальное приближение имеет вид:  " << endl;
	printVector(x0test, 2);

	 //double *resulttest = NewtonMethodSystem(system1function1, system1function2, x0test, eps, system1function1derivedanalyticsx,
	 //system1function1derivedanalyticsy, system1function2derivedanalyticsx, system1function2derivedanalyticsy);
	double *resulttest = NewtonMethodSystem(system1function1, system1function2, x0test, eps, nullptr,
		nullptr, nullptr, nullptr);

	cout << "Ответ:  " << endl;
	printVector(resulttest, 2);
	cout << endl;

	cout << "Решение системы нелинейных уравнений (3) методом Ньютона: " << endl;
	//решения (-3,1) (1,-3) (1,2) (2,1)
	double *x0test1;
	x0test1 = new double[2];
	x0test1[0] = -4;
	x0test1[1] = 2;
	cout << "Начальное приближение имеет вид:  " << endl;
	printVector(x0test1, 2);

	//double *resulttest1 = NewtonMethodSystem(system2function1, system2function2, x0test1, eps, system2function1derivedanalyticsx,
		//system2function1derivedanalyticsy, system2function2derivedanalyticsx, system2function2derivedanalyticsy);
	double *resulttest1 = NewtonMethodSystem(system3function1, system3function2, x0test1, eps, nullptr,
		nullptr, nullptr, nullptr);

	cout << "Ответ:  " << endl;
	printVector(resulttest1, 2);
	cout << endl;

	//открываем файл для построения области сходимости системы |x|<=10, |y|<=10
	int testcnt = 0;
	double *x0testiter;
	x0testiter = new double[2];
	double *resulttestiter=nullptr;
	ofstream fout;
	//fout.open("system3analytics.txt"); // связываем объект с файлом
	fout.open("system3numerical.txt"); // связываем объект с файлом
	
	for (int i = -10; i < 11; i++)
	{

		for (int j = -10; j < 11; j++)
		{
			x0testiter[0] = i;
			x0testiter[1] = j;
			
			//система 2
			 //resulttestiter = NewtonMethodSystem(system3function1, system3function2, x0testiter, eps, system3function1derivedanalyticsx,
			//system3function1derivedanalyticsy, system3function2derivedanalyticsx, system3function2derivedanalyticsy, &testcnt);
			resulttestiter = NewtonMethodSystem(system3function1, system3function2, x0testiter, eps, nullptr,
				nullptr, nullptr, nullptr, &testcnt);
			
				 //система 1
			 /*resulttestiter = NewtonMethodSystem(system1function1, system1function2, x0testiter, eps, system1function1derivedanalyticsx,
			 system1function1derivedanalyticsy, system1function2derivedanalyticsx, system1function2derivedanalyticsy, &testcnt);*/
			/* resulttestiter = NewtonMethodSystem(system1function1, system1function2, x0testiter, eps, nullptr,
				 nullptr, nullptr, nullptr, &testcnt);*/


			fout << i << " " << j << " " << testcnt << endl; // запись строки в файл
		}
	}
	fout.close(); // закрываем файл


	cout << endl;
	system("pause");
	return 0;
}

