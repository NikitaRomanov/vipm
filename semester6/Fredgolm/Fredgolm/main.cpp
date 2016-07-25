#include <iostream>
#include <iomanip>
#include <vector>
#include <unordered_map>

using namespace std;

//нули полинома Лежандра
const static vector<double> lejRoots = { -sqrt((35 + 2 * sqrt(70.0)) / 63.0), -sqrt((35 - 2 * sqrt(70.0)) / 63.0), 0, sqrt((35 - 2 * sqrt(70.0)) / 63.0),sqrt((35 + 2 * sqrt(70.0)) / 63.0) };
//веса формулы Гаусса
const static vector<double> weights = { (322 - 13 * sqrt(70.0)) / 900.0, (322 + 13 * sqrt(70.0)) / 900.0, 128 / 225.0, (322 + 13 * sqrt(70.0)) / 900.0, (322 - 13 * sqrt(70.0)) / 900.0 };
//число точек в формуле Гаусса
const static size_t r = weights.size();
const double pi = 3.14159265358979323846;

//инетгрирование по Гауссу. напрямую здесь не используется
template <class _Function>
double integrate(double a, double b, const _Function& f) {
	double sum = 0;
	for (size_t i = 0; i < 5; i++)
	{
		double t = (a + b) / 2.0 + (b - a)*lejRoots[i] / 2.0;
		t = f(t) * weights[i];
		sum += t;
	}
	sum *= (b - a) / 2.0;
	return sum;
}

//полином Лагранжа. напрямую здесь не используется
template <class _Function>
double lagr(double p, double a, double b, size_t N, const _Function& f) {
	double sum = 0;
	for (size_t i = 0; i < N; i++)
	{
		double ti = a + i*(b - a) / N;
		double mult = 1;
		for (size_t j = 0; j < N; j++)
		{
			double tj = a + j*(b - a) / N;
			if (i != j)
				mult *= (p - tj) / (ti - tj);
		}
		mult *= f(ti);
		//wcout << x.at(tkj[k][i]) << endl;
		sum += mult;
	}
	return sum;
}

//правая часть
double f(double x) {
	return 2 * x;
}
//ядро 
double h(double x, double s) {
	return x*s;
}
//точное решение
double exact(double x) {
	return 3 * x;
}

int main() {
	wcout << "Hi!" << endl;
#ifdef _DEBUG
	wcout << "Select Release instead of Debug. The program will run much faster." << endl;
#endif

	//левая граница
	double a = 0;
	//правая граница
	double b = 1;
	//лямбда
	double lmb = 1;
	//число точек разбиения
	size_t N = 10;
	//число итераций
	size_t iters_count = 40;

	vector<double> t;
	for (size_t i = 0; i <= N; i++)
	{
		double ti = a + (b - a)*i / N;
		t.push_back(ti);
	}
	vector<vector<double>> tkj(t.size());
	tkj[0].assign(r + 2, 0);
	for (size_t i = 1; i < t.size(); i++)
	{
		tkj[i].assign(r + 2, 0);
		tkj[i][0] = t[i - 1];
		for (size_t j = 1; j <= r; j++)
		{
			tkj[i][j] = (t[i - 1] + t[i]) / 2 + (t[i] - t[i - 1])*lejRoots[j - 1] / 2;
		}
		tkj[i][r + 1] = t[i];
	}

	unordered_map<double, double> x;
	for (size_t i = 1; i < t.size(); i++)
	{
		for (size_t j = 0; j <= r + 1; j++)
		{
			x[tkj[i][j]] = f(tkj[i][j]);
		}
	}

	for (size_t iterations = 0; iterations < iters_count; iterations++)
	{
		unordered_map<double, double> x_temp(x);
		for (size_t i = 1; i < t.size(); i++)
		{
			for (size_t m = 0; m <= r + 1; m++)
			{
				double p = tkj[i][m];
				double sum = 0;
				for (size_t k = 1; k < t.size(); k++)
				{
					double part_sum = 0;
					for (size_t j = 1; j <= r; j++)
					{
						part_sum += weights[j - 1] * h(p, tkj[k][j])*x.at(tkj[k][j]);
					}
					part_sum *= (t[k] - t[k - 1]) / 2.0;
					sum += part_sum;
				}
				sum *= lmb;
				sum += f(p);
				x_temp.at(p) = sum;
			}
		}
		x = move(x_temp);
	}

	auto spline = [&](double p) {
		size_t k = 0;
		for (size_t i = 1; i < t.size(); i++) {
			if (p < t[i]) {
				k = i;
				break;
			}
		}

		double sum = 0;
		for (size_t i = 0; i <= r + 1; i++)
		{
			double mult = 1;
			for (size_t j = 0; j <= r + 1; j++)
			{
				if (i != j)
					mult *= (p - tkj[k][j]) / (tkj[k][i] - tkj[k][j]);
			}
			mult *= x.at(tkj[k][i]);
			//wcout << x.at(tkj[k][i]) << endl;
			sum += mult;
		}
		return sum;
	};

	wios state(nullptr);
	state.copyfmt(wcout);
	//точность 8 знаков при выводе на экран
	wcout << setprecision(8) << fixed;
	double eps1 = abs(x.at(tkj[1][0]) - exact(tkj[1][0]));
	for (size_t i = 1; i < t.size(); i++)
	{
		for (size_t j = 0; j <= r; j++)
		{
			double diff = abs(x.at(tkj[i][j]) - exact(tkj[i][j]));
			//wcout << tkj[i][j] << setw(16) << x.at(tkj[i][j]) << setw(16) << exact(tkj[i][j]) << setw(20) << scientific << diff << fixed << endl;
			if (diff > eps1)
				eps1 = diff;
		}
	}
	wcout << L"Max error at points: " << endl;
	wcout << scientific << eps1 << fixed << endl;
	wcout << L"==========" << endl;

	//величина шага при поиске погрешности
	double step = 0.0001;
	double p = a;
	double eps2 = abs(spline(p) - exact(p));
	while (p <= b) {
		double diff = abs(spline(p) - exact(p));
		//wcout << p << setw(16) << spline(p) << setw(16) << exact(p) << setw(20) << scientific << spline(p) - exact(p) << fixed << endl;
		if (diff > eps2)
			eps2 = diff;
		p += step;
	}
	wcout << L"Max error on the interval (spline): " << endl;
	wcout << scientific << eps2 << fixed << endl;
	wcout.copyfmt(state);

	system("pause");
}
