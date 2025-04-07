#include"user_funs.h"
#include"solution.h"
#include <cmath>

#include <fstream>
extern std::ofstream file;

#define M_PI 3.14

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) { //zmiana objetosci wody w zbiorniku
    double fOutA = 0, fOutB = 0;
    double fInB = 0.01, tInB = 10.0;
    double tA = 90.0;
    matrix dY(3, 1);
    if (Y(0) > 0)fOutA = 0.98 * 0.63 * m2d(ud2) * sqrt(2.0 * 9.81 * Y(0) / 0.5);
    if (Y(1) > 0)fOutB = 0.98 * 0.63 * 0.00365665 * sqrt(2.0 * 9.81 * Y(1) / 1.0);
    dY(0) = -fOutA;
    dY(1) = (fOutA + fInB - fOutB);
    dY(2) = (fInB / Y(1)) * (tInB - Y(2)) + (fOutA / Y(1)) * (tA - Y(2));
    return dY;
}
matrix ff1R(matrix x, matrix ud1, matrix ud2) { // zmiana temperatury w zbiorniku
    matrix y;
    matrix Y0 = matrix(3, new double[3] { 5, 1, 20 });
    matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);
    int n = get_len(Y[0]);
    double max = Y[1](0, 2);
    for (int i = 0; i < n; i++) {
        if (max < Y[1](i, 2)) max = Y[1](i, 2);
    }
    y = abs(max - 50);
    return y;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    double arg = m2d(x);
    double result = -cos(0.1 * arg) * exp(-pow(0.1 * arg - 2 * 3.14, 2)) + 0.002 * pow(0.1 * arg, 2);
    matrix res = matrix(result);
    return res;
}

//matrix ff3T(matrix x1, matrix x2, matrix ud2) {
//    return (sin(M_PI * sqrt(pow(x1() / M_PI, 2) + pow(x2() / M_PI, 2)))) / (M_PI * sqrt(pow(x1() / M_PI, 2) + pow(x2() / M_PI, 2)));
//}

matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
	matrix Y;

	double x_1 = m2d(x(0));
	double x_2 = m2d(x(1));
	double numerator = sin(M_PI * sqrt(pow((x_1 / M_PI), 2) + pow((x_2 / M_PI), 2)));
	double denominator = M_PI * sqrt(pow((x_1 / M_PI), 2) + pow((x_2 / M_PI), 2));
	Y = numerator / denominator;

	double a = ud1(0);
	double c = ud2(0); //c>0

	//ograniczenia funckja kary wewnêtrzna
	if (1 > x_1) {
		Y = 1e05;
	}
	else {
		Y = Y - c / (1 - x(0));
	}

	if (1 > x_2) {
		Y = 1e05;
	}
	else {
		Y = Y - c / (1 - x_2);
	}

	if (norm(x) > a) {
		Y = 1e05;
	}
	else {
		Y = Y - c / (norm(x) - a);
	}

	return Y;

	//ograniczenia fucnkja kary zewnêtrzna 
	if (ud2(1) > 1) {
		//g1
		if (1 > x_1) {
			Y = Y + c * pow((-x_1) + 1, 2);
		}

		//g2
		if (1 > x_2) {
			Y = Y + c * pow((-x_2) + 1, 2);
		}

		//g3
		if (norm(x) > a) {
			Y = Y + c * pow(norm(x) - a, 2);
		}

		return Y;
	}
}

//matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
//	double m = 0.6;  // masa [kg]
//	double r = 0.12; // promieñ pi³ki [m]
//
//	double g = 9.81; // si³a grawitacji [m/s^2]
//	double C = 0.47; // wspolczynnik oporu
//	double ro = 1.2; // gestosc powietrza [kg/m^3] 
//
//	double S = (M_PI * pow(r, 2)); // najwiêkszy przekrój bry³y (pole ko³a) [m^2]
//	double V_x = m2d(Y(1));        // prêdkoœæ wzglêdem osi x [m/s]
//	double V_y = m2d(Y(3));        // prêdkoœæ wglêdem osi y [m/s]
//	double omega = m2d(ud2);       // predkosc k¹towa [rad/s]
//
//	//si³y oporu powietrza 
//	double D_x = (0.5 * C * ro * S * V_x * abs(V_x));
//	double D_y = (0.5 * C * ro * S * V_y * abs(V_y));
//
//	//si³y Magnusa
//	double F_Mag_x = (ro * V_y * omega * M_PI * pow(r, 3));
//	double F_Mag_y = (ro * V_x * omega * M_PI * pow(r, 3));
//
//	//równania ruchu pi³ki (pochodne cz¹stkowe drugiego stopnia po czasie)
//	matrix dY(4, 1);
//	//dx
//	dY(0) = V_x;
//	dY(1) = (((-D_x) - F_Mag_x) / m);
//	//dy
//	dY(2) = V_y;
//	dY(3) = (((-m) * g - D_y - F_Mag_y) / m);
//
//	return dY;
//}
//
//matrix ff3R(matrix x, matrix ud1, matrix ud2) {
//	matrix y;
//	matrix Y0(4, new double[4] { 0, x(0), 100, 0 });
//	matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, ud1, x(1));
//	int N = get_len(Y[0]);
//
//	int closestIndex50 = 1e10;
//	int closestIndex0 = 1e10;
//	double minHeight50 = 1e10;
//	double minHeight = 1e10;
//
//	for (int i = 0; i < N - 1; ++i) {
//		double y_i = Y[1](i, 2);
//		if (abs(y_i) < minHeight) { // Szukanie takiego x, dla ktorego y = 0
//			minHeight = abs(y_i);
//			closestIndex0 = i;
//		}
//		if (abs(y_i - 50) < minHeight50) { // Szukanie takiego x, dla ktorego y = 50
//			minHeight50 = abs(y_i - 50);
//			closestIndex50 = i;
//		}
//	}
//
//	if (minHeight50 == 1e10 && minHeight == 1e10)
//		return matrix(1e10);
//
//	y = -Y[1](closestIndex0, 0); // Znalezienie minimum
//
//	if (abs(x(0)) > 10)
//		y = y + ud2() * pow(abs(x(0)) - 10, 2); //Kary dla x, omega, 50
//	if (abs(x(1)) > 15)
//		y = y + ud2() * pow(abs(x(1)) - 15, 2);
//	if (abs(Y[1](closestIndex50, 0) - 5) > 1)
//		y = y + ud2() * pow(abs(Y[1](closestIndex50, 0) - 5) - 0.5, 2);
//
//	ofstream file("symulacja.csv"); // Zapis do csv
//	if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
//	file << "t,x,y" << endl;
//
//
//	for (int i = 0; i < N; ++i) {
//		file << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 2) << endl; \
//	}
//
//	//file.close();
//
//	return y;
//}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
	double C = 0.47;
	double r = 0.12;
	double m = 0.6;
	double ro = 1.2;
	double g = 9.81;

	double S = M_PI * pow(r, 2);

	double Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1));
	double Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3));

	double Fmx = M_PI * ro * Y(3) * m2d(ud2) * pow(r, 3);
	double Fmy = M_PI * ro * Y(1) * m2d(ud2) * pow(r, 3);

	matrix dY(4, 1);
	dY(0) = Y(1);                        //x0 = v0x
	dY(1) = (-Dx - Fmx) / m;            //V0x = ax
	dY(2) = Y(3);                        //y0 = V0y
	dY(3) = (-m * g - Dy - Fmy) / m;    //V0y = ay

	return dY;

}

void print_csv(matrix*& Y, int N) {
	ofstream file("symulacja.csv");
	if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
	file << "t,x,y" << endl;

	for (int i = 0; i < N; ++i)
		file << Y[0](i) << "," << Y[1](i, 0) << "," << Y[1](i, 2) << endl;

	file.close();
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0(4, new double[4] { 0, x(0), 100, 0 });
	matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, ud1, x(1));
	int N = get_len(Y[0]);

	int minHeight50Index = 1e10, minHeight0Index = 1e10;
	double minHeight50 = 1e10;
	double minHeight0 = 1e10;

	for (int i = 0; i < N - 1; ++i) {
		double y_i = Y[1](i, 2);
		if (abs(y_i - 50) < minHeight50) { //Szukanie x, takiego, ze y = 0
			minHeight50 = abs(y_i - 50);
			minHeight50Index = i;
		}
		if (abs(y_i) < minHeight0) { //Szukanie x, takeigo, ze y = 50
			minHeight0 = abs(y_i);
			minHeight0Index = i;
		}
	}

	if (minHeight50Index == 1e10 || minHeight0Index == 1e10)
		return matrix(1e10);

	y = -Y[1](minHeight0Index, 0); //Szukanie minimum

	if (abs(x(0)) - 10 > 0) { //Kary dla y = 0, omega i y =0
		y = y + ud2() * pow(abs(x(0)) - 10, 2);
	}
	if (abs(x(1)) - 15 > 0) {
		y = y + ud2() * pow(abs(x(1)) - 15, 2);
	}
	if (abs(Y[1](minHeight50Index, 0) - 5) - 0.5 > 0) {
		y = y + ud2() * pow(abs(Y[1](minHeight50Index, 0) - 5) - 0.5, 2);
	}

	print_csv(Y, N);

	return y;
}

matrix gradient(matrix x, matrix ud1, matrix ud2) {
	matrix G(2, 1);

	G(0) = x(0) * 10 + 8 * x(1) - 34; //pochodna cz¹stkowa wzglêdem x1
	G(1) = x(0) * 8 + 10 * x(1) - 38; //pochodna cz¹stkowa wzglêdem x2

	return G;
}

matrix hesjan(matrix x, matrix ud1, matrix ud2) {
	matrix H(2, 2); //macierz pochodnych cz¹stkowych drugiego rzêdu

	H(0, 0) = 10;
	H(1, 0) = 8;
	H(0, 1) = 8;
	H(1, 1) = 10;

	return H;
}

matrix ffT4(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	if (isnan(ud2(0, 0))) {
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	}
	else {
		matrix trans_x = x;

		while (!isnan(ud2(0, 0))) {
			trans_x = ud2[0] + trans_x * ud2[1];

			if (isnan(ud1(0, 0))) {
				break;
			}
			ud2 = ud1;
		}
		y = pow(trans_x(0) + 2 * trans_x(1) - 7, 2) + pow(2 * trans_x(0) + trans_x(1) - 5, 2);
	}
	return y;
}

matrix ff4R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	int m = 100;
	int n = get_len(x);
	static matrix X(n, m), Y(1, m);
	if (solution::f_calls == 1) {
		ifstream in("XData.txt");
		in >> X;
		in.close();
		in.open("YData.txt");
		in >> Y;
		in.close();
	}
	int P = 0;
	double h;
	y = 0;
	for (int i = 0; i < m; i++) {
		h = m2d(trans(x) * X[i]);
		h = 1.0 / (1.0 + exp(-h));
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
	}
	y = y / m;
	return y;
}

matrix gf(matrix x, matrix ud1, matrix ud2) {
	int m = 100;
	int n = get_len(x);
	matrix g(n, 1);
	static matrix X(n, m), Y(1, m);
	if (solution::g_calls == 1) {
		ifstream in("XData.txt");
		in >> X;
		in.close();
		in.open("YData.txt");
		in >> Y;
		in.close();
	}

	double h;
	for (int j = 0; j < n; ++j) {
		for (int i = 0; i < m; ++i) {
			h = m2d(trans(x) * X[i]);
			h = 1 / (1 + exp(-h));
			g(j) = g(j) + X(j, i) * (h - Y(0, i));
		}
		g(j) = g(j) / m;
	}
	return g;
}

matrix ff5T(matrix x, matrix ud1, matrix ud2) {

	//cout << ud1 << endl << ud2 << endl;
	//cout << x(0) << endl;
	//cout << x(1) << endl;
	//cout << ud1 * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2)) << endl;
	if (isnan(ud2(0, 0))) {
		matrix f1;
		matrix f2;

		f1 = ud1(0) * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
		f2 = (1.0 / ud1(0)) * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));

		matrix y;

		y = ud1(1) * f1 + (1 - ud1(1)) * f2;

		return y;

	}
	else {
		matrix y;

		y = ff5T(ud2[0] + x * ud2[1], ud1);

		return y;

	}

}

matrix ff5T_f1(matrix x, matrix ud1, matrix ud2) {
	matrix f1 = ud1(0) * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
	return f1;
}

matrix ff5T_f2(matrix x, matrix ud1, matrix ud2) {
	matrix f2 = (1.0 / ud1(0)) * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
	return f2;
}