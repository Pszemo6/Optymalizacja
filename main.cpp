/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

//nie trzrba w sproawizdaniu opisaywaæ dzia³ania fucnkji- to na egzamin
//na koniec trzeba wrzucic kod trzech fucnkji i funckji lab_1 i funcknji którue uzycwam te ff chyba)
//spraozdanie wysyla jedna osoba
//nazwa pliku od naszych nazwisk
//do czwartku o 8:00


#include"opt_alg.h"
//#include <fstream>
#include <cstdlib>
#include <ctime>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab5();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{

	// Otwieranie pliku CSV do zapisu wyników
	ofstream file ("lab1_tab2_fib.csv");


	// Sprawdzenie, czy plik zosta³ poprawnie otwarty
	if (!file.is_open ()) {
		cout << "file error" << endl;
		return;
	}

	// Nag³ówki w pliku CSV, dostosowane do struktury tabeli w Excelu
	file << "Da;y;T\n";

	solution profile = fib(ff1R, 0.0001, 0.01, 0.0000001);
	matrix YO = matrix(3, new double[3] {5, 1, 20});
	matrix* Y = solve_ode(df1, 0, 1, 1000, YO, NAN, profile.x(0));

	cout << profile.x(0) << endl << profile.y(0) << endl << profile.f_calls << endl;

	cout << Y[1] << endl;



	file << Y[1] << "\n";

	// Zamkniêcie pliku po zakoñczeniu
	file.close ();



	// Otwieranie pliku CSV do zapisu wyników
//	ofstream file ("lab1_tab1.csv");

	// Generator liczb losowych
//	srand (time (NULL));
	/*random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> distribution(-100.0, 100.0);*/

	// Sprawdzenie, czy plik zosta³ poprawnie otwarty
	/*if (!file.is_open ()) {
		cout << "file error" << endl;
		return;
	}*/

	// Nag³ówki w pliku CSV, dostosowane do struktury tabeli w Excelu
//file << "Wspó³czynnik ekspansji,Lp.,x(0),a,b,Liczba wywo³añ funkcji celu,"
//	<< "Fib_x*,Fib_y*,Fib_Liczba wywo³añ funkcji celu,Fib_Minimum,"
//	<< "Lag_x*,Lag_y*,Lag_Liczba wywo³añ funkcji celu,Lag_Minimum\n";

	//double wspolczynnik_ekspansji = 50; // Przyk³ad wartoœci wspó³czynnika ekspansji, mo¿na dostosowaæ

	// Pêtla dla 100 optymalizacji
	//for (int i = 0; i < 1; i++) {
		// Losowanie punktu startowego
	//	double startPoint = rand () % 201 - 100;

		// Ekspansja wstêpna, zawê¿enie przedzia³u
	//	double* part = expansion (ff1T, startPoint, 0.5, wspolczynnik_ekspansji, 2000);

		// Wynik metody Fibonacciego
//		solution minimumFib = fib (ff1T, -100, 100, 0.0000001);

		// Wynik metody Lagrange'a
		//solution minimumLag = lag (ff1T, -100, 100, 0.00001, 0.0000001, 1000);

		// Okreœlenie, czy minimum jest lokalne czy globalne (Fibonacciego)
	//	string fib_minimum_type = (minimumFib.x (0) < 50) ? "lokalne" : "globalne";

		// Okreœlenie, czy minimum jest lokalne czy globalne (Lagrange'a)
	//	string lag_minimum_type = (minimumLag.x (0) < 50) ? "lokalne" : "globalne";

		// Zapis wyników do pliku CSV zgodnie z szablonem
	//	file << wspolczynnik_ekspansji << ","   // Wspó³czynnik ekspansji
	//		<< i + 1 << ","                    // Lp.
	//		<< startPoint << ","               // Punkt startowy x(0)
	//		<< part[0] << ","                  // Dolny zakres przedzia³u (a)
	//		<< part[1] << ","                  // Górny zakres przedzia³u (b)
	//		<< solution::f_calls << ","        // Liczba wywo³añ funkcji celu dla ekspansji
	//		<< minimumFib.x (0) << ","          // Wynik optymalizacji Fibonacciego (x*)
	//		<< minimumFib.y (0) << ","          // Wartoœæ funkcji celu dla Fibonacciego (y*)
	//		<< minimumFib.f_calls << ","       // Liczba wywo³añ funkcji celu dla Fibonacciego
	//		<< fib_minimum_type << ","         // Typ minimum (lokalne/globalne) dla Fibonacciego
	//		<< minimumLag.x (0) << ","          // Wynik optymalizacji Lagrange'a (x*)
	//		<< minimumLag.y (0) << ","          // Wartoœæ funkcji celu dla Lagrange'a (y*)
	//		<< minimumLag.f_calls << ","       // Liczba wywo³añ funkcji celu dla Lagrange'a
	//		<< lag_minimum_type << "\n";       // Typ minimum (lokalne/globalne) dla Lagrange'a
	//}

	//// Zamkniêcie pliku po zakoñczeniu
	//file.close ();

	//}
}

void lab2()
{

}

void lab3()
{
	//// Sekcja Tostowa

	//ofstream file("lab_03_tostowe.csv");
	//if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
	//file << "x0_1, x0_2, xz_1, xz_2, norm_z, y_z, f_calls_z"
	//	<< "xw_1, xw_2, norm_w, y_w, f_calls_w\n";

	//matrix x0 = matrix(2, 1, 1.0);

	//double c0 = 2;
	//matrix a[3] = { 4, 4.4934, 5 };

	//const double epsilon = 1e-3;
	//const int Nmax = 10000;

	//double x0_1[100] = { 0 };
	//double x0_2[100] = { 0 };

	////a[0]

	//for (int i = 0; i < 100; i++) {
	//	do
	//		x0 = 5 * rand_mat(2, 1) + 1;
	//	while (norm(x0) > a[0]);

	//	x0_1[i] = x0(0);
	//	x0_2[i] = x0(1);
	//}

	//for (int i = 0; i < 100; i++) {
	//	x0(0) = x0_1[i];
	//	x0(1) = x0_2[i];

	//	cout << x0(0) << ";" << x0(1) << ";";

	//	solution curry_zewnetrzne1 = pen(ff3T, x0, c0, 2, epsilon, Nmax, a[0]);
	//	cout << curry_zewnetrzne1.x(0) << ";" << curry_zewnetrzne1.x(1) << ";" << norm(curry_zewnetrzne1.x) << ";" << curry_zewnetrzne1.y[0] << solution::f_calls << ";";

	//	file
	//		<< x0(0) << ","
	//		<< x0(1) << ","
	//		<< curry_zewnetrzne1.x(0) << ","
	//		<< curry_zewnetrzne1.x(1) << ","
	//		<< norm(curry_zewnetrzne1.x) << ","
	//		<< curry_zewnetrzne1.y(0) << ","
	//		<< solution::f_calls << ",";

	//	solution curry_wewnetrzne1 = pen(ff3T, x0, c0, 0.5, epsilon, Nmax, a[0]);
	//	cout << curry_wewnetrzne1.x(0) << ";" << curry_wewnetrzne1.x(1) << ";" << norm(curry_wewnetrzne1.x) << ";" << curry_wewnetrzne1.y[0] << solution::f_calls << endl;
	//	file
	//		<< curry_wewnetrzne1.x(0) << ","
	//		<< curry_wewnetrzne1.x(1) << ","
	//		<< norm(curry_wewnetrzne1.x) << ","
	//		<< curry_wewnetrzne1.y(0) << ","
	//		<< solution::f_calls << ","
	//		<< "\n";

	//}

	//// a[1]

	//for (int i = 0; i < 100; i++) {
	//	do
	//		x0 = 5 * rand_mat(2, 1) + 1;
	//	while (norm(x0) > a[0]);

	//	x0_1[i] = x0(0);
	//	x0_2[i] = x0(1);
	//}

	//for (int i = 0; i < 100; i++) {
	//	x0(0) = x0_1[i];
	//	x0(1) = x0_2[i];

	//	cout << x0(0) << ";" << x0(1) << ";";

	//	solution curry_zewnetrzne2 = pen(ff3T, x0, c0, 2, epsilon, Nmax, a[1]);
	//	cout << curry_zewnetrzne2.x(0) << ";" << curry_zewnetrzne2.x(1) << ";" << norm(curry_zewnetrzne2.x) << ";" << curry_zewnetrzne2.y[0] << solution::f_calls << ";";

	//	file
	//		<< x0(0) << ","
	//		<< x0(1) << ","
	//		<< curry_zewnetrzne2.x(0) << ","
	//		<< curry_zewnetrzne2.x(1) << ","
	//		<< norm(curry_zewnetrzne2.x) << ","
	//		<< curry_zewnetrzne2.y(0) << ","
	//		<< solution::f_calls << ",";

	//	solution curry_wewnetrzne2 = pen(ff3T, x0, c0, 0.5, epsilon, Nmax, a[1]);
	//	cout << curry_wewnetrzne2.x(0) << ";" << curry_wewnetrzne2.x(1) << ";" << norm(curry_wewnetrzne2.x) << ";" << curry_wewnetrzne2.y[0] << solution::f_calls << endl;
	//	file
	//		<< curry_wewnetrzne2.x(0) << ","
	//		<< curry_wewnetrzne2.x(1) << ","
	//		<< norm(curry_wewnetrzne2.x) << ","
	//		<< curry_wewnetrzne2.y(0) << ","
	//		<< solution::f_calls << ","
	//		<< "\n";

	//}

	//// a[2]

	//for (int i = 0; i < 100; i++) {
	//	do
	//		x0 = 5 * rand_mat(2, 1) + 1;
	//	while (norm(x0) > a[0]);

	//	x0_1[i] = x0(0);
	//	x0_2[i] = x0(1);
	//}

	//for (int i = 0; i < 100; i++) {
	//	x0(0) = x0_1[i];
	//	x0(1) = x0_2[i];

	//	cout << x0(0) << ";" << x0(1) << ";";

	//	solution curry_zewnetrzne3 = pen(ff3T, x0, c0, 2, epsilon, Nmax, a[2]);
	//	cout << curry_zewnetrzne3.x(0) << ";" << curry_zewnetrzne3.x(1) << ";" << norm(curry_zewnetrzne3.x) << ";" << curry_zewnetrzne3.y[0] << solution::f_calls << ";";

	//	file
	//		<< x0(0) << ","
	//		<< x0(1) << ","
	//		<< curry_zewnetrzne3.x(0) << ","
	//		<< curry_zewnetrzne3.x(1) << ","
	//		<< norm(curry_zewnetrzne3.x) << ","
	//		<< curry_zewnetrzne3.y(0) << ","
	//		<< solution::f_calls << ",";

	//	solution curry_wewnetrzne3 = pen(ff3T, x0, c0, 0.5, epsilon, Nmax, a[2]);
	//	cout << curry_wewnetrzne3.x(0) << ";" << curry_wewnetrzne3.x(1) << ";" << norm(curry_wewnetrzne3.x) << ";" << curry_wewnetrzne3.y[0] << solution::f_calls << endl;
	//	file
	//		<< curry_wewnetrzne3.x(0) << ","
	//		<< curry_wewnetrzne3.x(1) << ","
	//		<< norm(curry_wewnetrzne3.x) << ","
	//		<< curry_wewnetrzne3.y(0) << ","
	//		<< solution::f_calls << ","
	//		<< "\n";

	//}


	//file.close();

	//Problem rzeczywisty 
	matrix x0 = matrix(2, 1);
	double c0 = 2;
	matrix a[3] = { 4, 4.4934, 5 };

	const double epsilon = 1e-5;
	const int Nmax = 10000;

	x0(0) = 7;        //V0x 
	x0(1) = 7;    //Omega [rad/s]
	solution wynik_real = pen(ff3R, x0, c0, 2, epsilon, Nmax);
	cout << wynik_real;
}

void lab4()
{
	////SEKCJA TOSTOWA

	//ofstream file("lab_04_tostowe.csv");
	//if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
	//file << "x0_1, x0_2, x1_sd, x2_sd, y_sd, f_calls_sd, g_calls_sd, "
	//	<< "x1_cg, x2_cg, y_cg, f_calls_cg, g_calls_cg, "
	//	<< "x1_newton, x2_newton, y_newton, f_calls_newton, g_calls_newton, H_calls_newton\n";

	//matrix x0 = matrix(2, 1, 0.0);
	//const double epsilon = 1e-3;
	//const int Nmax = 10000;
	//matrix H[3] = { 0.05, 0.12, 2137 };

	//double point1[100] = { 0 };
	//double point2[100] = { 0 };

	//for (int i = 0; i < 100; i++) {
	//	x0 = 20 * rand_mat(2, 1) - 10;
	//	point1[i] = x0(0);
	//	point2[i] = x0(1);
	//}

	//for (int j = 0; j < 3; j++) {
	//	double h = m2d(H[j]);
	//	cout << h << endl;
	//	for (int i = 0; i < 100; i++) {
	//		x0(0) = point1[i];
	//		x0(1) = point2[i];

	//		file
	//			<< x0(0) << ","
	//			<< x0(1) << ",";

	//		solution sd = SD(ff4T, gradient, x0, h, epsilon, Nmax);
	//		//cout << "SD:" << endl << sd << endl;

	//		file
	//			<< sd.x(0) << ","
	//			<< sd.x(1) << ","
	//			<< sd.y(0) << ","
	//			<< solution::f_calls << ","
	//			<< solution::g_calls << ",";

	//		solution cg = CG(ff4T, gradient, x0, h, epsilon, Nmax);
	//		//cout << "CG:" << endl << cg << endl;

	//		file
	//			<< cg.x(0) << ","
	//			<< cg.x(1) << ","
	//			<< cg.y(0) << ","
	//			<< solution::f_calls << ","
	//			<< solution::g_calls << ",";

	//		solution newton = Newton(ff4T, gradient, hesjan, x0, h, epsilon, Nmax);
	//		//cout << "Newton:" << endl << newton << endl;

	//		file
	//			<< newton.x(0) << ","
	//			<< newton.x(1) << ","
	//			<< newton.y(0) << ","
	//			<< solution::f_calls << ","
	//			<< solution::g_calls << ","
	//			<< solution::H_calls << ","
	//			<< "\n";
	//	}
	//}

	//file.close();

	////SEKCJA WYKRESOWA
	//matrix x0 = matrix(2, 1, 0.0);
	//const double epsilon = 1e-3;
	//const int Nmax = 10000;
	//x0(0) = 0.5;
	//x0(1) = 1;
	
	//solution sd1 = SD(ff4T, gradient, x0, 0.05, epsilon, Nmax);
	//solution sd2 = SD(ff4T, gradient, x0, 0.12, epsilon, Nmax);
	//solution sd3 = SD(ff4T, gradient, x0, 2137, epsilon, Nmax);
	//solution cg1 = CG(ff4T, gradient, x0, 0.05, epsilon, Nmax);
	//solution cg2 = CG(ff4T, gradient, x0, 0.12, epsilon, Nmax);
	//solution cg3 = CG(ff4T, gradient, x0, 2137, epsilon, Nmax);
	//solution newton1 = Newton(ff4T, gradient, hesjan, x0, 0.05, epsilon, Nmax);
	//solution newton2 = Newton(ff4T, gradient, hesjan, x0, 0.12, epsilon, Nmax);
	//solution newton3 = Newton(ff4T, gradient, hesjan, x0, 2137, epsilon, Nmax);


	////SEKCJA REALNA
	//matrix x1(3, 1, 0.0); // wyliczone z wczeœniejszych

	//double h = 0.01;
	////double h = 0.001;
	////double h = 0.0001;
	//int Nmax = 10000;

	//solution Real_sol = SD(ff4R, gf, x1, h, 0.01, Nmax);

	//int m = 100;
	//static matrix X(3, m), Y(1, m);
	//ifstream in("XData.txt");
	//in >> X;
	//in.close();
	//in.open("YData.txt");
	//in >> Y;
	//in.close();

	//double P = 0.0;

	//for (int i = 0; i < 100; i++) {
	//	double h = 1.0 / (1 + exp(-(trans(Real_sol.x) * X[i])()));
	//	if (lroundf(h) == Y(0, i))
	//		h = 1;
	//	else
	//		h = 0;

	//	P += h;
	//}
	//P /= m;

	//cout << Real_sol.x(0, 0) << " " << Real_sol.x(1, 0) << " " << Real_sol.x(2, 0) << " " << Real_sol.y(0, 0) << " " << P << " " << Real_sol.g_calls << endl;
}

void lab5()
{
	////TESTOWANIE FUNKCJI
	//matrix x0 = matrix(2, 1, 5.0);
	//matrix a = matrix(1, 1, 10.0);
	//solution X0;
	//X0.x = x0;
	//X0.fit_fun(ff5T, a, 0.7);

	//cout << X0.x(0) << "--------" << X0.x(1) << "----------" << X0.y(0) << endl;

	//SEKCJA TOSTOWA
	ofstream file("lab_05_tostowe.csv");
	if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
	file << "x0_1, x0_2, x1_a1, x2_a1, f1_a1, f2_a1, f_calls_a1, "
		<< "x1_a2, x2_a2, f1_a2, f2_a2, f_calls_a2, "
		<< "x1_a3, x2_a3, f1_a3, f2_a3, f_calls_a3, \n";

	double a1 = 1.0;
	double a2 = 10.0;
	double a3 = 100.0;

	double epsilon = 1e-3;
	int Nmax = 10000;

	matrix aw = matrix(2, 1);
	matrix x0(2, 1);

	solution powell;
	matrix f1, f2;

	for (double w = 0.00; w <= 1.01; w += 0.01) {
		aw(1) = w;

		x0 = 20 * rand_mat(2, 1) - 10;
		file
			<< x0(0) << ","
			<< x0(1) << ",";

		aw(0) = a1;
		//cout << x0 << endl << endl;
		powell = Powell(ff5T, x0, epsilon, Nmax, aw);
		//cout << powell << endl << endl;
		f1 = ff5T_f1(powell.x, aw);
		f2 = ff5T_f2(powell.x, aw);
		//cout << f1 << endl;
		//cout << f2 << endl;

		file
			<< powell.x(0) << ","
			<< powell.x(1) << ","
			<< f1(0) << ","
			<< f2(0) << ","
			<< solution::f_calls << ",";

		aw(0) = a2;
		//cout << x0 << endl << endl;
		powell = Powell(ff5T, x0, epsilon, Nmax, aw);
		//cout << powell << endl << endl;
		f1 = ff5T_f1(powell.x, aw);
		f2 = ff5T_f2(powell.x, aw);
		//cout << f1 << endl;
		//cout << f2 << endl;

		file
			<< powell.x(0) << ","
			<< powell.x(1) << ","
			<< f1(0) << ","
			<< f2(0) << ","
			<< solution::f_calls << ",";

		aw(0) = a3;
		//cout << x0 << endl << endl;
		powell = Powell(ff5T, x0, epsilon, Nmax, aw);
		//cout << powell << endl << endl;
		f1 = ff5T_f1(powell.x, aw);
		f2 = ff5T_f2(powell.x, aw);
		//cout << f1 << endl;
		//cout << f2 << endl;

		file
			<< powell.x(0) << ","
			<< powell.x(1) << ","
			<< f1(0) << ","
			<< f2(0) << ","
			<< solution::f_calls << ","
			<< "\n";

	}

	file.close();
}

void lab6()
{

}