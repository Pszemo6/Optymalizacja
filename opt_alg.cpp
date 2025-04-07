#include"opt_alg.h"
#include <cmath>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
		//Tu wpisz kod funkcji
        solution X0(x0), X1(x0 + d);
        X0.fit_fun(ff, ud1, ud2);
        X1.fit_fun(ff, ud1, ud2);
        if(X0.y == X1.y){
            p[0] = m2d(X0.x);
            p[1] = m2d(X1.x);
            return p;
        }
        
        if(X1.y > X0.y){
            d = -d;
            X1.x = X0.x + d;
            X1.fit_fun(ff, ud1, ud2);
            if (X1.y>=X0.y){
                p[0]=m2d(X1.x);
                p[1]=m2d(X0.x-d);
                return p;
            }
        }
        int i = 0;
        solution Xi = X0;
        while(true){
            if(X0.f_calls > Nmax){
                return NULL;
            }
            i = i+1;
            if(i>1){
                Xi=X1;
            }
            X1.x = X0.x + pow(alpha, i)*d;
            X1.fit_fun(ff, ud1, ud2);
            Xi.fit_fun(ff,ud1,ud2);
            if (Xi.y<=X1.y){
                break;
            }
        }
        solution Xim;
        Xim.x = X0.x + pow(alpha, i-2)*d;
        if (d>0){
            p[0] = m2d(Xim.x);
            p[1] = m2d(X1.x);
            return p;
        }
        p[0]= m2d(X1.x);
        p[1] = m2d(Xim.x);
        return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

double oblicz_fib(int k){
    double x = 0;
    double y = 1;
    
    for(int i = 0; i < k; i++){
        y = y + x;
        x = y - x;
    }
    return y;
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
        solution a0(a), b0(b), c0, d0;
        int k = 0;
        while(oblicz_fib(k) <= (b-a)/epsilon){
            k++;
        }
        c0.x = b0.x - ((oblicz_fib(k-1)/(oblicz_fib(k)))*(b0.x-a0.x));
        d0.x = a0.x + b0.x - c0.x;
        for(int i = 0; i<k-3;i++){
            c0.fit_fun(ff,ud1,ud2);
            d0.fit_fun(ff,ud1,ud2);
            if(c0.y<d0.y){
                a0.x = a0.x;
                b0.x = d0.x;
            }
            else{
                b0.x = b0.x;
                a0.x = c0.x;
            }
            c0.x = b0.x - ((oblicz_fib(k-i-2)/(oblicz_fib(k-i-1)))*(b0.x-a0.x));
            d0.x = a0.x + b0.x - c0.x;
            
        }
        Xopt = c0;
        Xopt.fit_fun(ff, ud1, ud2);
        return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2) {
	try {
		solution::clear_calls();
		solution Ai = a,
			Bi = b,
			Ci = (a + b) / 2,
			Xopt = Ci,
			Di = 0;
		int i = 0;
		double l, m;
		do {
			Ai.fit_fun(ff, ud1, ud2); Bi.fit_fun(ff, ud1, ud2); Ci.fit_fun(ff, ud1, ud2);
			l = Ai.y(0) * (pow(Bi.x(0), 2) - pow(Ci.x(0), 2))
				+ Bi.y(0) * (pow(Ci.x(0), 2) - pow(Ai.x(0), 2))
				+ Ci.y(0) * (pow(Ai.x(0), 2) - pow(Bi.x(0), 2));
			m = Ai.y(0) * (Bi.x(0) - Ci.x(0))
				+ Bi.y(0) * (Ci.x(0) - Ai.x(0))
				+ Ci.y(0) * (Ai.x(0) - Bi.x(0));
			if (m <= 0) {
				Xopt.flag = 0;
				break;
			}
			Di = l / 2 / m;
			Di.fit_fun(ff, ud1, ud2);
			if (Ai.x(0) < Di.x(0) && Di.x(0) < Ci.x(0)) {
				if (Di.y(0) < Ci.y(0)) {
					Bi = Ci;
					Ci = Di;
				}
				else Ai = Di;
			}
			else {
				if (Ci.x(0) < Di.x(0) && Di.x(0) < Bi.x(0)) {
					if (Di.y(0) < Ci.y(0)) {
						Ai = Ci;
						Ci = Di;
					}
					else Bi = Di;
				}
				else {
					Xopt.flag = 0;
					break;
				}
			}
			if (Xopt.f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}
			Xopt = Di;
		} while (Bi.x(0) - Ai.x(0) >= epsilon || abs(Di.x(0) - Xopt.x(0)) >= gamma);
		return Xopt;
	}
	catch (string ex_info) {
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double a, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution::clear_calls();
		double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
		solution X(x0), X_i;
		matrix c0(2, new double[2] { c, a });
		while (true) {
			X_i = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c0);
			if (norm(X_i.x - X.x) < epsilon) {
				return X_i;
			}
			c0(0) = c0(0) * a;
			if (solution::f_calls > Nmax) {
				X_i.flag = -1;
				cout << "Error..." << endl;
				break;
			}
			X = X_i;

		}
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(x0);
		matrix e = ident_mat(n);
		int N = n + 1;
		solution* p = new solution[N];
		solution p_odb, p_z, p_e;
		matrix p_podloga(n, 1);
		int max;
		int min;

		p[0].x = x0;
		p[0].fit_fun(ff, ud1, ud2);

		for (int i = 1; i < N; ++i) {
			p[i].x = p[0].x + s * e[i - 1];
			p[i].fit_fun(ff, ud1, ud2);
		}

		while (true) {
			max = 0;
			min = 0;
			p_podloga = matrix(n, 1);

			for (int i = 1; i < N; ++i) {
				if (p[i].y < p[min].y)
					min = i;
				if (p[i].y > p[max].y)
					max = i;
			}


			for (int i = 0; i < N; ++i) {
				if (i != max)
					p_podloga = p_podloga + p[i].x;
			}

			p_podloga = p_podloga / n;
			p_odb = p_podloga + alpha * (p_podloga - p[max].x);
			p_odb.fit_fun(ff, ud1, ud2);

			if (p_odb.y < p[min].y) {
				p_e = p_podloga + gamma * (p_odb.x - p_podloga);
				p_e.fit_fun(ff, ud1, ud2);
				if (p_e.y < p_odb.y) {
					p[max] = p_e;
				}
				else {
					p[max] = p_odb;
				}
			}
			else {
				if (p[min].y <= p_odb.y && p_odb.y < p[max].y) {
					p[max] = p_odb;
				}
				else {
					p_z = p_podloga + beta * (p[max].x - p_podloga); // Punkt kontrakcji
					p_z.fit_fun(ff, ud1, ud2);
					if (p_z.y >= p[max].y) {
						for (int i = 0; i < N; ++i) {
							if (i != min) {
								p[i] = p[min].x + delta * (p[i].x - p[min].x);
								p[i].fit_fun(ff, ud1, ud2);
							}
						}
					}
					else {
						p[max] = p_z;
					}
				}
			}

			if (solution::f_calls > Nmax) {
				p[min].flag = -1;
				cout << "Error: f_calls > Nmax" << endl;
				break;

			}
			double max_s = 0.0;
			for (int i = 1; i < N; ++i) {
				if (max_s < norm(p[min].x - p[i].x))
					max_s = norm(p[min].x - p[i].x);
			}
			if (max_s < epsilon) {
				return p[min];
			}
		}
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		ofstream file("lab_04_wykresowe_sd.csv");
		if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
		file << "x0, x1\n";


		solution::clear_calls();
		solution X0, Xi;
		X0.x = x0;
		int n = get_len(x0);
		matrix di(n, 1), P(n, 2);
		solution h;
		double* exp_result;
		while (true) {
			di = -X0.grad(gf, ud1, ud2);

			if (h0 == 2137) {
				P.set_col(X0.x, 0);
				P.set_col(di, 1);
				exp_result = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, exp_result[0], exp_result[1], epsilon, Nmax, ud1, P);
				Xi.x = X0.x + h.x * di;
			}
			else {
				Xi.x = X0.x + h0 * di;
			}

			file << Xi.x(0) << "," << Xi.x(1) << "\n";
			//cout << X1.x(0) << endl;
			//cout << X1.x(1) << endl;
			
			if (solution::g_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: f_calls > Nmax" << endl;
				break;
			}
			if (solution::g_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: g_calls > Nmax" << endl;
				break;
			}

			if (norm(Xi.x - X0.x) < epsilon) {
				Xi.fit_fun(ff, ud1, ud2);
				Xi.flag = 1;
				file.close();
				return Xi;
			}
			X0 = Xi;
		}
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		ofstream file("lab_04_wykresowe_cg.csv");
		if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
		file << "x0, x1\n";


		//cout << x0(0) << endl;
		//cout << x0(1) << endl;
		solution::clear_calls();
		int n = get_len(x0);
		solution X0, Xi;
		X0.x = x0;
		matrix d0(n, 1), di(n, 1), P(n, 2);
		solution h;
		double* exp_result, beta;
		d0 = -X0.grad(gf, ud1, ud2);
		di = d0;
		while (true) {
			if (h0 == 2137) {
				P.set_col(X0.x, 0);
				P.set_col(di, 1);
				exp_result = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, exp_result[0], exp_result[1], epsilon, Nmax, ud1, P);
				Xi.x = X0.x + h.x * di;
			}
			else {
				Xi.x = X0.x + h0 * di;
			}

			file << Xi.x(0) << "," << Xi.x(1) << "\n";
			//cout << X1.x(0) << endl;
			//cout << X1.x(1) << endl;

			if (solution::g_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: f_calls > Nmax" << endl;
				break;
			}
			if (solution::g_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: g_calls > Nmax" << endl;
				break;
			}

			if (norm(Xi.x - X0.x) < epsilon) {
				Xi.fit_fun(ff, ud1);
				Xi.flag = 1;
				file.close();
				return Xi;
			}
			Xi.grad(gf);
			beta = pow(norm(Xi.g), 2) / pow(norm(X0.g), 2);
			di = -Xi.g + beta * d0;
			d0 = di;
			X0 = Xi;
		}
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		ofstream file("lab_04_wykresowe_newton.csv");
		if (!file.is_open())throw string("Nie da sie otworzyc pliku csv");
		file << "x0, x1\n";

		//cout << x0(0) << endl;
		//cout << x0(1) << endl;
		solution::clear_calls();
		int n = get_len(x0);
		solution X0, Xi;
		X0.x = x0;
		matrix d(n, 1), P(n, 2);
		solution h;
		double* exp_result;
		while (true) {
			X0.grad(gf);
			X0.hess(Hf);
			d = -inv(X0.H) * X0.g;
			if (h0 == 2137) {
				P.set_col(X0.x, 0);
				P.set_col(d, 1);
				exp_result = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
				h = golden(ff, exp_result[0], exp_result[1], epsilon, Nmax, ud1, P);
				Xi.x = X0.x + h.x * d;
			}
			else {
				Xi.x = X0.x + h0 * d;
			}

			file << Xi.x(0) << "," << Xi.x(1) << "\n";
			//cout << X1.x(0) << endl;
			//cout << X1.x(1) << endl;

			if (solution::f_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: f_calls > Nmax" << endl;
				break;
			}
			if (solution::g_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: g_calls > Nmax" << endl;
				break;
			}
			if (solution::H_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: H_calls > Nmax" << endl;
				break;
			}

			if (norm(Xi.x - X0.x) < epsilon) {
				Xi.fit_fun(ff, ud1);
				Xi.flag = 1;
				file.close();
				return Xi;
			}
			X0 = Xi;
		}
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution A, B, C, D;

		double alfa = (sqrt(5) - 1) / 2;
		A.x = a;
		B.x = b;
		C.x = B.x - alfa * (B.x - A.x);
		D.x = A.x + alfa * (B.x - A.x);

		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);

		while (true) {
			if (C.y < D.y) {
				B = D;
				D = C;
				C.x = B.x - alfa * (B.x - A.x);
				C.fit_fun(ff, ud1, ud2);
			}
			else {
				A = C;
				C = D;
				D.x = A.x + alfa * (B.x - A.x);
				D.fit_fun(ff, ud1, ud2);
			}

			if (solution::f_calls > Nmax) {
				A.flag = -1;
				cout << "Error: f_calls > Nmax" << endl;
				break;
			}

			if (B.x - A.x < epsilon) {
				A.x = (A.x + B.x) / 2;
				A.fit_fun(ff, ud1, ud2);
				A.flag = 1;
				return A;
			}
		}

	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution::clear_calls();

		int n = get_len(x0);
		matrix d = ident_mat(n);
		solution Xi, h;
		solution p;
		Xi.x = x0;
		matrix aw = ud1;
		matrix con(n, 2);
		double* exp_result;

		while (true) {
			p = Xi;
			for (int j = 0; j < n; j++) {
				con.set_col(p.x, 0);
				con.set_col(d[j], 1);
				exp_result = expansion(ff, 0, 1, 1.2, Nmax, aw, con);
				h = golden(ff, exp_result[0], exp_result[1], epsilon, Nmax, aw, con);
				p.x = p.x + h.x * d[j];
			}
			if (norm(p.x - Xi.x) < epsilon) {
				Xi.fit_fun(ff, aw, ud2);
				Xi.flag = 1;
				return Xi;
			}
			for (int j = 0; j < n - 1; ++j) {
				d.set_col(d[j + 1], j);
			}
			d.set_col(p.x - Xi.x, n - 1);
			con.set_col(p.x, 0);
			con.set_col(d[n-1], 1);
			exp_result = expansion(ff, 0, 1, 1.2, Nmax, aw, con);
			h = golden(ff, exp_result[0], exp_result[1], epsilon, Nmax, aw, con);
			Xi.x = p.x + h.x * d[n - 1];


			if (solution::f_calls > Nmax) {
				Xi.flag = -1;
				cout << "Error: f_calls > Nmax" << endl;
				break;
			}
		}


	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
