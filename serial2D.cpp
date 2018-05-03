#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <cmath>
#define N 8
#define K1 1000
#define K2 1

using namespace Eigen;
using namespace std;

int main(){
	int size = N*N;
	MatrixXf M_x(size,size), K_x(size,size), K_stiff(N,N), K_fluid(N,N);
	K_x.setZero(); M_x.setZero(); K_stiff.setZero(); K_fluid.setZero();

	VectorXf K_upper_diag_stiff(N-1), K_diag_stiff(N);
	K_upper_diag_stiff.setConstant(-1*K1);
	K_diag_stiff.setConstant(2*K1);
	K_stiff.diagonal() = K_diag_stiff;
	K_stiff.diagonal(1) = K_upper_diag_stiff;
	K_stiff.diagonal(-1) = K_upper_diag_stiff;

	VectorXf K_upper_diag_fluid(N-1), K_diag_fluid(N);
	K_upper_diag_fluid.setConstant(-1*K1);
	K_diag_fluid.setConstant(2*K1);
	K_diag_fluid.segment(3,N-6).setConstant(2*K2);
	K_upper_diag_fluid.segment(2,N-5).setConstant(-1*K2);
	K_diag_fluid(2) = K1 + K2;
	K_diag_fluid(N-3) = K1 + K2;
	K_fluid.diagonal() = K_diag_fluid;
	K_fluid.diagonal(-1) = K_upper_diag_fluid;
	K_fluid.diagonal(1) = K_upper_diag_fluid;

	for(int i = 0; i < N-1; i++){
		if(i < 3 || i > N-4)
			K_x.block(i*N,i*N,N,N) =  K_stiff;
		else
			K_x.block(i*N,i*N,N,N) = K_fluid;
	}

/*

	float delta_t = 0.05;

	Vector3f u_old_e, v_old_e, a_old_e, 
	u_new_e, v_new_e, a_new_e,
	u_old_i, v_old_i, a_old_i, 
	u_new_i, v_new_i, a_new_i, Fe, Fi ;

	u_old_e << 0, 0, 0;
	u_old_i << 0, 0, 0;
	v_old_i << 0, 0, 0;
	v_old_e << 0, 0, 0;
	a_old_e << 0, 0, 0;
	a_old_i << 0, 0, 0;

	Fe << 0, 0.1, 0;
	Fi << 0, 0.2, 0;
		
		
	M << 1, 0 , 0 , 0 , 0 , 0,
		 0, 1 , 0 , 0 , 0 , 0,
		 0, 0 , 1 , 0 , 0 , 0,
		 0, 0 , 0 , 1 , 0 , 0,
		 0, 0 , 0 , 0 , 1 , 0,
		 0, 0 , 0 , 0 , 0 , 1;

	K << 200, -100, 0, 0, 0, 0,
		 -100, 200, -100, 0, 0, 0,
		 0, -100, 101, -1, 0, 0,
		 0, 0, -1, 2, -1, 0,
		 0, 0, 0, -1, 2, -1,
		 0, 0, 0, 0, -1, 2;
	
	Mi = M.block(0,0,3,3);
	Mie = M.block(3,0,3,3);
	Mei = M.block(0,3,3,3);
	Me = M.block(3,3,3,3);
	Ki = K.block(0,0,3,3);
	Kie = K.block(3,0,3,3);
	Kei = K.block(0,3,3,3);
	Ke = K.block(3,3,3,3);

	A = Mi + Ki*(delta_t*delta_t/4);
	A = A.inverse();

	

	for (int i = 0; i < 10000; i++) {
		a_new_e = (Fe - Ke*u_old_e - Kei*u_old_i);
		v_new_e = v_old_e + (delta_t/2) * ( a_old_e + a_new_e);
		u_new_e = u_old_e + delta_t*v_new_e + (delta_t*delta_t/2)*a_new_e;

		a_new_i = A*(Fi - Kie*u_new_e - Ki*(u_old_i + delta_t*v_old_i + (delta_t*delta_t/4)*a_old_i));
		u_new_i = u_old_i + delta_t*v_old_i + (delta_t*delta_t/4)*(a_old_i + a_new_i);
		v_new_i = v_old_i + (delta_t/2)*(a_old_i+a_new_i);

		u_old_e = u_new_e;
		v_old_e = v_new_e;
		a_old_e = a_new_e;

		u_old_i = u_new_i;
		v_old_i = v_new_i;
		a_old_i = a_new_i;

	}
	cout << a_new_e << endl;
*/
}