#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <cmath>
#define N 8
#define K1 1000
#define K2 1
#define OMEGA 10

using namespace Eigen;
using namespace std;

MatrixXf Create_Stiffness_Matrix() {
	MatrixXf A(3,3);
	A << 1, 0, 0, 
		 0, 1, 0, 
		 0, 0, 1;
	return A;
}

int main(){
	int size = N*N, count_expl = 0, count_impl = 0;
	MatrixXf M(size,size), K(size,size), K_stiff(N,N), K_fluid(N,N), M_stiff(N,N), M_fluid(N,N);
	MatrixXf Z = Create_Stiffness_Matrix();

	K.setZero(); M.setZero(); K_stiff.setZero(); K_fluid.setZero();

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

	VectorXf M_diag_stiff(N);
	M_diag_stiff.setConstant(1);
	M_stiff.diagonal() = M_diag_stiff;

	VectorXf M_diag_fluid(N);
	M_diag_fluid.setConstant(1);
	M_fluid.diagonal() = M_diag_fluid;

	for(int i = 0; i < N; i++){
		if(i < 3 || i > N-4){
			K.block(i*N,i*N,N,N) = K_stiff;
			M.block(i*N,i*N,N,N) = M_stiff;
		}
		else{
			K.block(i*N,i*N,N,N) = K_fluid;
			M.block(i*N,i*N,N,N) = M_fluid;
		}
	}

	vector<int> ordering;
	int permutation[size];

	for(int i=0; i<size; i++) {
		if(K(i,i) > 200) {
			ordering.push_back(i);
			count_impl++;
		}		
	}

	for(int i=0; i < size; i++) {
		if(K(i,i)<=200) {
			ordering.push_back(i);
			count_expl++;
		}
	}

	for(int i=0; i<size; i++){
		permutation[i] = ordering[i];
	}


	Map<VectorXi> v1(permutation,size);
	PermutationMatrix<Dynamic,Dynamic> perm(v1);
	K = (perm.transpose())*K*perm;
	M = (perm.transpose())*M*perm;


	MatrixXf K_e(count_expl, count_expl), K_i(count_impl, count_impl),
	K_ei(count_expl, count_impl), K_ie(count_impl, count_expl), 
	M_i(count_impl, count_impl), A(count_impl,count_impl), B(count_impl, count_impl);
	
	K_i = K.block(0,0,count_impl,count_impl);
	K_e = K.block(count_impl,count_impl, count_expl, count_expl);
	K_ie = K.block(0, count_impl, count_impl, count_expl);
	K_ei = K.block(count_impl, 0, count_expl, count_impl);
	M_i = M.block(0,0,count_impl,count_impl);

	VectorXf u_old_e_x(count_expl), v_old_e_x(count_expl), a_old_e_x(count_expl),
			 u_new_e_x(count_expl), v_new_e_x(count_expl), a_new_e_x(count_expl),
			 u_old_i_x(count_impl), v_old_i_x(count_impl), a_old_i_x(count_impl),
			 u_new_i_x(count_impl), v_new_i_x(count_impl), a_new_i_x(count_impl), 
			 Fi_x(count_impl), Fe_x(count_expl);
				 
	VectorXf u_old_e_y(count_expl), v_old_e_y(count_expl), a_old_e_y(count_expl),
			 u_new_e_y(count_expl), v_new_e_y(count_expl), a_new_e_y(count_expl),
			 u_old_i_y(count_impl), v_old_i_y(count_impl), a_old_i_y(count_impl),
			 u_new_i_y(count_impl), v_new_i_y(count_impl), a_new_i_y(count_impl), 
			 Fi_y(count_impl), Fe_y(count_expl);
	
	float delta_t = 1;

	u_old_e_x.setZero();
	v_old_e_x.setZero();
	a_old_e_x.setZero();
	u_old_i_x.setZero();
	v_old_i_x.setZero();
	a_old_i_x.setZero();

	Fi_x.setZero();
	Fe_x.setZero();
		
	A = M_i + K_i*(delta_t*delta_t/4);
	B = A.inverse();



	for (int i = 0; i < 10000; i++) {
		a_new_e_x = (Fe_x - K_e*u_old_e_x - K_ei*u_old_i_x);
		v_new_e_x = v_old_e_x + (delta_t/2) * ( a_old_e_x + a_new_e_x);
		u_new_e_x = u_old_e_x + delta_t*v_new_e_x + (delta_t*delta_t/2)*a_new_e_x;
		
		a_new_i_x = B*(Fi_x - K_ie*u_new_e_x - K_i*(u_old_i_x + delta_t*v_old_i_x + (delta_t*delta_t/4)*a_old_i_x));
		u_new_i_x = u_old_i_x + delta_t*v_old_i_x + (delta_t*delta_t/4)*(a_old_i_x + a_new_i_x);
		v_new_i_x = v_old_i_x + (delta_t/2)*(a_old_i_x+a_new_i_x);
		
		u_old_e_x = u_new_e_x;
		v_old_e_x = v_new_e_x;
		a_old_e_x = a_new_e_x;
		
		u_old_i_x = u_new_i_x;
		v_old_i_x = v_new_i_x;
		a_old_i_x = a_new_i_x;
		
		Fe_x.setConstant(OMEGA*i*delta_t);
		Fe_x = cos(Fe_x.array());
		Fe_x = 10*Fe_x;

		a_new_e_y = (Fe_y - K_e*u_old_e_y - K_ei*u_old_i_y);
		v_new_e_y = v_old_e_y + (delta_t/2) * ( a_old_e_y + a_new_e_y);
		u_new_e_y = u_old_e_y + delta_t*v_new_e_y + (delta_t*delta_t/2)*a_new_e_y;
		
		a_new_i_y = B*(Fi_y - K_ie*u_new_e_y - K_i*(u_old_i_y + delta_t*v_old_i_y + (delta_t*delta_t/4)*a_old_i_y));
		u_new_i_y = u_old_i_y + delta_t*v_old_i_y + (delta_t*delta_t/4)*(a_old_i_y + a_new_i_y);
		v_new_i_y = v_old_i_y + (delta_t/2)*(a_old_i_y + a_new_i_y);
		
		u_old_e_y = u_new_e_y;
		v_old_e_y = v_new_e_y;
		a_old_e_y = a_new_e_y;
		
		u_old_i_y = u_new_i_y;
		v_old_i_y = v_new_i_y;
		a_old_i_y = a_new_i_y;
		
		Fe_y.setConstant(OMEGA*i*delta_t);
		Fe_y = sin(Fe_y.array());
		Fe_y = 10*Fe_y;
		
		
	}
	//cout << u_new_i_x << endl;
	cout << u_new_e_x << endl;
}