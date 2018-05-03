#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <cmath>
#include <fstream>

#define N 20
#define K1 200
#define K2 1
#define OMEGA 10
#define NUM_OF_ITER 50

using namespace Eigen;
using namespace std;

MatrixXf Create_Stiffness_Matrix(int n, int size, float k1, float k2) {

	MatrixXf K(size,size), K_fluid(n,n), K_stiff(n,n);
	K_fluid.setZero(); K_stiff.setZero();

	VectorXf K_upper_diag_stiff(n-1), K_diag_stiff(n);
	K_upper_diag_stiff.setConstant(-1*k1);
	K_diag_stiff.setConstant(2*k1);
	K_stiff.diagonal() = K_diag_stiff;
	K_stiff.diagonal(1) = K_upper_diag_stiff;
	K_stiff.diagonal(-1) = K_upper_diag_stiff;

	VectorXf K_upper_diag_fluid(n-1), K_diag_fluid(n);
	K_upper_diag_fluid.setConstant(-1*k1);
	K_diag_fluid.setConstant(2*k1);
	K_diag_fluid.segment(3,n-6).setConstant(2*k2);
	K_upper_diag_fluid.segment(2,n-5).setConstant(-1*k2);
	K_diag_fluid(2) = k1 + k2;
	K_diag_fluid(N-3) = k1 + k2;
	K_fluid.diagonal() = K_diag_fluid;
	K_fluid.diagonal(-1) = K_upper_diag_fluid;
	K_fluid.diagonal(1) = K_upper_diag_fluid;

	for(int i = 0; i < n; i++){
		if(i < 3 || i > N-4)
			K.block(i*n,i*n,n,n) = K_stiff;
		else
			K.block(i*n,i*n,n,n) = K_fluid;
	}

	return K;	
}

MatrixXf Create_Mass_Matrix(int n, int size) {

	MatrixXf M(size,size), M_fluid(n,n), M_stiff(n,n);
	M_fluid.setZero(); M_stiff.setZero();	

	VectorXf M_diag_stiff(n);
	M_diag_stiff.setConstant(1);
	M_stiff.diagonal() = M_diag_stiff;

	VectorXf M_diag_fluid(n);
	M_diag_fluid.setConstant(1);
	M_fluid.diagonal() = M_diag_fluid;

	for(int i = 0; i < n; i++){
		if(i < 3 || i > n-4)
			M.block(i*n,i*n,n,n) = M_stiff;
		else
			M.block(i*n,i*n,n,n) = M_fluid;
	}

	return M;	
}

int main(){
	ofstream output_x, output_y;
	output_x.open("u_x.txt");
	output_y.open("u_y.txt");
	int size = N*N, count_expl = 0, count_impl = 0;
	float delta_t = 1;

	MatrixXf M(size,size), K(size,size);

	K.setZero(); M.setZero();

	K = Create_Stiffness_Matrix(N, size, K1, K2);
	M = Create_Mass_Matrix(N, size);


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

	u_old_e_x.setZero();
	v_old_e_x.setZero();
	a_old_e_x.setZero();
	u_old_i_x.setZero();
	v_old_i_x.setZero();
	a_old_i_x.setZero();

	u_old_e_y.setZero();
	v_old_e_y.setZero();
	a_old_e_y.setZero();
	u_old_i_y.setZero();
	v_old_i_y.setZero();
	a_old_i_y.setZero();

	Fi_x.setZero();
	Fe_x.setZero();
		
	A = M_i + K_i*(delta_t*delta_t/4);
	B = A.inverse();

	MatrixXf u_x(N,N), u_y(N,N);
	MatrixXf buffer = u_new_e_x;
	u_x.setZero(); u_y.setZero();

	for (int i = 0; i < NUM_OF_ITER; i++) {
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
		Fe_x = 0.5*Fe_x;

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
		Fe_y = 0.5*Fe_y;

		for(int i=0; i<N; i++) {
			if(i<3){
				u_x.block(i,0,1,N) = u_new_i_x.segment(i*N,N).transpose();
				u_y.block(0,i,N,1) = u_new_i_y.segment(i*N,N);
			}else if(i < N-3) {
				u_x.block(i,0,1,3) = u_new_i_x.segment(3*N + (i-3)*6,3).transpose();
				u_x.block(i,N-3,1,3) = u_new_i_x.segment(3*N +(i-3)*6+3,3).transpose();
				u_y.block(0,i,3,1) = u_new_i_y.segment(3*N + (i-3)*6, 3);
				u_y.block(N-3,i,3,1) = u_new_i_y.segment(3*N +(i-3)*6+3,3);
			}else{
				u_x.block(i,0,1,N) = u_new_i_x.segment(3*N + (N-6)*6 + (i-N+3)*N, N).transpose();
				u_y.block(0,i,N,1) = u_new_i_y.segment(3*N + (N-6)*6 + (i-N+3)*N, N);
			}
		}
		buffer.resize(N-6,N-6);
		u_x.block(3,3,N-6,N-6) = buffer;
		buffer = u_new_e_y;
		buffer.resize(N-6,N-6);
		u_y.block(3,3,N-6,N-6) = buffer.transpose();

		output_x << u_x <<'\n';
		output_y << u_y <<'\n';
	}
}