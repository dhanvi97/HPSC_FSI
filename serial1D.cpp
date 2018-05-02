#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <cmath>
#define N 1024

using namespace Eigen;
using namespace std;

int main(){
	MatrixXf M(6,6), K(6,6), 
	Me(3,3), Mei(3,3), Mie(3,3), Mi(3,3), 
	Ke(3,3), Kei(3,3), Kie(3,3), Ki(3,3), A(3,3);

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

	

	for (int i = 0; i < 2; i++) {
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

	


}