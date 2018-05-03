#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include <cmath>
#define N 1024
#define SIZE 6
#define SIZE_I 3
#define SIZE_E 3

using namespace Eigen;
using namespace std;
typedef Matrix<float, Dynamic, Dynamic, RowMajor> RowMatrixXi;
typedef Eigen::VectorXd RealVector;



// __global__ void interaction(){
//    int ROW = blockIdx.y*blockDim.y+threadIdx.y;
//    int COL = blockIdx.x*blockDim.x+threadIdx.x;
// }



void print_arr(float print_arr[][SIZE], int size_row, int size_column){
	for (int i = 0 ; i <size_row; i++){
		for (int j =0; j<size_column; j++){
			cout<<print_arr[i][j];
		}
			cout<<endl;
		}
	}

int main(){




	MatrixXf M(6,6), K(6,6), 
	Me(3,3), Mei(3,3), Mie(3,3), Mi(3,3), 
	Ke(3,3), Kei(3,3), Kie(3,3), Ki(3,3), A(3,3);

	float delta_t = 0.005;

	Vector3f u_old_e, v_old_e, a_old_e, 
	u_new_e, v_new_e, a_new_e,
	u_old_i, v_old_i, a_old_i, 
	u_new_i, v_new_i, a_new_i, Fe, Fi ;

	float u_old_e_host[1][SIZE_I];
	float u_old_i_host[1][SIZE_E];
	float v_old_e_host[1][SIZE_I];
	float v_old_i_host[1][SIZE_E];
	float a_old_e_host[1][SIZE_E];
	float a_old_i_host[1][SIZE_I];
  



	float u_new_e_host[1][SIZE_I];
	float u_new_i_host[1][SIZE_E];
	float v_new_e_host[1][SIZE_I];
	float v_new_i_host[1][SIZE_E];
	float a_new_e_host[1][SIZE_E];
	float a_new_i_host[1][SIZE_I];





	u_old_e << 0, 0, 0;
	u_old_i << 0, 0, 0;
	v_old_i << 0, 0, 0;
	v_old_e << 0, 0, 0;
	a_old_e << 0, 0, 0;
	a_old_i << 0, 0, 0;


	Map<RowMatrixXi>(&u_old_e_host[0][0], 1, SIZE_E) = u_old_e;
	Map<RowMatrixXi>(&u_old_i_host[0][0], 1, SIZE_I) = u_old_i;
	Map<RowMatrixXi>(&v_old_e_host[0][0], 1, SIZE_E) = v_old_e;
	Map<RowMatrixXi>(&v_old_i_host[0][0], 1, SIZE_I) = v_old_i;
	Map<RowMatrixXi>(&a_old_e_host[0][0], 1, SIZE_E) = a_old_e;
	Map<RowMatrixXi>(&a_old_i_host[0][0], 1, SIZE_I) = a_old_i;

	Map<RowMatrixXi>(&u_new_e_host[0][0], 1, SIZE_E) = u_new_e;
	Map<RowMatrixXi>(&u_new_i_host[0][0], 1, SIZE_I) = u_new_i;
	Map<RowMatrixXi>(&v_new_e_host[0][0], 1, SIZE_E) = v_new_e;
	Map<RowMatrixXi>(&v_new_i_host[0][0], 1, SIZE_I) = v_new_i;
	Map<RowMatrixXi>(&a_new_e_host[0][0], 1, SIZE_E) = a_new_e;
	Map<RowMatrixXi>(&a_new_i_host[0][0], 1, SIZE_I) = a_new_i;





	float Fe_host[SIZE_E];
	float Fi_host[1][SIZE_I];



	Fe << 0, 0.1, 0;
	Fi << 0, 0.2, 0;

	Map<RowMatrixXi>(&Fe_host[0], SIZE_E ,1) = Fe;
	Map<RowMatrixXi>(&Fi_host[0][0], SIZE_I,1) = Fi;

	for (int i = 0; i<3 ;i++){
		cout<<Fe_host[0][i]<<endl;
	}



		
		
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

	float M_host[SIZE][SIZE];
    float K_host[SIZE][SIZE];

    float Mi_host[SIZE_I][SIZE_I];
    float Mie_host[SIZE_I][SIZE - SIZE_I];
    float Mei_host[SIZE_E][SIZE - SIZE_E];
    float Me_host[SIZE_E][SIZE_E];

    float Ki_host[SIZE_I][SIZE_I];
    float Kie_host[SIZE_I][SIZE - SIZE_I];
    float Kei_host[SIZE_E][SIZE - SIZE_E];
    float Ke_host[SIZE_E][SIZE_E];
	
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

	Map<RowMatrixXi>(&M_host[0][0], SIZE, SIZE) = M;
	Map<RowMatrixXi>(&K_host[0][0], SIZE, SIZE) = K;

    Map<RowMatrixXi>(&Mi_host[0][0], SIZE_I, SIZE_I) = Mi;
    Map<RowMatrixXi>(&Mei_host[0][0], SIZE_E, SIZE - SIZE_E) = Mei;
    Map<RowMatrixXi>(&Mie_host[0][0], SIZE_I, SIZE - SIZE_I) = Mie;
    Map<RowMatrixXi>(&Me_host[0][0], SIZE_E, SIZE_E) = Me;

    Map<RowMatrixXi>(&Ki_host[0][0], SIZE_I, SIZE_I) = Ki;
    Map<RowMatrixXi>(&Kei_host[0][0], SIZE_E, SIZE - SIZE_E) = Kei;
    Map<RowMatrixXi>(&Kie_host[0][0], SIZE_I, SIZE - SIZE_I) = Kie;
    Map<RowMatrixXi>(&Ke_host[0][0], SIZE_E, SIZE_E) = Ke;

    size_t bytes = SIZE;
    size_t bytes_i = SIZE_I*sizeof(float);
    size_t bytes_e = SIZE_E*sizeof(float);
  

    int *u_new_i_device;
    cudaMalloc((void**)&u_new_i_device, 3*sizeof(float));

    cudaMemcpy(u_new_i_device,u_new_i_host,cudaMemcpyHostToDevice);



	

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