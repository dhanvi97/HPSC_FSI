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





__global__ void interaction(float *Mi, float *Mie, float *Mei, float *Me, float *Ki, float *Kie, float *Kei, float *Ke, float *u_old_i, float *u_old_e, float *v_old_i, float *v_old_e,float *a_old_i,float *a_old_e, float *Fe, float *Fi, float * A){


int tid = threadIdx.x + blockIdx.x * blockDim.x;

float delta_t = 0.05;
float dot = 0;

float a_new_e[3];
float u_new_e[3];
float v_new_e[3];

float a_new_i[3];
float u_new_i[3];
float v_new_i[3];


for (int j = 0;j <3;j++){
a_new_e[j] = 0;
a_new_i[j] = 0;

u_new_e[j] = 0;
u_new_e[j] = 0;

v_new_e[j] = 0;
v_new_e[j] = 0;
}


for (int iter = 0; iter <2; iter ++){

for (int j = 0 ; j <3 ; j ++){
dot += Ke[tid*3 + j]*u_old_e[j] -Kei[tid*3 +j]*u_old_i[j];
}


a_new_e[tid] = Fe[tid] - dot;
v_new_e[tid] = v_old_e[tid] + (delta_t/2 * (a_old_e[tid] + a_new_e[tid]));
u_new_e[tid] = u_old_e[tid] + delta_t*v_new_e[tid] + (delta_t*delta_t/2)*a_new_e[tid];


float mult = 0;

__shared__ float temp[3];

for (int j = 0; j<3 ;j++){
mult +=  Kie[tid*3+j]*u_new_e[j] +Ki[tid*3+j]*(u_old_i[j] + delta_t*v_old_i[j] +
(delta_t*delta_t/4)*a_old_i[j]) ;
}

temp[tid] = Fi[tid] - mult;

__syncthreads();

for (int j =0; j <3; j++){
a_new_i[tid] = A[tid*3+j]*temp[j];
}

u_new_i[tid] = u_old_i[tid] + delta_t*v_old_i[tid] + (delta_t*delta_t/4)*(a_old_i[tid] + a_new_i[tid]);
v_new_i[tid] = v_old_i[tid] + (delta_t/2)*(a_old_i[tid] + a_new_i[tid]);


for ( int i = 0 ; i <3 ; i ++){
        u_old_e[i] = u_new_e[i];
        v_old_e[i] = v_new_e[i];
		a_old_e[i] = a_new_e[i];

		u_old_i[i] = u_new_i[i];
		v_old_i[i] = v_new_i[i];
		a_old_i[i] = a_new_i[i];
		}
		
__syncthreads();

   }
   
   

}


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

	float u_old_e_host[SIZE_I];
	float u_old_i_host[SIZE_E];
	float v_old_e_host[SIZE_I];
	float v_old_i_host[SIZE_E];
	float a_old_e_host[SIZE_E];
	float a_old_i_host[SIZE_I];
  



	float u_new_e_host[SIZE_I];
	float u_new_i_host[SIZE_E];
	float v_new_e_host[SIZE_I];
	float v_new_i_host[SIZE_E];
	float a_new_e_host[SIZE_E];
	float a_new_i_host[SIZE_I];





	u_old_e << 0, 0, 0;
	u_old_i << 0, 0, 0;
	v_old_i << 0, 0, 0;
	v_old_e << 0, 0, 0;
	a_old_e << 0, 0, 0;
	a_old_i << 0, 0, 0;


	Map<RowMatrixXi>(&u_old_e_host[0], 1, SIZE_E) = u_old_e;
	Map<RowMatrixXi>(&u_old_i_host[0], 1, SIZE_I) = u_old_i;
	Map<RowMatrixXi>(&v_old_e_host[0], 1, SIZE_E) = v_old_e;
	Map<RowMatrixXi>(&v_old_i_host[0], 1, SIZE_I) = v_old_i;
	Map<RowMatrixXi>(&a_old_e_host[0], 1, SIZE_E) = a_old_e;
	Map<RowMatrixXi>(&a_old_i_host[0], 1, SIZE_I) = a_old_i;

	// Map<RowMatrixXi>(&u_new_e_host[0], 1, SIZE_E) = u_new_e;
	// Map<RowMatrixXi>(&u_new_i_host[0], 1, SIZE_I) = u_new_i;
	// Map<RowMatrixXi>(&v_new_e_host[0], 1, SIZE_E) = v_new_e;
	// Map<RowMatrixXi>(&v_new_i_host[0], 1, SIZE_I) = v_new_i;
	// Map<RowMatrixXi>(&a_new_e_host[0], 1, SIZE_E) = a_new_e;
	// Map<RowMatrixXi>(&a_new_i_host[0], 1, SIZE_I) = a_new_i;





	float Fe_host[SIZE_E];
	float Fi_host[SIZE_I];



	Fe << 0, 0.1, 0;
	Fi << 0, 0.2, 0;

	Map<RowMatrixXi>(&Fe_host[0], SIZE_E ,1) = Fe;
	Map<RowMatrixXi>(&Fi_host[0], SIZE_I,1) = Fi;

	//for (int i = 0; i<3 ;i++){
		//cout<<Fe_host[i]<<endl;
	//}



		
		
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

	//float M_host[SIZE][SIZE];
    //float K_host[SIZE][SIZE];

    float Mi_host[9];
    float Mie_host[9];
    float Mei_host[9];
    float Me_host[9];

    float Ki_host[9];
    float Kie_host[9];
    float Kei_host[9];
    float Ke_host[9];
    
    float A_host[9];
	
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
	
	Mi.resize(1,9);
	Mie.resize(1,9);
	Mei.resize(1,9);
	Me.resize(1,9);
	
	Ki.resize(1,9);
	Kie.resize(1,9);
	Kei.resize(1,9);
	Ke.resize(1,9);



	 
	//cout << A<<endl;
	A.resize(1,9);

	//Map<RowMatrixXi>(&M_host[0], SIZE, SIZE) = M;
	//Map<RowMatrixXi>(&K_host[0], SIZE, SIZE) = K;

    Map<RowMatrixXi>(&Mi_host[0], 9, 1) = Mi;
    Map<RowMatrixXi>(&Mei_host[0], 9, 1) = Mei;
    Map<RowMatrixXi>(&Mie_host[0], 9, 1) = Mie;
    Map<RowMatrixXi>(&Me_host[0], 9, 1) = Me;

    Map<RowMatrixXi>(&Ki_host[0], 9, 1) = Ki;
    Map<RowMatrixXi>(&Kei_host[0], 9, 1) = Kei;
    Map<RowMatrixXi>(&Kie_host[0], 9, 1) = Kie;
    Map<RowMatrixXi>(&Ke_host[0], 9,1) = Ke;
    
    Map<RowMatrixXi>(&A_host[0], 9, 1) = A;
    
//    for (int i = 0; i <9; i++){
//  cout<<A_host[i]<<endl;}



    float *u_old_i_device;
    float *u_old_e_device;
    float *v_old_i_device;
    float *v_old_e_device;
    float *a_old_i_device;
    float *a_old_e_device;
    
    float *A_device;
    
    float *Mi_device;
    float *Mie_device;
    float *Mei_device;
    float *Me_device;
    
    float *Ki_device;
    float *Kie_device;
    float *Kei_device;
    float *Ke_device;
    
    float *Fe_device;
    float *Fi_device; 
    
    //cudaMalloc((void**)&Fe_device, 3*sizeof(float));
    //cudaMalloc((void**)&Fe_trial, 3*sizeof(float));
        
    cudaMalloc((void**)&u_old_i_device, 3*sizeof(float));
    cudaMalloc((void**)&u_old_e_device, 3*sizeof(float));
    cudaMalloc((void**)&v_old_i_device, 3*sizeof(float));
    cudaMalloc((void**)&v_old_e_device, 3*sizeof(float));
    cudaMalloc((void**)&a_old_i_device, 3*sizeof(float));
    cudaMalloc((void**)&a_old_e_device, 3*sizeof(float));
    
    
    cudaMalloc((void**)&Mi_device, 9*sizeof(float));
    cudaMalloc((void**)&Mie_device, 9*sizeof(float));
    cudaMalloc((void**)&Mei_device, 9*sizeof(float));
    cudaMalloc((void**)&Me_device, 9*sizeof(float));
    
    
    cudaMalloc((void**)&Ki_device, 9*sizeof(float));
    cudaMalloc((void**)&Kie_device, 9*sizeof(float));
    cudaMalloc((void**)&Kei_device, 9*sizeof(float));
    cudaMalloc((void**)&Ke_device, 9*sizeof(float));
    
    cudaMalloc((void**)&A_device, 9*sizeof(float));
    
    cudaMalloc((void**)&Fe_device, 3*sizeof(float));
    cudaMalloc((void**)&Fi_device, 3*sizeof(float));
          
               
   
    cudaMemcpy(u_old_i_device, &u_old_i_host[0],3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(u_old_e_device, &u_old_e_host[0],3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(v_old_i_device, &v_old_i_host[0],3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(v_old_e_device, &v_old_e_host[0],3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(a_old_i_device, &a_old_i_host[0],3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(a_old_e_device, &a_old_e_host[0],3*sizeof(float),cudaMemcpyHostToDevice);
    
    
    
    cudaMemcpy(Ki_device, &Ki_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Kie_device, &Kie_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Kei_device, &Kei_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Ke_device, &Ke_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
     
        
    cudaMemcpy(Mi_device, &Mi_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Mie_device, &Mie_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Mei_device, &Mei_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Me_device, &Me_host[0],9*sizeof(float),cudaMemcpyHostToDevice);
    
    cudaMemcpy(A_device,&A_host[0], 9*sizeof(float),cudaMemcpyHostToDevice);
    

    cudaMemcpy(Fe_device,&Fe_host[0], 3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(Fi_device,&Fi_host[0], 3*sizeof(float),cudaMemcpyHostToDevice);
    
    
    //interaction<<<1,1>>>(Fe_device,Fe_trial);
    //float trial[3];
    
    
    
    //cudaMemcpy(trial,Fe_trial, 3*sizeof(float),cudaMemcpyDeviceToHost);
    
    //cout<<trial[0]<<endl;
    
    //cout<<trial[1]<<endl;
    //cout<<trial[2]<<endl;
    
    interaction <<<1,3>>>(Mi_device, Mie_device, Mei_device, Me_device, Ki_device, Kie_device, Kei_device, Ke_device, u_old_i_device, u_old_e_device, v_old_i_device, v_old_e_device,
    a_old_i_device,a_old_e_device,Fe_device, Fi_device, A_device);
    
    cudaMemcpy(&a_new_e_host[0], a_old_e_device, 3*sizeof(float), cudaMemcpyDeviceToHost);
    
    for ( int k = 0 ; k < 3; k++)
    {cout <<a_new_e_host[k]<<endl;}



	

/*	for (int i = 0; i < 2; i++) {
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
	*/
	//cout << a_new_e << endl;

	


}
