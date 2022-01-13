double Jpsi_the = g_xithe[index];
double alpha_Jpsi = g_para.alpha_jpsi;
double delta_phi = g_para.phi_jpsi;
double *thC = new double [16];;

double v  = sqrt(1-alpha_Jpsi*alpha_Jpsi);
double st = sin(Jpsi_the);
double ct = cos(Jpsi_the);
double cp = cos(delta_phi);
double sp = sin(delta_phi);
for(int k=0;k<4;k++){
		for(int l=0;l<4;l++){
				double c=0;
				if(k==0&&l==0)c=(1+alpha_Jpsi*ct*ct);
				if(k==0&&l==2)c=v*st*ct*sp;
				if(k==1&&l==1)c=st*st;
				if(k==1&&l==3)c=v*st*ct*cp;
				if(k==2&&l==0)c=-v*st*ct*sp;
				if(k==2&&l==2)c=alpha_Jpsi*st*st;
				if(k==3&&l==1)c=-v*st*ct*cp;
				if(k==3&&l==3)c=-alpha_Jpsi - ct*ct;
				*(thC + k*4 + l)=c;  // jpsi
		}
}

double *tHa =  new double [16];
double alpha = g_para.alpha_xi;
double xiphi = g_para.phi_xi;
double theta = g_lthe[index];
double phi = g_lphi[index];

double gamma = sqrt(1-alpha*alpha)*cos(xiphi);
double beta = sqrt(1-alpha*alpha)*sin(xiphi);
cp = cos(phi);
st = sin(theta);
sp = sin(phi);
ct = cos(theta);
for(int k=0;k<4;k++){
		for(int l=0;l<4;l++){
				double c=0;
				if(k==0&& l == 0) c=1;
				if(k==0&& l == 3) c=alpha;
				if(k==1&& l == 0) c=alpha*cp*st;
				if(k==1&& l == 1) c=gamma*cp*ct - beta*sp;
				if(k==1&& l == 2) c=-(beta*cp*ct) - gamma*sp;
				if(k==1&& l == 3) c=cp*st;
				if(k==2&& l == 0) c=alpha*sp*st;
				if(k==2&& l == 1) c=beta*cp + gamma*ct*sp;
				if(k==2&& l == 2) c=gamma*cp - beta*ct*sp;
				if(k==2&& l == 3) c=sp*st;
				if(k==3&& l == 0) c=alpha*ct;
				if(k==3&& l == 1) c=-(gamma*st);
				if(k==3&& l == 2) c=beta*st;
				if(k==3&& l == 3) c=ct;
				*(tHa + k*4 + l)=c;   //  lambda
		}
}


double *tHb  =  new double [16];
alpha = g_para.alpha_xibar;
xiphi = g_para.phi_xibar;
theta = g_lbthe[index];
phi = g_lbphi[index];

gamma = sqrt(1-alpha*alpha)*cos(xiphi);
beta = sqrt(1-alpha*alpha)*sin(xiphi);
cp = cos(phi);
st = sin(theta);
sp = sin(phi);
ct = cos(theta);
for(int k=0;k<4;k++){
		for(int l=0;l<4;l++){
				double c=0;
				if(k==0&& l == 0) c=1;
				if(k==0&& l == 3) c=alpha;
				if(k==1&& l == 0) c=alpha*cp*st;
				if(k==1&& l == 1) c=gamma*cp*ct - beta*sp;
				if(k==1&& l == 2) c=-(beta*cp*ct) - gamma*sp;
				if(k==1&& l == 3) c=cp*st;
				if(k==2&& l == 0) c=alpha*sp*st;
				if(k==2&& l == 1) c=beta*cp + gamma*ct*sp;
				if(k==2&& l == 2) c=gamma*cp - beta*ct*sp;
				if(k==2&& l == 3) c=sp*st;
				if(k==3&& l == 0) c=alpha*ct;
				if(k==3&& l == 1) c=-(gamma*st);
				if(k==3&& l == 2) c=beta*st;
				if(k==3&& l == 3) c=ct;
				*(tHb + k*4 + l)=c;  //  lambda bar
		}
}



/////////// proton /////////////////////

double *tHc  =  new double [16];
if(g_flag < 2){
alpha = g_para.alpha1_lambda;
}
else {
alpha = g_para.alpha2_lambda;
}
xiphi = 0;
theta = g_pthe[index];
phi = g_pphi[index];

gamma = sqrt(1-alpha*alpha)*cos(xiphi);
beta = sqrt(1-alpha*alpha)*sin(xiphi);
cp = cos(phi);
st = sin(theta);
sp = sin(phi);
ct = cos(theta);
for(int k=0;k<4;k++){
		for(int l = 0; l < 4; l++){
				double c=0;
				if(k==0&&l==0) c=1;
				if(k==1&&l==0) c=alpha*cp*st;
				if(k==2&&l==0) c=alpha*sp*st;
				if(k==3&&l==0) c=alpha*ct;
				 *(tHc + k*4 + l)=c;   // proton
		}
}



/////////// proton bar /////////////////////

double *tHd= new  double [16];
if(g_flag < 2){
alpha = g_para.alpha1_lambdabar;
}
else {
alpha = g_para.alpha2_lambdabar;
}
xiphi = 0;
theta = g_apthe[index];
phi = g_apphi[index];

gamma = sqrt(1-alpha*alpha)*cos(xiphi);
beta = sqrt(1-alpha*alpha)*sin(xiphi);
cp = cos(phi);
st = sin(theta);
sp = sin(phi);
ct = cos(theta);
for(int k=0;k<4;k++){
		for(int l = 0; l < 4; l++){
				double c=0;
				if(k==0&&l==0) c=1;
				if(k==1&&l==0) c=alpha*cp*st;
				if(k==2&&l==0) c=alpha*sp*st;
				if(k==3&&l==0) c=alpha*ct;
				 *(tHd + k*4 + l)=c;   // proton bar
		}
}

	//	double * tep = new double [16];
		for(int mu=0; mu<4;mu++){// Xi loop
				for(int nu=0;nu<4;nu++){// Xibar loop
						for(int k=0;k<4;k++){
								for(int j=0;j<4;j++){
										g_eval[index]+=thC[mu*4 +nu]*
												tHa[mu*4 + k]*tHb[nu*4+j]*tHc[k*4+0]*tHd[j*4 +0];
								}
						}
				}
		}
//	printf("from gpu: %f\n", g_eval[index]);
		delete [] tHa;
		delete [] tHb;
		delete [] thC;
		delete [] tHc;
		delete [] tHd;
//for(int i =0 ; i < 16; i++){
//		g_eval[index] += *(tep + i)*tHa[i%4]*tHb[i/4];
//}
		
//delete [] tep;


/*

		double *a_c = new double [16];

		for(int j=0; j<4;j++){// Xi loop
				for(int nu=0;nu<4;nu++){// Xibar loop
						for(int mu=0;mu<4;mu++){
								*(a_c+j*4 +nu) += *(tHa + j*4 +mu)* *(thC + mu *4 + nu );
						}
				}
		}

		double *a_c_ab = new double [16];
		for(int j=0; j<4;j++){// Xi loop
				for(int k=0;k<4;k++){// Xibar loop
						for(int nu=0;nu<4;nu++){
								*(a_c_ab + j*4 + k) += *(a_c + j*4 +nu)* (*(tHb + nu *4 + k ));
						}
				}
		}
		delete [] tHa;
		delete [] tHb;
		delete [] thC;


for(int k=0;k<4;k++){
		for(int j=0;j<4;j++){
			tep[k*4 +j] = 	tHc[k*4+0]*a_c_ab[k*4 + j]*tHd[j*4 +0];
		}
}

		delete [] tHc;
		delete [] tHd;
		delete [] a_c;
		delete [] a_c_ab;
for(int i =0 ; i < 16; i++){
		g_eval[index] += *(tep + i);
}
delete [] tep;
*/
