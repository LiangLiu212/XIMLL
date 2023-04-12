//  Constraint alpha_xi

void fcnMLLGISO(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
		double llf = 0.0;

		double pp[10]={0};
		pp[0] = par[0];
		pp[1] = par[1];

		pp[2] = par[2] * par[3] + par[2];
		pp[3] = par[5] + par[4];

		pp[4] = par[2] * par[3] - par[2];
		pp[5] = par[5] - par[4];

		pp[6] = (par[6] - par[8]) * par[7] + (par[6] - par[8]);
		pp[7] = (par[6] + par[8]) * par[9] - (par[6] + par[8]);
		pp[8] = (par[6] + par[8]) * par[9] + (par[6] + par[8]);
		pp[9] = (par[6] - par[8]) * par[7] - (par[6] - par[8]);

		llf = rf->fcnmll(pp);
		if(llf == 0){ cout << "ERROR" << endl; return;}
		f =  llf;
}
// using minuit2 to minimize the likelihood function
double RosenbrockISO(const double * par){
		double pp[10]={0};
		pp[0] = par[0];
		pp[1] = par[1];

		pp[2] = par[2] * par[3] + par[2];
		pp[3] = par[5] + par[4];

		pp[4] = par[2] * par[3] - par[2];
		pp[5] = par[5] - par[4];

		pp[6] = (par[6] - par[8]) * par[7] + (par[6] - par[8]);
		pp[7] = (par[6] + par[8]) * par[9] - (par[6] + par[8]);
		pp[8] = (par[6] + par[8]) * par[9] + (par[6] + par[8]);
		pp[9] = (par[6] - par[8]) * par[7] - (par[6] - par[8]);

		double llf = 0.0;
		llf = rf->fcnmll(pp);
		if(llf == 0){ cout << "ERROR" << endl; exit(0);}
		return llf;
}

//=====================================================================
// using minuit2 to minimize the likelihood function
void Minuit2XiXiMLLISO(int index, int MM, TString outfile_name = "out.txt"){

		ofstream out;
		outfile_name = "Minuit2" + outfile_name;
		cout << outfile_name << endl;
		out.open(outfile_name, ios::out | ios::app);

		cout << "OK 11111111113" << endl;
		double Jpsi_alpha   = 0.586;  	  // alpha_J/Psi 
		double Jpsi_phi     = 1.121;		  //-TMath::Pi()/4.; // relative phase, Dphi_J/Psi
		double xi_alpha     = -0.3756;  		  // alpha (Sgm->p pi0)
		double xi_Acp       = 0.0;   			  // alpha (Sgm->pbar pi0)
		double xi_phi      = 0.012;   
		double xi_phicp    = 0.0;  
		double L_alpha0     = 0.700;   
		double L_A0cp       = 0.0;   
		double L_alpham     = 0.0;   
		double L_Amcp       = 0.0;

		cout << "OK" << endl;
		// cout << argv[1] << endl;
		// fit nr is used to tell which analysis cuts that are used

		const char * algoName = "migrad";
		const int printlevel = 1;
	//	ROOT::Math::Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer(algoName);
		ROOT::Minuit2::Minuit2Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer(algoName);
		min->SetMaxFunctionCalls(1000000);
		min->SetTolerance(0.001);
		min->SetPrintLevel(printlevel);
		ROOT::Math::Functor f(&RosenbrockISO,10);
		min->SetFunction(f);
		min->SetErrorDef(0.5);

		double step[10] = {0};
		for(int i = 0; i < 10;i++){
				step[i] = 0.001;
		}
		min->SetVariable(0,"alpha_jpsi", Jpsi_alpha, step[0]);
		min->SetVariable(1,"dphi_jpsi",  Jpsi_phi,   step[1]);
		min->SetVariable(2,"xi_alpha",   xi_alpha,   step[2]);
		min->SetVariable(3,"xi_Acp",     xi_Acp,     step[3]);
		min->SetVariable(4,"xi_phi",     xi_phi,     step[4]);
		min->SetVariable(5,"xi_phicp",   xi_phicp,   step[5]);
		min->SetVariable(6,"L_alpha0",   L_alpha0,   step[6]);
		min->SetVariable(7,"L_A0cp",     L_A0cp,     step[7]);
		min->SetVariable(8,"L_alpham",   L_alpham,   step[8]);
		min->SetVariable(9,"L_Amcp",     L_Amcp,     step[9]);
		min->SetVariableLimits(0, -1., 1.);
		min->SetVariableLimits(1, -TMath::Pi(), TMath::Pi());
		min->SetVariableLimits(2, -1., 0.);
		min->SetVariableLimits(3, -0.1, 0.1);
		min->SetVariableLimits(4, -TMath::Pi(), TMath::Pi());
		min->SetVariableLimits(5, -0.1, 0.1);
		min->SetVariableLimits(6, 0., 1.);
		min->SetVariableLimits(7, -0.1, 0.1);
		min->SetVariableLimits(8, -1., 1.);
		min->SetVariableLimits(9, -0.1, 0.1);
		
		if(fixAcpXi){
				min->FixVariable(3);
		}
		if(fixphicpXi){
				min->FixVariable(5);
		}
		if(fixAcp0){
				min->FixVariable(7);
		}
		if(fixAcpm){
				min->FixVariable(9);
		}

		min->FixVariable(8);

		min->Minimize();

		
		cout << "Cov Matrix" << endl;
		for(int i = 0; i < 10; i++){
				for(int j = 0; j< 10; j++){
						cout << Form("%10.3e", min->CovMatrix(i, j)) << "  ";
				}
				cout << endl;
		}
	
		cout << "Corr. Matrix" << endl;
		for(int i = 0; i < 10; i++){
				for(int j = 0; j< 10; j++){
						cout << Form("%10.3f", min->Correlation(i, j)) << "  ";
				}
				cout << endl;
		}



//#define minuit2_test
#ifndef minuit2_test

		cout << "OK 11111111113" << endl;
		double res[10], err_res[10];
		double low[10], high[10];
		for(int i = 0; i < 10; i++){
				low[i] = 0;
				high[i] = 0;
				min->GetMinosError(i, low[i], high[i]);
				res[i] = min->State().Parameter(i).Value();
				err_res[i] = (fabs(low[i]) + fabs(high[i]))/2.0;
		}

		double fmin = min->State().Fval();
		int istat = min->Status();
		out << fmin << "," << istat << ","; 
		out << res[0]<< "," << err_res[0] << "," << res[1]<< "," << err_res[1]<< ","; 
		out << res[2]<< "," << err_res[2] << "," << res[3]<< "," << err_res[3]<< ",";
		out << res[4]<< "," << err_res[4] << "," << res[5]<< "," << err_res[5]<< ",";
		out << res[6]<< "," << err_res[6] << "," << res[7]<< "," << err_res[7]<< ",";
		out << res[8]<< "," << err_res[8] << "," << res[9]<< "," << err_res[9]<< endl;
	
		for(int i = 0; i < 10; i++){
				for(int j = 0; j< 10; j++){
						out << min->CovMatrix(i, j) << ",";
				}
				out << endl;
		}
		out.close();
		//	return 0;
#endif
}

//=====================================================================
// input [1] =  0; [2] =  type; [3] = step; [4] = output file
void XiXiMLLISO(int index, int MM, TString outfile_name = "out.txt"){

		ofstream out;
		cout << outfile_name << endl;
		outfile_name = "Minuit" + outfile_name;
		out.open(outfile_name, ios::out | ios::app);

		cout << "OK 11111111113" << endl;

		double Jpsi_alpha   = 0.586;  	  // alpha_J/Psi 
		double Jpsi_phi     = 1.121;		  //-TMath::Pi()/4.; // relative phase, Dphi_J/Psi
		double xi_alpha     = -0.3756;  		  // alpha (Sgm->p pi0)
		double xi_Acp       = 0.0;   			  // alpha (Sgm->pbar pi0)
		double xi_phi      = 0.012;   
		double xi_phicp    = 0.0;  
		double L_alpha0     = 0.70;   
		double L_A0cp       = 0.0;   
		double L_alpham     = 0.0;   
		double L_Amcp       = 0.0;

		cout << "OK" << endl;
		// cout << argv[1] << endl;
		// fit nr is used to tell which analysis cuts that are used
		myMinuit *minuit=new myMinuit(10);
	//	if(!isIO){
	//			minuit->setRandomSeed(3423);
	//	}
		Int_t ierflag=0; 
		Double_t arglist[100];
		cout << "OK 11111111111" << endl;
		minuit->SetFCN(fcnMLLGISO);
		cout << "OK 11111111111" << endl;
		arglist[0]= 0;
		minuit->mnexcm("SET PRINT",arglist,1,ierflag);
		arglist[0]= 0.5;
		minuit->mnexcm("SET ERR",arglist,1,ierflag);
		minuit->mnparm(0, "alpha_jpsi" ,Jpsi_alpha, 0.001, -1., 1., ierflag);
		minuit->mnparm(1, "dphi_jpsi", Jpsi_phi, 0.001, -TMath::Pi(), TMath::Pi(), ierflag);
		cout << "OK 11111111111" << endl;
		minuit->mnparm(2, "xi_alpha" , xi_alpha, 0.001, -1., 0., ierflag);
		minuit->mnparm(3, "xi_Acp" , xi_Acp, 0.001,  -0.1, 0.1, ierflag);
		minuit->mnparm(4, "xi_phi" , xi_phi, 0.001, -TMath::Pi(), TMath::Pi(), ierflag);
		minuit->mnparm(5, "xi_phicp" , xi_phicp, 0.001,  0.5, 0.5, ierflag);
		minuit->mnparm(6, "L_alpha0" , L_alpha0, 0.001, 0., 1., ierflag);
		minuit->mnparm(7, "L_A0cp" , L_A0cp, 0.001, -0.1, 0.1, ierflag);
		minuit->mnparm(8, "L_alpham" , L_alpham, 0.001, -1., 1., ierflag);
		minuit->mnparm(9, "L_Amcp" , L_Amcp, 0.001, -0.1, 0.1, ierflag);
		cout << "OK 11111111113" << endl;

		if(fixAcpXi){
				minuit->FixParameter(3);
		}
		if(fixphicpXi){
				minuit->FixParameter(5);
		}
		if(fixAcp0){
				minuit->FixParameter(7);
		}
		if(fixAcpm){
				minuit->FixParameter(9);
		}
		minuit->FixParameter(8);

		// 	 		minuit->FixParameter(7);
		//		 minuit->FixParameter(12);
		//		 minuit->FixParameter(13);
		cout << "OK 11111111114" << endl;
		minuit->mnexcm("MINI",arglist,0,ierflag); //minimization using the migrag
		cout << "OK 11111111115" << endl;
		//limits both 0 implies no limit 
		minuit->mnexcm("MINOS",arglist,0,ierflag);
		cout << "OK 11111111115" << endl;
		minuit->mnmatu(1);
		cout << "OK 11111111115" << endl;
		Double_t fmin, fedm, errdef;
		Int_t   npari, nparx, istat; 
		minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
		double res[10], err_res[10];
		for(int p = 0; p < 10; p++)
				minuit->GetParameter(p, res[p], err_res[p]);
		out << fmin << "," << istat << ","; 
		out << res[0]<< "," << err_res[0] << "," << res[1]<< "," << err_res[1]<< ","; 
		out << res[2]<< "," << err_res[2] << "," << res[3]<< "," << err_res[3]<< ",";
		out << res[4]<< "," << err_res[4] << "," << res[5]<< "," << err_res[5]<< ",";
		out << res[6]<< "," << err_res[6] << "," << res[7]<< "," << err_res[7]<< ",";
		out << res[8]<< "," << err_res[8] << "," << res[9]<< "," << err_res[9]<< endl;
		out.close();
		//	return 0;
}


