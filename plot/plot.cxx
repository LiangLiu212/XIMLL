/*
 * =====================================================================================
 *
 *       Filename:  plot.cxx
 *
 *    Description: plot 
 *
 *        Version:  1.0
 *        Created:  01/24/2022 20:48:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liu Liang (LL), liangzy@mail.ustc.edu.cn
 *   Organization:  USTC
 *
 * =====================================================================================
 */

void plot(){
	std::ifstream input( "../run4/out.txt" );
	double mllv;
	int status, nevent1, nevent2;
	double ipara[10];
	double iparaerr[10];
	double para[10][30];
	double paraerr[10][30];
	int j = 0;
	double paramean[10] = {0.586, 1.213, -0.375, 0.02, 0.375, -0.02, 0.692, -0.757, 0.757, -0.692};
	TH1D *hnom[10];
	for(int i = 0; i< 10; i++){
		TString name = Form("hnom%d", i+1);
		if(i == 1){
		hnom[i] = new TH1D(name, name, 9, -5, 5);
		}
		else{
		hnom[i] = new TH1D(name, name, 11, -5, 5);
		}
	}

	TString title[10] = {"#alpha_{J/#psi}", 
	"#Delta#Phi_{J/#psi}",
	"#alpha_{#Xi}",
	"#phi_{#Xi}",
	"#alpha_{#bar{#Xi}}",
	"#phi_{#bar{#Xi}}",
	"#alpha_{0}",
	"#alpha_{+}",
	"#alpha_{-}",
	"#bar{#alpha}_{0}"
	};



	for( std::string line; getline( input, line ); ){
		istringstream iss(line);
	//	cout << line << endl;
		iss >> mllv >> status >> nevent1>> nevent2;
		for(int i = 0; i < 10; i++){
			iss >> ipara[i] >> iparaerr[i];
			para[i][j] = ipara[i];
			paraerr[i][j] = iparaerr[i];
		}
		j++;
	}
	TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
	c1->Divide(4, 3);
	TF1* f1 = new TF1("f1", "gaus",  -5, 5);
	for(int k = 0; k < 10; k++){
		for(int i=0; i < 29; i++){
			hnom[k]->Fill((para[k][i] - paramean[k])/paraerr[k][i]);
		}
		c1->cd(k+1);
		hnom[k]->Draw();
		hnom[k]->GetXaxis()->SetTitle(title[k]);
		hnom[k]->SetLineColor(2);
		hnom[k]->Fit("f1");
		hnom[k]->GetXaxis()->CenterTitle(true);
		hnom[k]->GetXaxis()->SetTitleSize(0.10);
		hnom[k]->GetXaxis()->SetTitleOffset(0.80);
		hnom[k]->GetYaxis()->SetTitleOffset(0.80);
		hnom[k]->GetYaxis()->CenterTitle(true);
		hnom[k]->GetYaxis()->SetTitleSize(0.10);
		hnom[k]->GetYaxis()->SetTitle("Samples");
		f1->SetLineColor(4);
	}

}
