int read_data_solid(char const * directory = "", int Q2_flag_bin = 4, int z_flag_bin = 4, int target_flag=3, int particle_flag=1, int flag_t=0){/*{{{*/
	//target_flag = 1->"p", 2->"d2" or 3->"3he"
	//particle_flag = 1->pip,2->pim
	//Q2_flag = 1 ~ 4
	//z_flag  = 1 ~ 4
	//flag_t: //1->Colllins, 2->Sivers,3->Pretzelosity

	TString target = "X";
	if(target_flag==1)
		target ="p";
	else if(target_flag==2)
		target ="d2";
	else if(target_flag==3)
		target ="3he";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString particle = "X";
	if(particle_flag==1)
		particle ="pip";
	else if(particle_flag==2)
		particle ="pim";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString filename;

	int ncount=0;

	for(int Q2_flag = 1; Q2_flag <=Q2_flag_bin; Q2_flag++){
		for(int z_flag = 1; z_flag <=z_flag_bin; z_flag++){

			filename.Form("%s/%s_%s_%d_%d.dat",directory,target.Data(),particle.Data(),z_flag,Q2_flag);
            cout <<"--- Reading data from "<< filename.Data()<<endl;
			ifstream infile(filename);



			gStyle->SetOptStat(0);
			int count1,count2;


			//double temp;
			int pt_flag = -1;

			infile >> count1;
			int pt_temp,x_temp;
			while(infile >> pt_flag >> count2){
				for (int j1=0;j1<count2;j1++){
					infile >> pt_temp >> x_temp >> zval[ncount] >> Q2val[ncount] >> ptval[ncount] >> xval[ncount] >> yval[ncount] >> Astat[ncount] >> N[ncount]>>
						coverage[ncount] >> coef[0][ncount] >> coef[1][ncount] >> coef[2][ncount];



					if(pt_temp==pt_flag && x_temp== j1){

						Astat[ncount] *= coef[2-flag_t][ncount]; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity

						//Only use the data points with err<5%

						//if (Astat[ncount]>0.&&Astat[ncount]<0.50){
							//if(VERBOSE1) cout<<Form("--- bin=%d z=%d, Q=%d pt = %f, x = %f, Astat= %f",ncount,z_flag,Q2_flag, ptval[ncount],xval[ncount], Astat[ncount])<<endl;
							ncount ++;
						//}

					}
				}
			}

			infile.close();
		}
	}

	return ncount; // total number of bins

}
/*}}}*/


