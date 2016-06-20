// Reading CLAS data
// TString directory -> directory path 
// int Q2_flag_bin -> number of Q2 bins
// int z_flag_bin-> number of z bins
// int target_flag -> target
// int particle_flag -> particle 
// int flag_t -> 1->Colllins, 2->Sivers,3->Pretzelosity
int read_data_clas(char const * directory = "", int particle_flag=1){
	//particle_flag = 1->pip,2->pim

	TString particle = "X";
	if(particle_flag==1)
		particle ="pip";
	else if(particle_flag==2)
		particle ="pim";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString filename;
	filename.Form("%s/projections_4D_%s_solid.out",directory,particle.Data());
	ifstream infile(filename);

	int ncount=0;
	int bin = -1;
	TString dum;
	infile >> dum >> dum >> dum >> dum >> dum >> dum >> dum >> dum;
	while(!(infile.eof())){
		infile >> bin >> xval[ncount] >> ptval[ncount] >> zval[ncount] >> Q2val[ncount] >> yval[ncount] >> N[ncount]>> Astat[ncount];
		if(Astat[ncount]>1.e-9){
			if(VERBOSE1) cout<<Form("---#%5d: bin=%d z=%f, Q=%f pt = %f, x = %f, Astat= %f",ncount,bin,zval[ncount],Q2val[ncount], ptval[ncount],xval[ncount], Astat[ncount])<<endl;
			ncount ++;
		}
	}
	infile.close();

	return ncount; // total number of bins
}



