void Ped1885(int cstart, int nx = 96)
{
// Pedestals for LRS 1885F
// Modified JRMA 8/1/12 to create config file for AcquDAQ
// The auto created file needs to be checked
// as this macro has not been tested for AcquDAQ
    char hist_name[16];
    char file_name[32];
    //    char h[8];
    TH1F* hist;
    TAxis* xaxis;
    Double_t mean;
    FILE* pfile;
    //
    sprintf(file_name,"lrs_1885.adc%d",cstart);
    printf("Pedestals saved to file %s\n",file_name);
    if( (pfile = fopen(file_name,"w")) == NULL ){
      printf("File open failed\n");
      return;
    }
    // Print the non-pedestal spec. part of the config file.
    // This will normally have a pretty standard format
    fprintf(pfile,"*     Auto-Generated Pedestal File for LRS 1885 QDC\n");
    fprintf(pfile,"*     ROOT Macro Ped1885, JRM Annand 08/01/12\n");
    fprintf(pfile,"BaseSetup: 0  3 VADC\n");
    fprintf(pfile,"BaseIndex: 0  96  12\n");
    fprintf(pfile,"Hardware-ID:	0x1045\n");
    fprintf(pfile,"Register: 0x0     80000000\n");
    fprintf(pfile,"Register: 0x0     08000000\n");
    fprintf(pfile,"Register: 0x1     5f\n");
    // Now print the pedestal specification lines
    for (Int_t t=cstart,i=0; t<(cstart+nx); t++,i++)
    {
	sprintf(hist_name,"ADC%d",t);
	hist = (TH1F)gROOT->Get(hist_name);
	xaxis = (TAxis)hist->GetXaxis();
//              xaxis->SetRange(0, 4096);
	mean = hist->GetMean();
	fprintf(pfile,"Register: pd%d     %d\n",i,mean);
    }
    // One further line...don't remember if this is necessary
    // needs to be checked
    fprintf(pfile,"RD-Mode:\n");
    fclose(pfile);
    return;
}
//                                           
