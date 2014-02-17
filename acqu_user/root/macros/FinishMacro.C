//#include <iostream>
//#include <string>
//#include <stdio.h>
#include <time.h>

// Get current date/time
/*const char currentDateTime()
{
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// http://www.cplusplus.com/reference/clibrary/ctime/strftime/
	strftime(buf, sizeof(buf), "%d.%m.%Y", &tstruct);

	return buf;
}*/

void saveAnalysis()
{
	time_t now = time(0);
	struct tm tstruct;
	char file[80];
	tstruct = *localtime(&now);
	// http://www.cplusplus.com/reference/clibrary/ctime/strftime/
	strftime(file, sizeof(file), "Analysis_%d.%m.%Y_autosave.root", &tstruct);

	printf("Executing End-of-Run macro, saving analysis to %s\n", file);
	gUAN->SaveAll(file);
  	printf("Done.\n", file);
	//std::cout << "currentDateTime()=" << currentDateTime() << std::endl;  <-- doesn't work in separate method ?!
	//getchar();  // wait for keyboard input
}

void FinishMacro(Char_t* file = NULL)
{
	printf("\nEnd-of-Run macro executing:\n");

	if (!file)
		file = "ARHistograms.root";
	TFile f1(file, "RECREATE");
	gROOT->GetList()->Write();
	f1.Close();
  	printf("done.\n", file);
  	printf("All histograms saved to %s\n", file);
}
