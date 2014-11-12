{
  gSystem->Load("libGui.so");

  if ( !gSystem->Load("../build/lib/libCaLib.so") )
    cout << "libCaLib.so is loaded" << endl;

}
