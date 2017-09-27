void makeGeom()
{
    new TGeoManager("ArgonCubeGeom", "ArgonCube geometry");

    TGeoMaterial *vacuumMat = new TGeoMaterial("vacuumMat", 0, 0, 0.);
    TGeoMedium *vacuum = new TGeoMedium("vacuum", 0, vacuumMat);
    TGeoMaterial *liquidArgonMat = new TGeoMaterial("liquidArgonMat", 40, 18, 1.4, 14., 83.7);
    TGeoMedium *liquidArgon = new TGeoMedium("liquidArgon", 1, liquidArgonMat);

    const double tpcRadius = 5.05; //cm
    const double tpcLength = 58.6; //cm
    const double worldDx = 4 * tpcRadius;
    const double worldDy = 4 * tpcRadius;
    const double worldDz = 2 * tpcLength;

    TGeoVolume *top = gGeoManager->MakeBox("world", liquidArgon, worldDx, worldDy, worldDz);
    gGeoManager->SetTopVolume(top); // mandatory !

    TGeoVolume *tpc = gGeoManager->MakeTube("tpc", liquidArgon, 0., tpcRadius, (0.5 * tpcLength));
    top->AddNode(tpc, 1, gGeoIdentity);

    //--- close the geometry
    gGeoManager->CloseGeometry();

    //--- draw the ROOT box
    gGeoManager->SetVisLevel(10);
    top->Draw("ogl");

    TFile *outfile = TFile::Open("ACDemoGeom.root", "RECREATE");
    gGeoManager->Write();
    outfile->Close();
}
