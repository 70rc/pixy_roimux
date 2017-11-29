#include <iostream>
#include <string>
#include <cerrno>

#include <TFile.h>
#include <TGeoManager.h>
#include <TTree.h>

#include <ConstField.h>
#include <EventDisplay.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <Track.h>


int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << program_invocation_short_name << " genfitTreeFileName geoFileName" << std::endl;
        return 1;
    }
    const std::string genfitTreeFileName = std::string(argv[1]);
    const std::string geoFileName = std::string(argv[2]);

    new TGeoManager("Geometry", "GenFit Geometry");
    TGeoManager::Import(geoFileName.c_str());
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., 0.));
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface);

    genfit::EventDisplay *display = genfit::EventDisplay::getInstance();

    TFile genfitTreeFile(genfitTreeFileName.c_str(), "READ");
    if (!genfitTreeFile.IsOpen()) {
        std::cerr << "ERROR: Failed to open genfitTreeFile " << genfitTreeFileName << '!' << std::endl;
        return 1;
    }
    TTree *genfitTree = nullptr;
    genfitTreeFile.GetObject("genfitTree", genfitTree);
    if (!genfitTree) {
        std::cerr << "ERROR: Failed to get genfitTree form file " << genfitTreeFileName << '!' << std::endl;
        return 1;
    }
    unsigned eventId;
    genfit::Track *track = nullptr;
    genfitTree->SetBranchAddress("eventId", &eventId);
    genfitTree->SetBranchAddress("Track", &track);

    for (unsigned eventIdx = 0; eventIdx < genfitTree->GetEntries(); ++eventIdx) {
        genfitTree->GetEntry(eventIdx);
        std::cout << "Adding event " << eventId << "...\n";
        display->addEvent(track);
    }

    genfitTreeFile.Close();

    display->open();

    return 0;
}