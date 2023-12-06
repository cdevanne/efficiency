// rootlogon.C
{
    // Chemin absolu du répertoire du script
    TString scriptDir = gSystem->ExpandPathName(gSystem->Getenv("PWD"));

    // Charge automatiquement tous les fichiers .h dans le répertoire include
    TString includeDir = scriptDir + "/../include/";
    gSystem->AddIncludePath("-I" + includeDir);
    
    // Charge automatiquement toutes les bibliothèques partagées
    TString libDir = scriptDir + "/../lib/";
    gSystem->Load(libDir + "libEfficiency.so");
    
    // Autres configurations spécifiques au projet
    R__ADD_INCLUDE_PATH(../include/)
}
