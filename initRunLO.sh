#initRun.sh
HomeFolder=/afs/ipp-garching.mpg.de/home/l/lscyboz
GenericFolder=/afs/ipp-garching.mpg.de/home/l/lscyboz/GenericLO
SettingsFolder=/afs/ipp-garching.mpg.de/home/l/lscyboz/Settings
HerwigInput=tT_matchbox_LO.in
GoSamInput=gosamtT.rc

source $HomeFolder/Herwig-7.0.3/bin/activate
export PATH=$HomeFolder/Herwig-7.0.3/bin/:$PATH

mkdir -p $GenericFolder
cd $GenericFolder
cp $SettingsFolder/$HerwigInput $SettingsFolder/$GoSamInput .
Herwig read $HerwigInput

