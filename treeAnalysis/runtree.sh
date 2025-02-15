input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_pTHard5")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB_pTHard5")'


input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_epMB.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-TREXTOUT_epMB")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-TREXTOUT_pTHard5")'
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_epMB_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","EVTTREE-ALLSILICON-TREXTOUT_epMB_pTHard5")'


if [ $1 = "DRpion" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_CALOSTANDALONE_SimplePhoton.root
    root -x -l -b -q 'treeProcessing.C("'$input'","CALOSTANDALONE_SimplePhoton",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
fi

if [ $1 = "calogamma" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_CALOSTANDALONE_SimplePhoton.root
    root -x -l -b -q 'treeProcessing.C("'$input'","CALOSTANDALONE_SimplePhoton",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
elif [ $1 = "calopion" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_CALOSTANDALONE_SimplePion.root
    root -x -l -b -q 'treeProcessing.C("'$input'","CALOSTANDALONE_SimplePion",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
elif [ $1 = "fhcalpion" ]; then
    input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/tree_FHCALSTANDALONE_SimplePion.root
    root -x -l -b -q 'treeProcessing.C("'$input'","FHCALSTANDALONE_SimplePion",true,true,true,true,true,true,true,true,true,true,true,true,false, true, true, -1, 0, 1, false)'
fi

# MODULAR
if [ $1 = "1" ]; then
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_epMB_pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "2" ]; then
input=/media/nschmidt/external3/EIC_outputs2/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh.root
root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "3" ]; then
input=/media/nschmidt/external2/EICsimulationOutputs/forFredi/EVTTREE-ALLSILICON-TREXTOUT_epMB_pTHard5.root
root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "4" ]; then
input=/media/nschmidt/external3/EIC_outputs2/Modular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT-HC2xEC2x_e10p250pTHard5.root
root -x -l -b -q 'treeProcessing.C("'$input'","MODULAR_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT-HC2xEC2x_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true, true, true, -1, 0, 2)'

# BEAST MODULAR
elif [ $1 = "5" ]; then
input=/media/nschmidt/external3/EIC_outputs2/BeastModular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5.root
#root -x -l -b -q 'treeProcessing.C("'$input'","BEAST_ALLSILICON-FTTLS3LC-ETTL-CTTL-ACLGAD-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "6" ]; then
input=/media/nschmidt/external3/EIC_outputs2/BeastModular/output_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh.root
root -x -l -b -q 'treeProcessing.C("'$input'","BEAST_ALLSILICON-FTTLS3LC-ETTL-CTTL-TREXTOUT_e10p250MBandpTh",true,true,true,true,true,true,true,true,true,true,true,true,true)'
elif [ $1 = "7" ]; then
input=/media/nschmidt/external3/EIC_outputs2/BeastModular/output_ALLSILICON-TREXTOUT_e10p250pTHard5.root
root -x -l -b -q 'treeProcessing.C("'$input'","BEAST_ALLSILICON-TREXTOUT_e10p250pTHard5",true,true,true,true,true,true,true,true,true,true,true,true,true)'
fi

if [ $1 = "newinput_inlay" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-ETTL-CTTL-INNERTRACKING-DRCALO-FwdSquare-FTTLDRC-DRTungsten-GEOMETRYTREE-TRACKEVALHITS_MBMC/G4EICDetector_eventtree.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'geometry.root'","NEWINPUT_GEOMETRY",true,true,false,  false,true,true,-1,0,false,0)'
fi

if [ $1 = "newinput_pythia" ]; then
    #input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-ETTL-CTTL-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-FTTLS3LC-GEOMETRYTREE-TEST/G4EICDetector_eventtree.root
    #geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-ETTL-CTTL-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-FTTLS3LC-GEOMETRYTREE-TEST/geometry.root
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-GEOMETRYTREE-PYTHIA-1x/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-ASYM-GEOMETRYTREE-PYTHIA-1x/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","NEWINPUT_TMSTUD_PYTHIA",true,false  ,true,true,-1,0,false,0)'
fi
#valgrind --tool=callgrind 

if [ $1 = "newdet_pythia" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-1/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-1/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","NEWDETS_PYTHIA",true,false  ,true,true,-1,0,false,0)'
fi


if [ $1 = "newdet_pythia_fix" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-EVALEIC-lll-tst2/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-PYTHIA-EVALEIC-lll-tst2/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","PYTHIA_EVALPROJFIX",true,false  ,true,true,-1,0,false,0)'
fi


if [ $1 = "newinput_singlep" ]; then
    input=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-SINGLEPART-EVALEIC-1/G4EICDetector_eventtree.root
    geometry=/media/nschmidt/local/EIC_running/ModularDetector/Modular_ALLSILICON-TTLF-INNERTRACKING-GEOMETRYTREE-TRACKEVALHITS-LFHCAL-GEOMETRYTREE-EEMCH-BECAL-SINGLEPART-EVALEIC-1/geometry.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","NEWINPUT_TMSTUD_SINGLEP_NEW",true,false  ,true,true,-1,0,false,0)'
fi




if [ $1 = "centralsim_singlepion" ]; then
maxevt=10000
    input=/media/nschmidt/local/EIC_running/singlePion/eval_00000/DST_General_particleGun_singlePion_merged_eval.root
    geometry=/media/nschmidt/local/EIC_running/singlePion/eval_00000/geometry2ndCampaign.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","CENTRALSIM_SINGLEPION",true,false  ,true,true,'$maxevt',0,false,0)'
fi

if [ $1 = "centralsim_singleelectron" ]; then
maxevt=10000
    input=/media/nschmidt/local/EIC_running/singleElectron/eval_00000/DST_General_particleGun_singleElectron_merged_eval.root
    geometry=/media/nschmidt/local/EIC_running/singlePion/eval_00000/geometry2ndCampaign.root
    root -x -l -b -q 'treeProcessing.C("'$input'","'$geometry'","CENTRALSIM_SINGLEELECTRON",true,false  ,true,true,'$maxevt',0,false,0)'
fi


# 7 10 15 20 25 30 50 75 100
basepath=/media/nschmidt/local/EIC_running/FoCal_TB
#filled="-Filled"
#filled="-Layered"
filled=""
if [ $1 = "focal_studies" ]; then
    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-FCGSI-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCGSI",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-FCGSI-SimpleElectron
    root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCGSI-Electron",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-FCRotated-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCRotated",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-FCTungsten-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCTungsten",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE-SimpleElectron
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-Electron",true,false  ,true,true,-1,0,false,0)'
#fi
#if [ $1 = "focal_studies_fulldet" ]; then
    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-FCFullGeo-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCFullGeo",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-FCFullGeo-FCTungsten-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCFullGeo-FCTungsten",true,false  ,true,true,-1,0,false,0)'


    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2$filled-FCFullGeo-FCTungsten-Lead-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCFullGeo-FCTungsten-Lead",true,false  ,true,true,-1,0,false,0)'
#fi

fi

#filled="-Filled"
filled="-Layered"
#filled=""
if [ $1 = "focal_studies_new" ]; then
    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-Layered-FCFullGeo-FCTungsten-Lead-SimplePionfwd-NEW
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCFullGeo-FCTungsten-Lead-NEW",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-Layered-FCFullGeo-FCTungsten-Lead-SimplePion-500GeVonly
    root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCFullGeo-FCTungsten-Lead-NEW-500G",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-FCGSI-Rot-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL-BH2-FCGSI-Rotated",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-FCGSI-SimplePion
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL-BH2-FCGSI",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-FCGSI-Rot-SimpleElectron
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL-BH2-FCGSI-Rotated-Electron",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-FCGSI-SimpleElectron
    #root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL-BH2-FCGSI-Electron",true,false  ,true,true,-1,0,false,0)'

fi

if [ $1 = "focal_studies_new2" ]; then
filled=""
    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-FCTB-SimpleElectron-NEW
    root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCTB-SimpleElectron-NEW",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-FCTB-FCTungsten-SimplePion-NEW
    root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCTB-FCTungsten-SimplePion-NEW",true,false  ,true,true,-1,0,false,0)'

    input=FoCal_TB_FOCAL-STANDALONE-GEOMETRYTREE-BLACKHOLE2-FCTB-SimplePion-NEW
    root -x -l -b -q 'treeProcessing.C("'$basepath/$input/G4EICDetector_eventtree.root'","'$basepath/$input/geometry.root'","FOCAL'$filled'-BH2-FCTB-SimplePion-NEW",true,false  ,true,true,-1,0,false,0)'

fi
