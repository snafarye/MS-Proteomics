library(readxl)
PatientInfo=read_excel("~/Downloads/PatientInfo.xlsx")
PatientInfo$AgeCategory=as.factor(PatientInfo$AgeCategory)
PatientInfo$LenType=as.factor(PatientInfo$LenType)
PatientInfo$"Sex 0=female 1=male"[PatientInfo$"Sex 0=female 1=male"==0]="female"
PatientInfo$"Sex 0=female 1=male"[PatientInfo$"Sex 0=female 1=male"==1]="male"
names(PatientInfo)[4]="Sex"
PatientInfo$Sex=as.factor(PatientInfo$Sex)
names(PatientInfo)[6]="Schirmer_collection_(mm)"
PatientInfo$"Previous red eye 0=no, 1=yes"[PatientInfo$"Previous red eye 0=no, 1=yes"==0]="no"
PatientInfo$"Previous red eye 0=no, 1=yes"[PatientInfo$"Previous red eye 0=no, 1=yes"==1]="yes"
names(PatientInfo)[7]="Previous_red_eye"
PatientInfo$Previous_red_eye=as.factor(PatientInfo$Previous_red_eye)
PatientInfo$Race=as.factor(PatientInfo$Race)
PatientInfo$"ethniciy 0=hisp, 1=nonhisp"[PatientInfo$"ethniciy 0=hisp, 1=nonhisp"==0]="hisp"
PatientInfo$"ethniciy 0=hisp, 1=nonhisp"[PatientInfo$"ethniciy 0=hisp, 1=nonhisp"==1]="nonhisp"
names(PatientInfo)[9]="Ethnicity"
names(PatientInfo)[10]="Medical_history"
PatientInfo$Ethnicity=as.factor(PatientInfo$Ethnicity)
PatientInfo$Brand=as.factor(PatientInfo$Brand)
PatientInfo$Solutions=as.factor(PatientInfo$Solutions)
names(PatientInfo)[14]="Duration_of_use_(yrs)"
names(PatientInfo)[15]="Lid_redness"
names(PatientInfo)[16]="Lid_Roughness"
PatientInfo$"Papillae 0=no, 1=yes"[PatientInfo$"Papillae 0=no, 1=yes"==0]="no"
PatientInfo$"Papillae 0=no, 1=yes"[PatientInfo$"Papillae 0=no, 1=yes"==1]="yes"
names(PatientInfo)[17]="Papillae"
PatientInfo$Papillae=as.factor(PatientInfo$Papillae)
PatientInfo$"Follicles 0=no, 1=yes"[PatientInfo$"Follicles 0=no, 1=yes"==0]="no"
PatientInfo$"Follicles 0=no, 1=yes"[PatientInfo$"Follicles 0=no, 1=yes"==1]="yes"
names(PatientInfo)[18]="Follicles"
PatientInfo$Follicles=as.factor(PatientInfo$Follicles)
names(PatientInfo)[19]="Limbal_Redness"
names(PatientInfo)[20]="Bulbar_redness"
PatientInfo$"Corneal Vascularization 0=absent, 1=present"[PatientInfo$"Corneal Vascularization 0=absent, 1=present"==0]="absent"
PatientInfo$"Corneal Vascularization 0=absent, 1=present"[PatientInfo$"Corneal Vascularization 0=absent, 1=present"==1]="present"
names(PatientInfo)[21]="Corneal_Vascularization"
PatientInfo$Corneal_Vascularization=as.factor(PatientInfo$Corneal_Vascularization)
names(PatientInfo)[23]="Corneal_staining"
names(PatientInfo)[24]="Conj_staining"
names(PatientInfo)[25]="TBUT"
FolcPat=PatientInfo[,-c(3,18)]
summary(lm(Corneal_staining~AgeCategory,data=FolcPat))

