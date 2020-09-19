# Name:     01_cna_pca.r
# Updated:  20200826
# Author:   Andrew Kleist

# packages, working directory
library(here) # sets wd
library(tidyverse)

##### FUNCTIONS ################################################################

  # FUNCTION LoadAtlas: 
  
  LoadAtlasInter <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                          LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    # load file
    pca.object <- read.table(PCA.OUTPUT, sep="\t", header=TRUE)
    pca.object <- pca.object %>% filter(Chain1 == "A" & Chain2 =="B")
    
    # substitute for BW and/or CCN names
    object.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.OBJ.NAME][match(pca.object[ ,"ResNum1"], lookup.bwccn[ ,LOOK.OBJ.RES])])
    colnames(object.comm.name) <- c("source_gnccn")
    ligand.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.LIG.NAME][match(pca.object[ ,"ResNum2"], lookup.bwccn[ , LOOK.LIG.RES])])
    colnames(ligand.comm.name) <- c("target_gnccn")
    pca.object <- cbind(pca.object, object.comm.name)
    pca.object <- cbind(pca.object, ligand.comm.name)
    
    # add file name and PDB "class" (eg ck_ckr, cc_dimer, cxc_dimer)
    pdb.des <- as.data.frame(rep(PDBOBJ, nrow(pca.object)))
    colnames(pdb.des) <- c("file")
    pca.object <- cbind(pca.object, pdb.des)
    class.des <- as.data.frame(rep(CLASS, nrow(pca.object)))
    colnames(class.des) <- c("class")
    pca.object <- cbind(pca.object, class.des)
    
    # return final df, remove other objects
    return(pca.object)
    rm(pca.object, object.comm.name, ligand.comm.name, pca.object, pdb.des, class.des)
  }
  
  LoadAtlasCK <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                             LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    # load file
    pca.object <- read.table(PCA.OUTPUT, sep="\t", header=TRUE)
    pca.object <- pca.object %>% filter(Chain1 == "A" & Chain2 =="A")
    
    # substitute for BW and/or CCN names
    object.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.OBJ.NAME][match(pca.object[ ,"ResNum1"], lookup.bwccn[ ,LOOK.OBJ.RES])])
    colnames(object.comm.name) <- c("source_gnccn")
    ligand.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.OBJ.NAME][match(pca.object[ ,"ResNum2"], lookup.bwccn[ , LOOK.OBJ.RES])])
    colnames(ligand.comm.name) <- c("target_gnccn")
    pca.object <- cbind(pca.object, object.comm.name)
    pca.object <- cbind(pca.object, ligand.comm.name)
    
    # add file name and PDB "class" (eg ck_ckr, cc_dimer, cxc_dimer)
    pdb.des <- as.data.frame(rep(PDBOBJ, nrow(pca.object)))
    colnames(pdb.des) <- c("file")
    pca.object <- cbind(pca.object, pdb.des)
    class.des <- as.data.frame(rep(CLASS, nrow(pca.object)))
    colnames(class.des) <- c("class")
    pca.object <- cbind(pca.object, class.des)
    
    # return final df, remove other objects
    return(pca.object)
    rm(pca.object, object.comm.name, ligand.comm.name, pca.object, pdb.des, class.des)
  }
  
  
  LoadAtlasCKR <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                             LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    # load file
    pca.object <- read.table(PCA.OUTPUT, sep="\t", header=TRUE)
    pca.object <- pca.object %>% filter(Chain1 == "B" & Chain2 =="B")
    
    # substitute for BW and/or CCN names
    object.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.LIG.NAME][match(pca.object[ ,"ResNum1"], lookup.bwccn[ ,LOOK.LIG.RES])])
    colnames(object.comm.name) <- c("source_gnccn")
    ligand.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.LIG.NAME][match(pca.object[ ,"ResNum2"], lookup.bwccn[ , LOOK.LIG.RES])])
    colnames(ligand.comm.name) <- c("target_gnccn")
    pca.object <- cbind(pca.object, object.comm.name)
    pca.object <- cbind(pca.object, ligand.comm.name)
    
    # add file name and PDB "class" (eg ck_ckr, cc_dimer, cxc_dimer)
    pdb.des <- as.data.frame(rep(PDBOBJ, nrow(pca.object)))
    colnames(pdb.des) <- c("file")
    pca.object <- cbind(pca.object, pdb.des)
    class.des <- as.data.frame(rep(CLASS, nrow(pca.object)))
    colnames(class.des) <- c("class")
    pca.object <- cbind(pca.object, class.des)
    
    # return final df, remove other objects
    return(pca.object)
    rm(pca.object, object.comm.name, ligand.comm.name, pca.object, pdb.des, class.des)
  }
  
  LoadAtlas <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                        LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    
    inter <- LoadAtlasInter(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                            LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
    ck <- LoadAtlasCK(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                            LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
    ckr <- LoadAtlasCKR(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                            LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
    master <- rbind(inter, ck, ckr)
  }
  
##### 1: CALCULATE INTERMOLECULAR CONTACTS FROM CK-CKR COMPLEXES ###############
  
  # load CCN and BW conversion file
  lookup.bwccn <- read.csv("input/lookup_ckckr_domain_20200827_unique.csv")
  
  # (1) CXCL8:CXCR1
  cna.1ilp <- LoadAtlas("input/1ILPCKCLEAN1538673510.txt",
                        "clean_1ilp_ck", "ccn_1ilp_ck",     # object resno, ccn
                        "clean_1ilp_ckr", "bw_1ilp_ckr",    # ligand resno, bw
                        "1ilp",                             # PDB ID
                        "nmr")                              # PDB "class"

  
  # (2) CCL11:CCR3
  cna.2mpm <- LoadAtlas("input/2MPMCKCLEAN1538684712.txt",
                          "clean_2mpm_ck", "ccn_2mpm_ck",   # object resno, ccn
                          "clean_2mpm_ckr", "bw_2mpm_ckr",  # ligand resno, bw
                          "2mpm",                           # PDB ID
                          "nmr")                            # PDB "class"
  
  # (3) CXCL12:CXCR4
  cna.2n55 <- LoadAtlas("input/2N55CKCLEAN1538684756.txt",
                          "clean_2n55_ck", "ccn_2n55_ck",     # object resno, ccn
                          "clean_2n55_ckr", "bw_2n55_ckr",  # ligand resno, bw
                          "2n55",                           # PDB ID
                          "nmr")                            # PDB "class"
  
  # (4) vMIPII:CXCR4
  cna.4rws <- LoadAtlas("input/4RWSCKCLEAN1538684823.txt",
                          "clean_4rws_ck", "ccn_4rws_ck",     # object resno, ccn
                          "clean_4rws_ckr", "bw_4rws_ckr",  # ligand resno, bw
                          "4rws",                           # PDB ID
                          "xray")                           # PDB "class"
  
  # (5) CX3CL1:US28
  cna.4xt1 <- LoadAtlas("input/4XT1CKCLEAN1538684872.txt",
                          "clean_4xt1_ck", "ccn_4xt1_ck",     # object resno, ccn
                          "clean_4xt1_ckr", "bw_4xt1_ckr",  # ligand resno, bw
                          "4xt1",                           # PDB ID
                          "xray")                           # PDB "class"
  
  # (6) CCL5[5P7]:CCR5
  cna.5uiw <- LoadAtlas("input/5UIWCKCLEAN1538684915.txt",
                          "clean_5uiw_ck", "ccn_5uiw_ck",     # object resno, ccn
                          "clean_5uiw_ckr", "bw_5uiw_ckr",  # ligand resno, bw
                          "5uiw",                           # PDB ID
                          "xray")                           # PDB "class"

  # (7) CXCL12:CXCR4
  cna.2k05 <- LoadAtlas("input/2K05CLEAN1538684566.txt",
                          "clean_2k05_ck", "ccn_2k05_ck",     # object resno, ccn
                          "clean_2k05_ckr", "bw_2k05_ckr",  # ligand resno, bw
                          "2k05",                           # PDB ID
                          "nmr")                            # PDB "class"
  
  # (8) CCL5:CCR5
  cna.6fgp <- LoadAtlas("input/6FGPCLEAN1538688371.txt",
                        "clean_6fgp_ck", "ccn_6fgp_ck",     # object resno, ccn
                        "clean_6fgp_ckr", "bw_6fgp_ckr",    # ligand resno, bw
                        "6fgp",                             # PDB ID
                        "nmr")                              # PDB "class"
  
  # (9) CX3CL1.35:US28
  cna.5wb2 <- LoadAtlas("input/5WB2CLEAN1538684959.txt",
                        "clean_5wb2_ck", "ccn_5wb2_ck",     # object resno, ccn
                        "clean_5wb2_ckr", "bw_5wb2_ckr",    # ligand resno, bw
                        "5wb2",                             # PDB ID
                        "xray")                             # PDB "class"
  
  # (10) CCL20:CCR6
  cna.6wwz <- LoadAtlas("input/6WWZCLEAN1598531253.txt",
                        "clean_6wwz_ck", "ccn_6wwz_ck",     # object resno, ccn
                        "clean_6wwz_ckr", "bw_6wwz_ckr",    # ligand resno, bw
                        "6wwz",                             # PDB ID
                        "xray")                             # PDB "class"
  
  # (11) CXCL8:CXCR2
  cna.6lfo <- LoadAtlas("input/6LFOCLEAN1598542701.txt",
                        "clean_6lfo_ck", "ccn_6lfo_ck",     # object resno, ccn
                        "clean_6lfo_ckr", "bw_6lfo_ckr",    # ligand resno, bw
                        "6lfo",                             # PDB ID
                        "xray")                             # PDB "class"
  
  # bind all contacts into a single df (9 PDBs), then remove
  cna.master <- rbind(cna.1ilp, cna.6fgp, cna.2mpm, cna.2n55, cna.2k05, 
                      cna.5uiw, cna.4rws, cna.4xt1, cna.5wb2, cna.6lfo, cna.6wwz)

  rm(cna.1ilp, cna.6fgp, cna.2mpm, cna.2n55, cna.2k05, 
                      cna.5uiw, cna.4rws, cna.4xt1, cna.5wb2, cna.6lfo, cna.6wwz)
  
  # add domain designations
  lookup.domain <- read.csv("input/bwccn_domain_conversion.csv")
  cna.master$dom1 <- lookup.domain$dom[match(unlist(cna.master[ ,"source_gnccn"]), lookup.domain$bwccn)]
  cna.master$dom2 <- lookup.domain$dom[match(unlist(cna.master[ ,"target_gnccn"]), lookup.domain$bwccn)]  
  
  # COMMENT OUT TO PRESERVE UNIQUE ATOMIC CONTACTS
  # make unique such that each residue can only make one unique contact
  cna.master <- cna.master[, -c(11:13)]
  cna.master <- cna.master[!duplicated(cna.master[c(7,9,13)]),]
  
  # write df
  write_csv(cna.master, "output/RIN.csv")
  rm(lookup.domain, lookup.bwccn)
