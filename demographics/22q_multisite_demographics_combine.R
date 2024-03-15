
rm(list = ls(all.names = TRUE))

# list of packages to load
#packages <- c("mgcv", "gratia", "longCombat", "lme4","lmerTest", "invgamma", "ciftiTools", "dplyr", "magrittr", "DescTools", "parallel", "tableone", "viridis", "ggplot2", "tibble", "scico", "reshape2", "data.table","readxl","gridExtra","ggpubr")
packages <- c("magrittr", "dplyr", "tidyr", "readxl","tableone","ggplot2")

# Install packages not yet installed
# Note: ciftiTools install fails if R is started without enough memory
installed_packages <- packages %in% rownames(installed.packages())
# comment out line below to skip install
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# path to demographics dir
demopath <- "/Users/charlie/Dropbox/github/22q_chr_fmri/demographics/"

# list of sessions with bold (from qunex preprocess scripts)
# NOTE: trio sessions on hoffman but marked as excluded in demo_mri:  "Q_0040_05142010" "Q_0040_06272011" "Q_0126_08012011"
# need to investigate and exclude on hoffman
# NOTE: prisma sessions on hoffman but not in demo_mri."Q_0432_02142020" "Q_0461_09112021" 

#SUNY
#sessions_suny=c("X002","X004","X005","X008","X010","X014","X018","X019","X020","X023","X024","X027","X028","X029","X030","X031","X032","X034","X041","X042","X043","X047","X061","X064","X066","X068","X073","X078","X084","X086","X108","X109","X111","X117","X119","X123","X124","X125","X126","X128","X132","X135","X137","X138","X139","X143","X145","X146","X149","X152","X155","X157","X166","X167","X170","X174","X176","X183","X189","X191","X192","X194","X196","X198","X203","X204","X205","X206","X214","X215","X218","X225","X228","X231","X232","X233","X236")
#sessions_suny=c("X002", "X004", "X005", "X008", "X010", "X014", "X018", "X019", "X020", "X023", "X024", "X027", "X028", "X029", "X030", "X031", "X032", "X034", "X041", "X042", "X043", "X047", "X061", "X064", "X066", "X068", "X073", "X078", "X084", "X086", "X108", "X109", "X111", "X117", "X119", "X123", "X124", "X125", "X126", "X128", "X132", "X135", "X137", "X138", "X139", "X143", "X145", "X146", "X149", "X152", "X155", "X157", "X166", "X167", "X170", "X174", "X176", "X183", "X189", "X191", "X192", "X194", "X196", "X198", "X203", "X204", "X205", "X206", "X214", "X215", "X218", "X225", "X228", "X231", "X232", "X233", "X236")
sessions_suny=c("X002", "X004", "X005", "X008", "X010", "X014", "X018", "X019", "X020", "X023", "X024", "X027", "X028", "X029", "X030", "X031", "X032", "X034", "X041", "X042", "X043", "X047", "X061", "X064", "X066", "X068", "X073", "X078", "X084", "X086", "X108", "X109", "X111", "X117", "X119", "X123", "X124", "X125", "X126", "X128", "X132", "X135", "X137", "X138", "X139", "X143", "X145", "X146", "X149", "X152", "X155", "X157", "X166", "X167", "X170", "X174", "X176", "X183", "X189", "X191", "X192", "X194", "X196", "X198", "X203", "X204", "X205", "X206", "X214", "X215", "X218", "X225", "X228", "X231", "X232", "X233", "X236")

#IoP
#sessions_iop=c("GQAIMS01","GQAIMS02","GQAIMS03","GQAIMS04","GQAIMS05","GQAIMS06","GQAIMS07","GQAIMS08","GQAIMS09","GQAIMS11","GQAIMS12","GQAIMS13","GQAIMS15","GQAIMS16","GQAIMS17","GQAIMS18","GQAIMS20","GQAIMS21","GQAIMS22","GQAIMS23","GQAIMS24","GQAIMS28","GQAIMS29","GQAIMS30","GQAIMS31","GQAIMS33","GQAIMS34","GQAIMS35","GQAIMS37","GQAIMS39","GQAIMS40","GQAIMS41","GQAIMS42","GQAIMS43","GQAIMS47","GQAIMS49","GQAIMS50","GQAIMS51","GQAIMS52","GQAIMS53","GQAIMS55","GQAIMS58","GQAIMS59","GQAIMS60","GQAIMS61","GQAIMS62","GQAIMS63","GQAIMS64","GQAIMS65","GQAIMS66","GQAIMS67","GQAIMS69","GQAIMS72","GQAIMS73","GQAIMS76","GQAIMS79","GQAIMS80","GQAIMS82","GQAIMS83","GQAIMS84","GQAIMS85","GQAIMS86","GQAIMS87","GQAIMS88","GQAIMS89","GQAIMS90","GQAIMS91","GQAIMS92","GQAIMS93","GQAIMS94")
#sessions_iop=c("GQAIMS01", "GQAIMS02", "GQAIMS03", "GQAIMS04", "GQAIMS05", "GQAIMS06", "GQAIMS07", "GQAIMS08", "GQAIMS09", "GQAIMS11", "GQAIMS12", "GQAIMS13", "GQAIMS15", "GQAIMS16", "GQAIMS17", "GQAIMS20", "GQAIMS21", "GQAIMS22", "GQAIMS23", "GQAIMS24", "GQAIMS28", "GQAIMS29", "GQAIMS30", "GQAIMS31", "GQAIMS33", "GQAIMS34", "GQAIMS35", "GQAIMS37", "GQAIMS39", "GQAIMS41", "GQAIMS42", "GQAIMS43", "GQAIMS47", "GQAIMS49", "GQAIMS50", "GQAIMS51", "GQAIMS52", "GQAIMS53", "GQAIMS58", "GQAIMS59", "GQAIMS61", "GQAIMS62", "GQAIMS63", "GQAIMS64", "GQAIMS65", "GQAIMS66", "GQAIMS67", "GQAIMS69", "GQAIMS72", "GQAIMS73", "GQAIMS76", "GQAIMS79", "GQAIMS80", "GQAIMS82", "GQAIMS84", "GQAIMS85", "GQAIMS86", "GQAIMS87", "GQAIMS88", "GQAIMS89", "GQAIMS90", "GQAIMS91", "GQAIMS92", "GQAIMS93", "GQAIMS94")
sessions_iop=c("GQAIMS01", "GQAIMS02", "GQAIMS03", "GQAIMS04", "GQAIMS05", "GQAIMS06", "GQAIMS07", "GQAIMS08", "GQAIMS09", "GQAIMS11", "GQAIMS12", "GQAIMS13", "GQAIMS15", "GQAIMS16", "GQAIMS17", "GQAIMS20", "GQAIMS21", "GQAIMS22", "GQAIMS23", "GQAIMS24", "GQAIMS28", "GQAIMS29", "GQAIMS30", "GQAIMS31", "GQAIMS33", "GQAIMS34", "GQAIMS35", "GQAIMS37", "GQAIMS39", "GQAIMS41", "GQAIMS42", "GQAIMS43", "GQAIMS47", "GQAIMS49", "GQAIMS50", "GQAIMS51", "GQAIMS52", "GQAIMS53", "GQAIMS58", "GQAIMS59", "GQAIMS61", "GQAIMS62", "GQAIMS63", "GQAIMS64", "GQAIMS65", "GQAIMS66", "GQAIMS67", "GQAIMS69", "GQAIMS72", "GQAIMS73", "GQAIMS76", "GQAIMS79", "GQAIMS80", "GQAIMS82", "GQAIMS84", "GQAIMS85", "GQAIMS86", "GQAIMS87", "GQAIMS88", "GQAIMS89", "GQAIMS90", "GQAIMS91", "GQAIMS92", "GQAIMS93", "GQAIMS94")

#Rome
#sessions_rome=c("C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11","C12","C13","C14","C15","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27","D01","D02","D03","D04","D05","D06","D07","D09","D10","D11","D12","D14","D15","D16","D17","D18","D20","D21","D22","D23","D24","P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23","P24","P25","P26","P27","P28","P29","P30","P31","P32","P33","P35","P36","P37","P38")
#sessions_rome=c("C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14", "C15", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D09", "D10", "D11", "D12", "D14", "D15", "D16", "D18", "D20", "D21", "D22", "D23", "P02", "P04", "P05", "P06", "P07", "P08", "P09", "P11", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20", "P21", "P22", "P23", "P24", "P25", "P26", "P27", "P28", "P29", "P30", "P31", "P32", "P33", "P35", "P36", "P37", "P38")
sessions_rome=c("C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14", "C15", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D09", "D10", "D11", "D12", "D14", "D15", "D16", "D18", "D20", "D21", "D22", "D23")
#"P02", "P04", "P05", "P06", "P07", "P08", "P09", "P11", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20", "P21", "P22", "P23", "P24", "P25", "P26", "P27", "P28", "P29", "P30", "P31", "P32", "P33", "P35", "P36", "P37", "P38")

# UCLA_trio
#sessions_trio=c("Q_0001_09242012","Q_0001_10152010","Q_0005_04182011","Q_0007_04012011","Q_0009_09292011","Q_0016_03292011","Q_0016_10082012","Q_0019_05052011","Q_0020_03262012","Q_0020_05052011","Q_0021_01132011","Q_0021_11052012","Q_0024_04152011","Q_0024_09212012","Q_0026_03232012","Q_0026_05032013","Q_0033_03192010","Q_0036_03302010","Q_0036_06212012","Q_0036_06282011","Q_0037_04202011","Q_0038_06042010","Q_0039_04132010","Q_0039_06012011","Q_0040_05142010","Q_0040_06272011","Q_0041_09262012","Q_0045_07132010","Q_0045_08302011","Q_0045_10052012","Q_0051_01302012","Q_0052_07302010","Q_0053_08022010","Q_0053_12212012","Q_0054_08062010","Q_0056_08252010","Q_0056_09102012","Q_0059_08032011","Q_0062_11302011","Q_0076_10212010","Q_0077_02292012","Q_0077_03112013","Q_0077_10252010","Q_0078_11102011","Q_0080_11222010","Q_0081_11242010","Q_0082_01062012","Q_0082_01112013","Q_0082_12162010","Q_0085_02242012","Q_0091_03052012","Q_0091_03072011","Q_0092_03102011","Q_0093_03082011","Q_0093_03122012","Q_0093_03122013","Q_0098_04012011","Q_0099_04012011","Q_0100_04082011","Q_0101_01032013","Q_0101_04122011","Q_0102_01032013","Q_0102_04122011","Q_0103_04182011","Q_0103_08272012","Q_0105_05022011","Q_0109_05262011","Q_0112_11142012","Q_0114_02062014","Q_0114_11092012","Q_0117_07062011","Q_0117_08312012","Q_0118_07072011","Q_0124_07292011","Q_0124_09172012","Q_0126_08012011","Q_0127_03052013","Q_0127_08012014","Q_0127_08082011","Q_0130_08132012","Q_0130_08162011","Q_0130_11042014","Q_0132_08222011","Q_0135_09012011","Q_0136_09012011","Q_0137_09022011","Q_0138_09032013","Q_0138_09042012","Q_0138_09062011","Q_0141_10072013","Q_0141_10092012","Q_0146_11082012","Q_0146_12072011","Q_0147_08112015","Q_0147_08142014","Q_0147_12122011","Q_0149_11192012","Q_0149_12192011","Q_0150_11192012","Q_0150_12202011","Q_0151_01122012","Q_0153_01202012","Q_0156_01312012","Q_0156_12102013","Q_0157_01312012","Q_0159_02122013","Q_0161_03012012","Q_0161_04192013","Q_0162_03202012","Q_0163_03262012","Q_0166_01092014","Q_0166_04052012","Q_0168_04052012","Q_0169_04062012","Q_0170_04102012","Q_0170_11252014","Q_0171_04102012","Q_0171_11252014","Q_0172_04122012","Q_0172_04262013","Q_0173_04122012","Q_0173_04172014","Q_0173_04262013","Q_0174_04182012","Q_0174_04272015","Q_0176_04242012","Q_0176_05072013","Q_0177_05012012","Q_0178_05012012","Q_0182_05182012","Q_0184_05242012","Q_0185_05242012","Q_0186_04092015","Q_0188_05242012","Q_0189_05242012","Q_0190_06012012","Q_0190_08132013","Q_0192_06112013","Q_0196_06192015","Q_0196_06202012","Q_0196_09122013","Q_0200_04032015","Q_0200_07092012","Q_0200_08052013","Q_0213_08032012","Q_0215_08062012","Q_0215_08172015","Q_0216_08062012","Q_0216_08172015","Q_0217_01142014","Q_0217_02042015","Q_0217_08092012","Q_0219_08122014","Q_0219_08162012","Q_0219_09012015","Q_0222_08212012","Q_0223_08212012","Q_0227_10032012","Q_0228_01292014","Q_0228_09272012","Q_0229_10022012","Q_0232_12122012","Q_0234_03252014","Q_0234_12132012","Q_0236_03072013","Q_0238_03202013","Q_0238_08042014","Q_0238_11132015","Q_0240_03232015","Q_0240_03272013","Q_0240_06262014","Q_0242_03272013","Q_0242_05162016","Q_0244_09162013","Q_0244_09232014","Q_0246_09242013","Q_0252_12102013","Q_0255_10212014","Q_0260_06092014","Q_0260_10262015","Q_0266_08122014","Q_0268_10132014","Q_0269_10142014","Q_0277_12082014","Q_0284_03182015","Q_0285_03182015","Q_0286_03182015","Q_0287_03182015","Q_0291_04172015","Q_0307_10132015","Q_0310_11182015","Q_0311_12082015","Q_0315_01262016","Q_0321_03232016","Q_0322_03312016","Q_0326_04142016","Q_0327_04142016","Q_0328_04262016")
#sessions_trio=c("Q_0001_09242012", "Q_0001_10152010", "Q_0005_04182011", "Q_0007_04012011", "Q_0009_09292011", "Q_0016_03292011", "Q_0016_10082012", "Q_0020_03262012", "Q_0020_05052011", "Q_0021_01132011", "Q_0021_11052012", "Q_0024_04152011", "Q_0024_09212012", "Q_0026_03232012", "Q_0026_05032013", "Q_0033_03192010", "Q_0036_03302010", "Q_0036_06212012", "Q_0036_06282011", "Q_0037_04202011", "Q_0038_06042010", "Q_0039_04132010", "Q_0039_06012011", "Q_0040_05142010", "Q_0040_06272011", "Q_0041_09262012", "Q_0045_07132010", "Q_0045_08302011", "Q_0045_10052012", "Q_0051_01302012", "Q_0052_07302010", "Q_0053_08022010", "Q_0053_12212012", "Q_0054_08062010", "Q_0056_08252010", "Q_0056_09102012", "Q_0059_08032011", "Q_0062_11302011", "Q_0076_10212010", "Q_0077_02292012", "Q_0077_03112013", "Q_0077_10252010", "Q_0078_11102011", "Q_0080_11222010", "Q_0081_11242010", "Q_0082_01062012", "Q_0082_01112013", "Q_0082_12162010", "Q_0085_02242012", "Q_0091_03052012", "Q_0091_03072011", "Q_0092_03102011", "Q_0093_03082011", "Q_0093_03122012", "Q_0093_03122013", "Q_0098_04012011", "Q_0099_04012011", "Q_0100_04082011", "Q_0101_01032013", "Q_0101_04122011", "Q_0102_01032013", "Q_0102_04122011", "Q_0103_04182011", "Q_0103_08272012", "Q_0105_05022011", "Q_0109_05262011", "Q_0112_11142012", "Q_0114_11092012", "Q_0117_07062011", "Q_0117_08312012", "Q_0118_07072011", "Q_0124_07292011", "Q_0124_09172012", "Q_0126_08012011", "Q_0127_03052013", "Q_0127_08012014", "Q_0127_08082011", "Q_0130_08132012", "Q_0130_08162011", "Q_0130_11042014", "Q_0132_08222011", "Q_0135_09012011", "Q_0141_10072013", "Q_0141_10092012", "Q_0146_11082012", "Q_0146_12072011", "Q_0147_08112015", "Q_0147_08142014", "Q_0147_12122011", "Q_0149_11192012", "Q_0149_12192011", "Q_0150_11192012", "Q_0150_12202011", "Q_0151_01122012", "Q_0153_01202012", "Q_0156_01312012", "Q_0156_12102013", "Q_0157_01312012", "Q_0159_02122013", "Q_0161_03012012", "Q_0161_04192013", "Q_0162_03202012", "Q_0163_03262012", "Q_0166_01092014", "Q_0166_04052012", "Q_0168_04052012", "Q_0169_04062012", "Q_0170_04102012", "Q_0170_11252014", "Q_0171_04102012", "Q_0171_11252014", "Q_0172_04122012", "Q_0172_04262013", "Q_0173_04122012", "Q_0173_04172014", "Q_0173_04262013", "Q_0174_04182012", "Q_0174_04272015", "Q_0176_04242012", "Q_0176_05072013", "Q_0177_05012012", "Q_0178_05012012", "Q_0182_05182012", "Q_0184_05242012", "Q_0185_05242012", "Q_0186_04092015", "Q_0188_05242012", "Q_0189_05242012", "Q_0190_06012012", "Q_0190_08132013", "Q_0196_06192015", "Q_0196_06202012", "Q_0196_09122013", "Q_0200_04032015", "Q_0200_08052013", "Q_0213_08032012", "Q_0215_08062012", "Q_0215_08172015", "Q_0216_08062012", "Q_0216_08172015", "Q_0217_01142014", "Q_0217_02042015", "Q_0217_08092012", "Q_0219_08122014", "Q_0219_08162012", "Q_0219_09012015", "Q_0222_08212012", "Q_0223_08212012", "Q_0227_10032012", "Q_0228_01292014", "Q_0228_09272012", "Q_0229_10022012", "Q_0232_12122012", "Q_0234_03252014", "Q_0234_12132012", "Q_0236_03072013", "Q_0238_03202013", "Q_0238_08042014", "Q_0238_11132015", "Q_0240_03232015", "Q_0240_03272013", "Q_0240_06262014", "Q_0242_03272013", "Q_0242_05162016", "Q_0244_09162013", "Q_0244_09232014", "Q_0246_09242013", "Q_0252_12102013", "Q_0255_10212014", "Q_0260_06092014", "Q_0260_10262015", "Q_0266_08122014", "Q_0268_10132014", "Q_0269_10142014", "Q_0277_12082014", "Q_0284_03182015", "Q_0285_03182015", "Q_0286_03182015", "Q_0291_04172015", "Q_0307_10132015", "Q_0310_11182015", "Q_0311_12082015", "Q_0315_01262016", "Q_0321_03232016", "Q_0322_03312016", "Q_0326_04142016", "Q_0327_04142016")
sessions_trio=c("Q_0001_09242012", "Q_0001_10152010", "Q_0005_04182011", "Q_0007_04012011", "Q_0009_09292011", "Q_0016_03292011", "Q_0016_10082012", "Q_0020_03262012", "Q_0020_05052011", "Q_0021_01132011", "Q_0021_11052012", "Q_0024_04152011", "Q_0024_09212012", "Q_0026_03232012", "Q_0026_05032013", "Q_0033_03192010", "Q_0036_03302010", "Q_0036_06212012", "Q_0036_06282011", "Q_0037_04202011", "Q_0038_06042010", "Q_0039_04132010", "Q_0039_06012011", "Q_0040_05142010", "Q_0040_06272011", "Q_0041_09262012", "Q_0045_07132010", "Q_0045_08302011", "Q_0045_10052012", "Q_0051_01302012", "Q_0052_07302010", "Q_0053_08022010", "Q_0053_12212012", "Q_0054_08062010", "Q_0056_08252010", "Q_0056_09102012", "Q_0059_08032011", "Q_0062_11302011", "Q_0076_10212010", "Q_0077_02292012", "Q_0077_03112013", "Q_0077_10252010", "Q_0078_11102011", "Q_0080_11222010", "Q_0081_11242010", "Q_0082_01062012", "Q_0082_01112013", "Q_0082_12162010", "Q_0085_02242012", "Q_0091_03052012", "Q_0091_03072011", "Q_0092_03102011", "Q_0093_03082011", "Q_0093_03122012", "Q_0093_03122013", "Q_0098_04012011", "Q_0099_04012011", "Q_0100_04082011", "Q_0101_01032013", "Q_0101_04122011", "Q_0102_01032013", "Q_0102_04122011", "Q_0103_04182011", "Q_0103_08272012", "Q_0105_05022011", "Q_0109_05262011", "Q_0112_11142012", "Q_0114_11092012", "Q_0117_07062011", "Q_0117_08312012", "Q_0118_07072011", "Q_0124_07292011", "Q_0124_09172012", "Q_0126_08012011", "Q_0127_03052013", "Q_0127_08012014", "Q_0127_08082011", "Q_0130_08132012", "Q_0130_08162011", "Q_0130_11042014", "Q_0132_08222011", "Q_0135_09012011", "Q_0141_10072013", "Q_0141_10092012", "Q_0146_11082012", "Q_0146_12072011", "Q_0147_08112015", "Q_0147_08142014", "Q_0147_12122011", "Q_0149_11192012", "Q_0149_12192011", "Q_0150_11192012", "Q_0150_12202011", "Q_0151_01122012", "Q_0153_01202012", "Q_0156_01312012", "Q_0156_12102013", "Q_0157_01312012", "Q_0159_02122013", "Q_0161_03012012", "Q_0161_04192013", "Q_0162_03202012", "Q_0163_03262012", "Q_0166_01092014", "Q_0166_04052012", "Q_0168_04052012", "Q_0169_04062012", "Q_0170_04102012", "Q_0170_11252014", "Q_0171_04102012", "Q_0171_11252014", "Q_0172_04122012", "Q_0172_04262013", "Q_0173_04122012", "Q_0173_04172014", "Q_0173_04262013", "Q_0174_04182012", "Q_0174_04272015", "Q_0176_04242012", "Q_0176_05072013", "Q_0177_05012012", "Q_0178_05012012", "Q_0182_05182012", "Q_0184_05242012", "Q_0185_05242012", "Q_0186_04092015", "Q_0188_05242012", "Q_0189_05242012", "Q_0190_06012012", "Q_0190_08132013", "Q_0196_06192015", "Q_0196_06202012", "Q_0196_09122013", "Q_0200_04032015", "Q_0200_08052013", "Q_0213_08032012", "Q_0215_08062012", "Q_0215_08172015", "Q_0216_08062012", "Q_0216_08172015", "Q_0217_01142014", "Q_0217_02042015", "Q_0217_08092012", "Q_0219_08122014", "Q_0219_08162012", "Q_0219_09012015", "Q_0222_08212012", "Q_0223_08212012", "Q_0227_10032012", "Q_0228_01292014", "Q_0228_09272012", "Q_0229_10022012", "Q_0232_12122012", "Q_0234_03252014", "Q_0234_12132012", "Q_0236_03072013", "Q_0238_03202013", "Q_0238_08042014", "Q_0238_11132015", "Q_0240_03232015", "Q_0240_03272013", "Q_0240_06262014", "Q_0242_03272013", "Q_0242_05162016", "Q_0244_09162013", "Q_0244_09232014", "Q_0246_09242013", "Q_0252_12102013", "Q_0255_10212014", "Q_0260_06092014", "Q_0260_10262015", "Q_0266_08122014", "Q_0268_10132014", "Q_0269_10142014", "Q_0277_12082014", "Q_0284_03182015", "Q_0285_03182015", "Q_0286_03182015", "Q_0291_04172015", "Q_0307_10132015", "Q_0310_11182015", "Q_0311_12082015", "Q_0315_01262016", "Q_0321_03232016", "Q_0322_03312016", "Q_0326_04142016", "Q_0327_04142016")

# UCLA_prisma
#sessions_prisma=c("Q_0246_10142016","Q_0263_11072016","Q_0271_10182016","Q_0291_11042016","Q_0005_11072017","Q_0017_10022017","Q_0019_03022020","Q_0036_08082017","Q_0036_08312018","Q_0036_09302019","Q_0037_08152017","Q_0041_10092017","Q_0051_02032017","Q_0105_03112020","Q_0105_12032018","Q_0114_12052017","Q_0141_06122018","Q_0141_08052019","Q_0147_12112017","Q_0147_12112018","Q_0196_01082020","Q_0213_05012017","Q_0214_05012017","Q_0217_01242017","Q_0217_02192020","Q_0235_05252017","Q_0238_02202018","Q_0240_09192017","Q_0240_11162018","Q_0246_10092018","Q_0246_10102017","Q_0260_06192017","Q_0260_06242019","Q_0260_06252018","Q_0263_02252019","Q_0263_02262018","Q_0277_12132016","Q_0278_11292017","Q_0278_12052019","Q_0278_12132016","Q_0279_11302017","Q_0279_12052019","Q_0279_12132016","Q_0285_03212017","Q_0285_06062018","Q_0286_03212017","Q_0287_06062018","Q_0288_06062018","Q_0289_03212017","Q_0289_06062018","Q_0291_11302018","Q_0304_12202016","Q_0310_01252018","Q_0310_02132017","Q_0310_04292019","Q_0319_03192018","Q_0319_03282017","Q_0321_03192018","Q_0321_03272017","Q_0324_04182018","Q_0326_10202017","Q_0326_12062018","Q_0327_10192017","Q_0327_12072018","Q_0331_06102019","Q_0331_06212018","Q_0331_06272017","Q_0333_04142017","Q_0334_12012016","Q_0336_01102017","Q_0345_04122017","Q_0345_08152018","Q_0345_09112019","Q_0346_04102017","Q_0346_04102018","Q_0348_04212017","Q_0348_08152018","Q_0350_04192017","Q_0350_09142018","Q_0353_04182018","Q_0353_05022017","Q_0355_05312018","Q_0356_05312018","Q_0361_10212019","Q_0361_11202018","Q_0369_04182018","Q_0369_06182019","Q_0371_08042020","Q_0374_05252018","Q_0381_08072018","Q_0381_09102019","Q_0382_08282018","Q_0383_08282018","Q_0387_08242018","Q_0387_12032019","Q_0390_09042018","Q_0390_09302019","Q_0391_09252018","Q_0395_11062018","Q_0397_10172019","Q_0402_01112019","Q_0404_03112019","Q_0407_06122019","Q_0408_05172019","Q_0414_07172019","Q_0415_07292019","Q_0416_07292019","Q_0425_11052019","Q_0429_01092020","Q_0432_02142020","Q_0433_02262020","Q_0443_06282021","Q_0446_06152021","Q_0459_08192021","Q_0461_09112021")
#sessions_prisma=c("Q_0005_11072017", "Q_0017_10022017", "Q_0019_03022020", "Q_0036_08082017", "Q_0036_08312018", "Q_0036_09302019", "Q_0037_08152017", "Q_0041_10092017", "Q_0051_02032017", "Q_0105_03112020", "Q_0105_12032018", "Q_0114_12052017", "Q_0141_06122018", "Q_0141_08052019", "Q_0147_12112017", "Q_0147_12112018", "Q_0196_01082020", "Q_0213_05012017", "Q_0214_05012017", "Q_0217_01242017", "Q_0217_02192020", "Q_0235_05252017", "Q_0238_02202018", "Q_0240_09192017", "Q_0240_11162018", "Q_0246_10092018", "Q_0246_10102017", "Q_0246_10142016", "Q_0260_06192017", "Q_0260_06242019", "Q_0260_06252018", "Q_0263_02252019", "Q_0263_02262018", "Q_0263_11072016", "Q_0271_10182016", "Q_0277_12132016", "Q_0278_11292017", "Q_0278_12052019", "Q_0278_12132016", "Q_0279_11302017", "Q_0279_12052019", "Q_0279_12132016", "Q_0285_03212017", "Q_0285_06062018", "Q_0286_03212017", "Q_0287_06062018", "Q_0288_06062018", "Q_0289_03212017", "Q_0289_06062018", "Q_0291_11042016", "Q_0291_11302018", "Q_0304_12202016", "Q_0310_01252018", "Q_0310_02132017", "Q_0310_04292019", "Q_0319_03192018", "Q_0319_03282017", "Q_0321_03192018", "Q_0321_03272017", "Q_0324_04182018", "Q_0326_10202017", "Q_0326_12062018", "Q_0327_10192017", "Q_0327_12072018", "Q_0331_06102019", "Q_0331_06212018", "Q_0331_06272017", "Q_0333_04142017", "Q_0334_12012016", "Q_0336_01102017", "Q_0345_04122017", "Q_0345_08152018", "Q_0345_09112019", "Q_0346_04102017", "Q_0346_04102018", "Q_0348_04212017", "Q_0348_08152018", "Q_0350_04192017", "Q_0350_09142018", "Q_0353_04182018", "Q_0353_05022017", "Q_0355_05312018", "Q_0356_05312018", "Q_0361_10212019", "Q_0361_11202018", "Q_0369_04182018", "Q_0369_06182019", "Q_0371_08042020", "Q_0374_05252018", "Q_0381_08072018", "Q_0381_09102019", "Q_0382_08282018", "Q_0383_08282018", "Q_0387_08242018", "Q_0387_12032019", "Q_0390_09042018", "Q_0390_09302019", "Q_0391_09252018", "Q_0395_11062018", "Q_0397_10172019", "Q_0402_01112019", "Q_0404_03112019", "Q_0407_06122019", "Q_0408_05172019", "Q_0414_07172019", "Q_0415_07292019", "Q_0416_07292019", "Q_0425_11052019", "Q_0429_01092020", "Q_0432_02142020", "Q_0433_02262020", "Q_0443_06282021", "Q_0446_06152021", "Q_0459_08192021", "Q_0461_09112021")
sessions_prisma=c("Q_0005_11072017", "Q_0017_10022017", "Q_0019_03022020", "Q_0036_08082017", "Q_0036_08312018", "Q_0036_09302019", "Q_0037_08152017", "Q_0041_10092017", "Q_0051_02032017", "Q_0105_03112020", "Q_0105_12032018", "Q_0114_12052017", "Q_0141_06122018", "Q_0141_08052019", "Q_0147_12112017", "Q_0147_12112018", "Q_0196_01082020", "Q_0213_05012017", "Q_0214_05012017", "Q_0217_01242017", "Q_0217_02192020", "Q_0235_05252017", "Q_0238_02202018", "Q_0240_09192017", "Q_0240_11162018", "Q_0246_10092018", "Q_0246_10102017", "Q_0246_10142016", "Q_0260_06192017", "Q_0260_06242019", "Q_0260_06252018", "Q_0263_02252019", "Q_0263_02262018", "Q_0263_11072016", "Q_0271_10182016", "Q_0277_12132016", "Q_0278_11292017", "Q_0278_12132016", "Q_0279_11302017", "Q_0279_12052019", "Q_0279_12132016", "Q_0285_03212017", "Q_0285_06062018", "Q_0286_03212017", "Q_0287_06062018", "Q_0288_06062018", "Q_0289_03212017", "Q_0289_06062018", "Q_0291_11042016", "Q_0291_11302018", "Q_0304_12202016", "Q_0310_01252018", "Q_0310_02132017", "Q_0310_04292019", "Q_0319_03192018", "Q_0319_03282017", "Q_0321_03192018", "Q_0321_03272017", "Q_0324_04182018", "Q_0326_10202017", "Q_0326_12062018", "Q_0327_10192017", "Q_0327_12072018", "Q_0331_06102019", "Q_0331_06212018", "Q_0331_06272017", "Q_0333_04142017", "Q_0334_12012016", "Q_0336_01102017", "Q_0345_04122017", "Q_0345_08152018", "Q_0345_09112019", "Q_0346_04102017", "Q_0346_04102018", "Q_0348_04212017", "Q_0348_08152018", "Q_0350_04192017", "Q_0350_09142018", "Q_0353_04182018", "Q_0353_05022017", "Q_0355_05312018", "Q_0356_05312018", "Q_0361_10212019", "Q_0361_11202018", "Q_0369_04182018", "Q_0369_06182019", "Q_0371_08042020", "Q_0374_05252018", "Q_0381_08072018", "Q_0381_09102019", "Q_0382_08282018", "Q_0383_08282018", "Q_0387_08242018", "Q_0387_12032019", "Q_0390_09042018", "Q_0390_09302019", "Q_0391_09252018", "Q_0395_11062018", "Q_0397_10172019", "Q_0402_01112019", "Q_0404_03112019", "Q_0407_06122019", "Q_0408_05172019", "Q_0414_07172019", "Q_0415_07292019", "Q_0416_07292019", "Q_0425_11052019", "Q_0429_01092020", "Q_0432_02142020", "Q_0433_02262020", "Q_0443_06282021", "Q_0446_06152021", "Q_0459_08192021", "Q_0461_09112021")

############################################################################################################
#UCLA


#dir <- "/Users/charlie/Dropbox/PhD/bearden_lab/22q/analyses/striatum_thalamus_fc"
#csvdir <- paste(dir,"csv",sep="/")
csvdir <- file.path(demopath,"ucla_sistat")
# get list of files in directory
files <- list.files(csvdir)
fpaths <- lapply(files, function(file) paste(csvdir,file,sep="/"))
# clean names
fnames <- gsub(".csv","",files)
fnames <- gsub("Re22Q_","",fnames)
fnames <- gsub("Form_","",fnames)
fnames <- gsub("Qry_","",fnames)
# read all, set to na: "-9999", "-9998","." 
input_all <- lapply(fpaths, read.csv, header=T, na.strings=c(".","-9999","-9998"), strip.white=T, sep=",")
names(input_all) <- fnames
df_all <- lapply(input_all, function(x) data.frame(x))

# filter based on lists above
ucla_demo <- filter(df_all$demo_mri, df_all$demo_mri$MRI_S_ID %in% c(sessions_trio,sessions_prisma))

# remove "FAMILY MEMBER" designation from subject identity
ucla_demo$SUBJECT_IDENTITY <- ucla_demo$SUBJECT_IDENTITY %>% sub("FAMILY MEMBER","",.) %>% sub(",","",.) %>% trimws(which="both") %>% as.factor
# change sex coding from 0/1 to F/M and set to factor
ucla_demo$SEX <- factor(ucla_demo$SEX,levels=c(0,1),labels=c("F","M"))
# set race=NA to 7 (unknown)
ucla_demo$RACE[is.na(ucla_demo$RACE)] <- 7
# set race as factor 1=American Indian/Alaska Native; 2=Asian; 3=Native Hawaiian/Pacific Islander; 4=Black or African American; 5=White; 6=Multiple; 7=Unknown
ucla_demo$RACE <- factor(ucla_demo$RACE,levels=c(1:7),labels=c("1_Native_American","2_Asian","3_Pacific_Island","4_Black","5_White","6_Multiple","7_Unknown"))
# ethnicity as factor with 0=N 1=Y
ucla_demo$HISPANIC[is.na(ucla_demo$HISPANIC)] <- "Unknown"
ucla_demo$HISPANIC <- factor(ucla_demo$HISPANIC,levels=c(0,1,"Unknown"),labels=c("N","Y","Unknown"))
# get more accurate age with AGEMONTH/12
ucla_demo$AGE <- as.numeric(ucla_demo$AGEMONTH)/12 

# function to add column to code timepoints relative to sample used (i.e. if visit 1 and 1.12 missing, then 1.24 is baseline)
# trio/prisma coded as T/P-visit_n where T-1 would be the subject's first trio scan and P-1 the first prisma, P-2 the second...
# function should be applied to the indicies of rows (r) in a subset of demo_mri
gettp <- function(r, df){
  sub <- df$SUBJECTID[[r]]
  visit <- df$CONVERTEDVISITNUM[[r]]
  all_visits <- df$CONVERTEDVISITNUM[which(df$SUBJECTID == sub)] %>% sort
  n_visits <- length(all_visits)
  nt_visits <-length(which(all_visits < 2))
  np_visits <- length(which(all_visits >= 2))
  visit_index <- which(all_visits == visit)
  if (visit < 2){
    label=paste("T-",visit_index,sep="")
  }else if (visit >= 2){
    p_visits <- all_visits[which(all_visits >= 2)] %>% sort
    p_visit_index <- which(p_visits == visit)
    label=paste("P-",p_visit_index,sep="")
  }
  return(c(sub,visit,label,n_visits,nt_visits,np_visits,visit_index))
}

# get timepoints
timepoints <- sapply(1:nrow(ucla_demo),function(r) gettp(r,ucla_demo)) %>% t %>% as.data.frame
colnames(timepoints) <- c("SUBJECTID","CONVERTEDVISITNUM","converted_timepoint","n_timepoints","n_trio","n_prisma","visit_index")
ucla_demo_tp <- cbind(ucla_demo,timepoints[,3:7])
ucla_demo_tp$visit_index %<>% as.factor

# add medication info from summPsych
ucla_summPsych <- df_all$summPsych
ucla_summPsych$PSYTYPE[is.na(ucla_summPsych$PSYTYPE)] <- 5
ucla_summPsych$medication <- factor(ucla_summPsych$PSYTYPE, levels=c(1,2,3,4,5), labels=c("antipsychotic","antidepressant_or_mood_stabilizer","stimulant","other","none"))
ucla_summPsych$apd_tf <- ucla_summPsych$medication == "antipsychotic"
ucla_summPsych$Med_Antipsychotic <- factor(ucla_summPsych$apd_tf, levels=c(T,F), labels=c("Y","N"))
# get psych dx 
ucla_summPsych$psych_dx <- factor(ucla_summPsych$PSYDIAGNOS, levels=c(1,0,3), labels=c("Y","N","N"))

# merge with demo
ucla_demo_use <- merge(x=ucla_demo_tp, y=ucla_summPsych[,c("SUBJECTID","CONVERTEDVISITNUM","Med_Antipsychotic","psych_dx")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T)

# get IQ
# WASI, WISC-IV, DKEFS and trail making all under df_all$DKEFS for trio data
# IQSS -- full scale WASI
ucla_neuro1 <- df_all$DKEFS[,c("SUBJECTID","CONVERTEDVISITNUM","VOCASS","MATRIXSS","IQSS")] %>% rename("WASI_verbal" = "VOCASS") %>% rename("WASI_matrix" = "MATRIXSS") %>% rename("IQ_full" = "IQSS")
# renewal neuro (prisma) under df_all$neurocogTest
ucla_neuro2 <- df_all$neurocogTest[,c("SUBJECTID","CONVERTEDVISITNUM","VOCA_TSCORE","MATRIX_TSCORE","IQ_SCORE")] %>% rename("WASI_verbal" = "VOCA_TSCORE") %>% rename("WASI_matrix" = "MATRIX_TSCORE") %>% rename("IQ_full" = "IQ_SCORE")
# combine 22q orig and renewal scores before merging with demo table
ucla_neuro <- rbind(ucla_neuro1, ucla_neuro2)
# merge neuro with demo table
ucla_demo_use <- merge(x=ucla_demo_use, y=ucla_neuro[,c("SUBJECTID","CONVERTEDVISITNUM","IQ_full")], by=c("SUBJECTID","CONVERTEDVISITNUM"), all.x=T) 
# record IQ instrument
ucla_demo_use$IQ_measure <- NA
ucla_demo_use$IQ_measure[!is.na(ucla_demo_use$IQ_full)] <- "WASI_full_scale"

# cardiac diagnoses 
pheno <- df_all$Phenotype
pheno_card <- pheno[complete.cases(pheno[,c("SUBJECTID","CONVERTEDVISITNUM", paste0("HEART",1:20))]),c("SUBJECTID","CONVERTEDVISITNUM", paste0("HEART",1:20))]
card_tf <- data.frame(SUBJECTID=unique(pheno_card$SUBJECTID) )
card_tf$cardiac <- lapply(card_tf$SUBJECTID,function(s)sum(filter(pheno_card, SUBJECTID==s)[,paste0("HEART",1:20)])>0) %>% do.call(rbind,.)

# merge cardiac info with dataframe and set NA to FALSE then convert to numeric
ucla_demo_use <- merge(x=ucla_demo_use, y=card_tf, by=c("SUBJECTID"), all.x=TRUE)
ucla_demo_use[which(is.na(ucla_demo_use$cardiac)),"cardiac"] <- FALSE
ucla_demo_use$cardiac <- factor(ucla_demo_use$cardiac, levels=c(FALSE,TRUE), labels=c("N","Y"))


# manually fix missing sex for Q_0381_09102019
ucla_demo_use[which(ucla_demo_use$MRI_S_ID == "Q_0381_09102019"),"SEX"] <- "F"

# get trio demographics
ages_trio_subset <- filter(ucla_demo_use, MRI_S_ID %in% sessions_trio)[,c("MRI_S_ID","SUBJECTID","SUBJECT_IDENTITY","AGE","SEX","RACE","HISPANIC","Med_Antipsychotic","psych_dx","IQ_full","IQ_measure","cardiac","visit_index")]
ages_trio_subset$Site <- rep("UCLAtrio",times=nrow(ages_trio_subset))

# get prisma demographics
ages_prisma_subset <- filter(ucla_demo_use, MRI_S_ID %in% sessions_prisma)[,c("MRI_S_ID","SUBJECTID","SUBJECT_IDENTITY","AGE","SEX","RACE","HISPANIC","Med_Antipsychotic","psych_dx","IQ_full","IQ_measure","cardiac","visit_index")]
ages_prisma_subset$Site <- rep("UCLAprisma",times=nrow(ages_prisma_subset))

########################################################################################################################
# SUNY
SUNY_demo <- read_xlsx(file.path(demopath,"multisite/SUNY_Demogs.xlsx"),trim_ws=T, na="", col_names=T)
# rename SUBJECT. to SUBJECTID, S[0-9] to X[0-9]
SUNY_demo$SUBJECTID <- SUNY_demo$SUBJECT. %>% gsub("S","X",.) 
# MRI_S_ID same as new subject ID
SUNY_demo$MRI_S_ID <- SUNY_demo$SUBJECTID
# rename GEND to SEX and set as factor
SUNY_demo$SEX <- as.factor(SUNY_demo$GEND)
# rename DAIG to SUBJECT_IDENTITY, changing VCFS to PATIENT-DEL and Control and Sibling to CONTROL
SUNY_demo$SUBJECT_IDENTITY <- SUNY_demo$DIAG %>% gsub("VCFS","PATIENT-DEL",.)  %>% gsub("Sibling","CONTROL",.) %>% gsub("Control","CONTROL",.)  %>% as.factor
# RACE as factor
SUNY_demo$RACE[is.na(SUNY_demo$RACE)] <- 7
SUNY_demo$RACE <- factor(SUNY_demo$RACE,levels=c(1:7),labels=c("1_Native_American","2_Asian","3_Pacific_Island","4_Black","5_White","6_Multiple","7_Unknown"))
# rename ETHN to HISPANIC and set to factor (1=y, 2=n)
SUNY_demo$HISPANIC <- SUNY_demo$ETHN
SUNY_demo$HISPANIC[is.na(SUNY_demo$HISPANIC)] <- "Unknown"
SUNY_demo$HISPANIC <- factor(SUNY_demo$HISPANIC,levels=c(1,2,"Unknown"),labels=c("Y","N","Unknown"))

# read scan dates (originally pulled from dicom headers)
SUNY_scan_dates <- read.csv(file.path(demopath,"multisite/SUNY_scan_dates.csv"),header=T, strip.white=T, sep=",")
# merge with demo
SUNY_demo_scan_dates <- merge(x=SUNY_scan_dates, y=SUNY_demo[,c("MRI_S_ID","SUBJECT_IDENTITY","DOB","SEX","RACE","HISPANIC")], all.x=T, by="MRI_S_ID")
# get age as difference between DOB and scan date
SUNY_age <- as.data.frame(as.numeric(difftime(strptime(SUNY_demo_scan_dates$scan_date, format="%m/%d/%y"),strptime(SUNY_demo_scan_dates$DOB, format="%Y-%m-%d"),units="days"))/365)
colnames(SUNY_age) <- "AGE"
df_SUNY_covars <- cbind(SUNY_demo_scan_dates,SUNY_age)

# get med info from T4_VCFS_Meds.xlsx
SUNY_meds <- read_xlsx(file.path(demopath,"multisite/SUNY_T4_VCFS_Meds.xlsx"),trim_ws=T, na="", col_names=T)
SUNY_meds$MRI_S_ID <- SUNY_meds$`SUBJECT#` %>% gsub("S","X",.)
SUNY_med_class <- SUNY_meds[,c("CurMED1Class","CurMED2Class","CurMED3Class","CurMED4Class")]
SUNY_apd_tf <- lapply(1:nrow(SUNY_med_class), function(r) "ANTIPSYCHOTIC" %in% (SUNY_med_class[r,] %>% toString %>% gsub('[[:punct:] ]+',' ',.) %>% strsplit(.," +") %>% unlist %>% toupper)) %>% do.call(rbind,.)
SUNY_meds_apd <- cbind(SUNY_meds, SUNY_apd_tf)
SUNY_meds_apd$Med_Antipsychotic <- factor(SUNY_meds_apd$SUNY_apd_tf, levels=c(T,F), labels=c("Y","N"))
# merge with demo
SUNY_demo_use <- merge(x=df_SUNY_covars, y=SUNY_meds_apd[,c("MRI_S_ID","Med_Antipsychotic")], by="MRI_S_ID")

# get IQ from T4_WAIS.xlsx
SUNY_neuro <- read_xlsx(file.path(demopath,"multisite/SUNY_T4_WAIS.xlsx"),trim_ws=T, na="", col_names=T)
SUNY_neuro$IQ_full <- SUNY_neuro$T4WAISFSIQ
SUNY_neuro$IQ_measure <- "WAIS_full_scale"
SUNY_neuro$MRI_S_ID <- SUNY_neuro$`SUBJECT#` %>% gsub("S","X",.)
SUNY_demo_use <- merge(x=SUNY_demo_use, y=SUNY_neuro[,c("MRI_S_ID","IQ_full","IQ_measure")], by="MRI_S_ID", all.x=T)

# get psych dx from scid
SUNY_scid1 <- read_xlsx(file.path(demopath,"multisite/SUNY_T4_SCID_Page2_6.xlsx"),trim_ws=T, na="", col_names=T) 
SUNY_scid1$MRI_S_ID <- SUNY_scid1$`SUBJECT#` %>% gsub("S","X",.)
# second scid page doesn't have SCZ spectrum diagnoses, ignoring for now
#SUNY_scid2 <- read_xlsx("~/Desktop/22q_multisite/SUNY_csv/SUNY_xlsx/T4_SCID_Page7.xlsx",trim_ws=T, na="", col_names=T)
# relevant variables: Schizophrenia_Curr, SZphreniformDis_Curr, SZaffectiveDis_Curr, DelusionalDis_Curr, BriefPsychDis_Curr, PsychDisNOS_Curr

# "3" indicates positive dx, "1" is no dx
#filter(SUNY_scid1, Schizophrenia_Curr==3 | SZphreniformDis_Curr==3 | SZaffectiveDis_Curr==3 | DelusionalDis_Curr==3 | BriefPsychDis_Curr==3 | PsychDisNOS_Curr==3)[,c("MRI_S_ID","Schizophrenia_Curr", "SZphreniformDis_Curr", "SZaffectiveDis_Curr", "DelusionalDis_Curr", "BriefPsychDis_Curr", "PsychDisNOS_Curr")]
suny_psychdx_subs <- filter(SUNY_scid1, Schizophrenia_Curr==3 | SZphreniformDis_Curr==3 | SZaffectiveDis_Curr==3 | DelusionalDis_Curr==3 | BriefPsychDis_Curr==3 | PsychDisNOS_Curr==3)$MRI_S_ID %>% as.array
SUNY_demo_use$psych_dx <- factor(SUNY_demo_use$MRI_S_ID %in% suny_psychdx_subs, levels=c(T,F), labels=c("Y","N")) 

# replicate MRI_S_ID as SUBJECTID
SUNY_demo_use$SUBJECTID <- SUNY_demo_use$MRI_S_ID

# cardiac
SUNY_dev <- read_xlsx(file.path(demopath,"multisite/SUNY_DEV-Q.xlsx"),trim_ws=T, na="", col_names=T)
SUNY_dev$SUBJECTID <- SUNY_dev$SUBJECT. %>% gsub("S","X",.)
SUNY_dev[which(is.na(SUNY_dev$Heart)),"Heart"] <- 0
SUNY_dev$cardiac <- factor(SUNY_dev$Heart, levels=c(0, 1), labels=c("N","Y"))
# merge cardiac info with dataframe and set NA to FALSE then convert to numeric
SUNY_demo_use <- merge(x=SUNY_demo_use, y=SUNY_dev[,c("SUBJECTID","cardiac")], by=c("SUBJECTID"), all.x=TRUE)
#SUNY_demo_use[which(is.na(SUNY_demo_use$cardiac)),"cardiac"] <- FALSE


# filter based on suny list
ages_suny_subset <- filter(SUNY_demo_use, SUNY_demo_use$MRI_S_ID %in% sessions_suny)[,c("MRI_S_ID","SUBJECTID","SUBJECT_IDENTITY","AGE","SEX","RACE","HISPANIC","Med_Antipsychotic","psych_dx","IQ_full","IQ_measure","cardiac")]
ages_suny_subset$visit_index <- rep(as.factor(1),times=nrow(ages_suny_subset))
ages_suny_subset$Site <- rep("SUNY",times=nrow(ages_suny_subset))
ages_suny_subset

########################################################################################################################
## IoP
# https://www.sistat.ucla.edu/22qiop/Login.asp
# Charlie; Bruin123
# read original demographics file
iop_demo <- read_xlsx(file.path(demopath,"multisite/22Qiop_Demographic.xlsx"),trim_ws=T, na=c("",-9999,-9998), col_names=T)

# read manually edited version that matches scan ID (GQAIMS*) to website ID (k_*)
iop_manual <- read.csv(file.path(demopath,"multisite/IoP_demo_manual_edit.csv"),header=T, strip.white=T, sep=",")
iop_manual$Website.ID <- gsub("k","K",iop_manual$Website.ID)
iop_demo_rename <- merge(x=iop_demo, y=iop_manual[,c("Scan.ID","Website.ID")], by.x="SUBJECTID", by.y="Website.ID")
iop_demo_rename$MRI_S_ID <- iop_demo_rename$Scan.ID
# filter based on lists above
iop_demo_use <- filter(iop_demo_rename, MRI_S_ID %in% sessions_iop)
# SUBJECT_IDENTITY from identity.  
iop_demo_use$SUBJECT_IDENTITY <- iop_demo_use$identity %>% gsub("1","PATIENT-DEL",.) %>% gsub("2","CONTROL",.) %>% as.factor
# change sex coding from 0/1 to F/M and set to factor
iop_demo_use$SEX <- factor(iop_demo_use$sex,levels=c(0,1),labels=c("F","M"))
# set race as factor
iop_demo_use$RACE <- iop_demo_use$race
iop_demo_use$RACE[is.na(iop_demo_use$RACE)] <- 7
iop_demo_use$RACE <- factor(iop_demo_use$RACE,levels=c(1:7),labels=c("1_Native_American","2_Asian","3_Pacific_Island","4_Black","5_White","6_Multiple","7_Unknown"))
# ethnicity as factor with 0=N 1=Y
iop_demo_use$HISPANIC <- iop_demo_use$hispanic
iop_demo_use$HISPANIC[is.na(iop_demo_use$HISPANIC)] <- "Unknown"
iop_demo_use$HISPANIC <- factor(iop_demo_use$HISPANIC,levels=c(0,1,"Unknown"),labels=c("N","Y","Unknown"))
# AGE
iop_demo_use$AGE <- iop_demo_use$age

# get antipsychotic status from 22Qiop_summPsych
iop_summPsych <- read_xlsx(file.path(demopath,"multisite/22Qiop_summPsych.xlsx"),trim_ws=T, na=c("",-9999,-9998), col_names=T)
# extract unique words from summPsych, order by decreasing length for convenience and print
spwords <- iop_summPsych$summPsychNotes %>% toString %>% gsub('[[:punct:] ]+',' ',.) %>% strsplit(.," +") %>% unlist %>% unique 
spwords2 <- cbind(spwords, as.integer(nchar(spwords))) %>% .[order(as.numeric(.[,2])),1] %>% rev %>% as.data.frame
print(spwords2)
# manually observe words and create a list of all strings corresponding to an antipsychotic medication name
iop_apd_strings_manual <- c("AMILSULPRIDE","ARIPRIPAZOLE","OLANZAPINE") 
# test each row of summPsych notes if it contains any of the above strings
iop_apd_tf <- lapply(iop_summPsych$summPsychNotes, function(r) iop_apd_strings_manual %in% unlist(strsplit(r," +")) %>% any) %>% do.call(rbind,.)
iop_summPsych$Med_Antipsychotic <- factor(iop_apd_tf, levels=c(T,F), labels=c("Y","N"))
# psychosis dx
iop_summPsych$psych_dx <- factor(iop_summPsych$psyDiagnos, levels=c(1,0,3), labels=c("Y","N","N"))
# merge with demo
iop_demo_use <- merge(x=iop_demo_use, y=iop_summPsych[,c("SUBJECTID","Med_Antipsychotic","psych_dx")], by="SUBJECTID")

# get IQ scaled score (WASI_full_scale) from 22Qiop_DKEFS.xlsx
iop_neuro <- read_xlsx(file.path(demopath,"multisite/22Qiop_DKEFS.xlsx"),trim_ws=T, na=c("",-9999,-9998), col_names=T)
iop_neuro$IQ_measure <- "WASI_full_scale"
iop_neuro$IQ_full <- iop_neuro$IQSS
iop_demo_use <- merge(x=iop_demo_use, y=iop_neuro[,c("SUBJECTID","IQ_full","IQ_measure")], by="SUBJECTID")

# cardiac diagnoses 
iop_pheno <- read_xlsx(file.path(demopath,"multisite/22Qiop_phenotype.xlsx"),trim_ws=T, na=c("",-9999,-9998), col_names=T)
iop_pheno_card <- iop_pheno[complete.cases(iop_pheno[,c("SUBJECTID", paste0("heart",1:20))]),c("SUBJECTID", paste0("heart",1:20))]
iop_card_tf <- data.frame(SUBJECTID=unique(iop_pheno_card$SUBJECTID) )
iop_card_tf$cardiac <- lapply(iop_card_tf$SUBJECTID,function(s)sum(filter(iop_pheno_card, SUBJECTID==s)[,paste0("heart",1:20)])>0) %>% do.call(rbind,.)

# merge cardiac info with dataframe and set NA to FALSE then convert to numeric
iop_demo_use <- merge(x=iop_demo_use, y=iop_card_tf, by=c("SUBJECTID"), all.x=TRUE)
iop_demo_use[which(is.na(iop_demo_use$cardiac)),"cardiac"] <- FALSE
iop_demo_use$cardiac <- factor(iop_demo_use$cardiac, levels=c(FALSE,TRUE), labels=c("N","Y"))

# take only desired variables 
ages_iop_subset <- iop_demo_use[,c("MRI_S_ID","SUBJECTID","SUBJECT_IDENTITY","AGE","SEX","RACE","HISPANIC","Med_Antipsychotic","psych_dx","IQ_full","IQ_measure","cardiac")]
# uncomment to check if sessions are missing
#iop_not_in_demo <- sessions_iop[which(sessions_iop %in% iop_demo_use$MRI_S_ID == F)]
# Add site and visit index
ages_iop_subset$visit_index <- rep(as.factor(1),times=nrow(ages_iop_subset))
ages_iop_subset$Site <- rep("IoP",times=nrow(ages_iop_subset))
ages_iop_subset

########################################################################################################################
## Rome
#rome_demo <- read.csv("~/Desktop/22q_multisite/DATA_ENIGMA_ROME_manual_edit.csv",header=T, strip.white=T, sep=",")
rome_demo <- read_xlsx(file.path(demopath,"multisite/ENIGMA_DB_ROME.xlsx"),trim_ws=T, na=c("","NA"), col_names=T)
rome_demo_use <- filter(rome_demo, SubjID %in% sessions_rome)
#sessions_rome %in% rome_demo_use$SubjID

# SUBJECT IDENTITY from Dx (1=PATIENT-DEL, 2=CONTROL)
rome_demo_use$SUBJECT_IDENTITY <- rome_demo_use$Dx %>% sub(1,"PATIENT-DEL",.) %>% sub(0,"CONTROL",.) %>% as.factor
# rename Age
rome_demo_use$AGE <- rome_demo_use$Age
# sex as factor (1=M, 2=F)
rome_demo_use$SEX <- factor(rome_demo_use$Sex,levels=c(1,2),labels=c("M","F"))
# race (all white)
rome_demo_use$RACE <- rep(5, times=nrow(rome_demo_use)) 
rome_demo_use$RACE <- factor(rome_demo_use$RACE,levels=c(1:7),labels=c("1_Native_American","2_Asian","3_Pacific_Island","4_Black","5_White","6_Multiple","7_Unknown"))
# ethnicity (all non-hispanic)
rome_demo_use$HISPANIC <- rep(0,times=nrow(rome_demo_use))
rome_demo_use$HISPANIC <- factor(rome_demo_use$HISPANIC,levels=c(0,1,"Unknown"),labels=c("N","Y","Unknown"))
# MRI_S_ID
rome_demo_use$MRI_S_ID <- rome_demo_use$SubjID
# Current Antisychotic (Y/N) -- none of the Rome subjects used are taking a typical antipsychotic, so just use Atypical column
rome_demo_use$Current_Atypical_Antipsychotic[is.na(rome_demo_use$Current_Atypical_Antipsychotic)] <- 0
rome_demo_use$Med_Antipsychotic <-factor(rome_demo_use$Current_Atypical_Antipsychotic, levels=c(0,1), labels=c("N","Y"))
rome_demo_use$IQ_full <- rome_demo_use$IQ
rome_demo_use$IQ_measure <- rome_demo_use$IQ_Method

# psychosis dx
rome_demo_use$Psychotic_Disorder[is.na(rome_demo_use$Psychotic_Disorder)] <- 0
rome_demo_use$psych_dx <- factor(rome_demo_use$Psychotic_Disorder,levels=c(0,1),labels=c("N","Y"))

# replicate MRI_S_ID as SUBJECTID
rome_demo_use$SUBJECTID <- rome_demo_use$MRI_S_ID

# cardiac

rome_heart <- read_xlsx(file.path(demopath,"multisite/CHD_del22_Rome_Enigma.xlsx"),trim_ws=T, na="", col_names=T) %>% rename("SUBJECTID"="SubjID")
rome_heart$cardiac <- factor(rome_heart$CHD, levels=c(1, 2), labels=c("Y","N"))
# merge cardiac info with dataframe and set NA to FALSE then convert to numeric
rome_demo_use <- merge(x=rome_demo_use, y=rome_heart[,c("SUBJECTID","cardiac")], by=c("SUBJECTID"), all.x=TRUE)
rome_demo_use[which(is.na(rome_demo_use$cardiac)),"cardiac"] <- "N"

# take only desired variables 
ages_rome_subset <- rome_demo_use[,c("MRI_S_ID","SUBJECTID","SUBJECT_IDENTITY","AGE","SEX","RACE","HISPANIC","Med_Antipsychotic","psych_dx","IQ_full","IQ_measure","cardiac")]
ages_rome_subset$visit_index <- rep(as.factor(1),times=nrow(ages_rome_subset))
ages_rome_subset$Site <- rep("Rome",times=nrow(ages_rome_subset))
ages_rome_subset %<>% as.data.frame
ages_rome_subset

############################################################################################################
# Combine all
multi_demo <- rbind(ages_trio_subset,ages_prisma_subset,ages_suny_subset,ages_iop_subset,ages_rome_subset)
multi_demo$AGE <- signif(as.numeric(multi_demo$AGE), digits=2)
multi_demo$Med_Antipsychotic[is.na(multi_demo$Med_Antipsychotic)] <- "N"
multi_demo$psych_dx[is.na(multi_demo$psych_dx)] <- "N"
multi_demo <- filter(multi_demo, SUBJECT_IDENTITY == "PATIENT-DEL" | SUBJECT_IDENTITY == "CONTROL")
multi_demo[which(is.na(multi_demo$cardiac)),"cardiac"] <- "N"

#View(multi_demo)
#write.csv(multi_demo, file=file.path(demopath,"multisite/22q_multisite_demo.csv"),row.names = FALSE)


