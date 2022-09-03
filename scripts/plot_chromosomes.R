require(tidyverse)

posFile = '/home/croux/Programmes/pipeline_guillaume/example/sweep/test_positions.txt'
statFile = '/home/croux/Programmes/pipeline_guillaume/example/sweep/test_sumStats.txt'
L = 100000 # length of the chromosome in nucleotides

##################################################################################################################
##################################################################################################################
##################################################################################################################

stat = tibble(read.table(statFile, h=T))
pos = as.numeric(read.table(posFile, h=T))

pi_avg = stat %>% dplyr::select(starts_with("pi_avg_")) %>% summarise(across(.cols=everything(), mean))
thetaW_avg = stat %>% dplyr::select(starts_with("thetaW_avg_")) %>% summarise(across(.cols=everything(), mean))
pi_std = stat %>% dplyr::select(starts_with("pi_std_")) %>% summarise(across(.cols=everything(), mean))
tajD = stat %>% dplyr::select(starts_with("tajD_")) %>% summarise(across(.cols=everything(), mean))
achazY = stat %>% dplyr::select(starts_with("achazY_")) %>% summarise(across(.cols=everything(), mean))
pearson_r = stat %>% dplyr::select(starts_with("pearson_r_")) %>% summarise(across(.cols=everything(), mean))
pearson_pval = stat %>% dplyr::select(starts_with("pearson_pval_")) %>% summarise(across(.cols=everything(), mean))
nHaplo = stat %>% dplyr::select(starts_with("nHaplo_")) %>% summarise(across(.cols=everything(), mean))
H1 = stat %>% dplyr::select(starts_with("H1_")) %>% summarise(across(.cols=everything(), mean))
H2 = stat %>% dplyr::select(starts_with("H2_")) %>% summarise(across(.cols=everything(), mean))
H12 = stat %>% dplyr::select(starts_with("H12_")) %>% summarise(across(.cols=everything(), mean))
H2overH1 = stat %>% dplyr::select(starts_with("H2overH1")) %>% summarise(across(.cols=everything(), mean))
D = stat %>% dplyr::select(starts_with("D_")) %>% summarise(across(.cols=everything(), mean))
r2 = stat %>% dplyr::select(starts_with("r2_")) %>% summarise(across(.cols=everything(), mean))

res = tibble(positions = L*pos, pi_avg=as.numeric(pi_avg), pi_std=as.numeric(pi_std), thetaW=as.numeric(thetaW_avg), tajimaD=as.numeric(tajD), achazY=as.numeric(achazY), pearson_r=as.numeric(pearson_r), pearson_pval=as.numeric(pearson_pval), nHaplo=as.numeric(nHaplo), H1=as.numeric(H1),  H2=as.numeric(H2), H12=as.numeric(H12), H2overH1=as.numeric(H2overH1), D=as.numeric(D), r2=as.numeric(r2))

P_pi = res %>% ggplot(aes(x=positions, y=pi_avg)) + geom_point() + theme_calc(base_size = 15) + ylab('pi (avg)')
P_thetaW = res %>% ggplot(aes(x=positions, y=thetaW)) + geom_point() + theme_calc(base_size = 15) + ylab("Watterson's theta")
P_pi_std = res %>% ggplot(aes(x=positions, y=pi_std)) + geom_point() + theme_calc(base_size = 15) + ylab('pi (std)')
P_Taj = res %>% ggplot(aes(x=positions, y=tajimaD)) + geom_point() + theme_calc(base_size = 15) + ylab("Tajima's D")
P_Y = res %>% ggplot(aes(x=positions, y=achazY)) + geom_point() + theme_calc(base_size = 15) + ylab("Achaz's Y")
P_D = res %>% ggplot(aes(x=positions, y=D)) + geom_point() + theme_calc(base_size = 15) + ylab("LD (D)")
P_r2 = res %>% ggplot(aes(x=positions, y=r2)) + geom_point() + theme_calc(base_size = 15) + ylab("LD (r2)")
P_pearson_r = res %>% ggplot(aes(x=positions, y=pearson_r)) + geom_point() + theme_calc(base_size = 15) + ylab('pearson r')
P_pearson_pval = res %>% ggplot(aes(x=positions, y=pearson_pval)) + geom_point() + theme_calc(base_size = 15) + ylab('pearson pval')
P_nHaplo = res %>% ggplot(aes(x=positions, y=nHaplo)) + geom_point() + theme_calc(base_size = 15) + ylab('nHaplo')
P_H1 = res %>% ggplot(aes(x=positions, y=H1)) + geom_point() + theme_calc(base_size = 15) + ylab('H1')
P_H2 = res %>% ggplot(aes(x=positions, y=H2)) + geom_point() + theme_calc(base_size = 15) + ylab('H2')
P_H12 = res %>% ggplot(aes(x=positions, y=H12)) + geom_point() + theme_calc(base_size = 15) + ylab('H12')
P_H2H1 = res %>% ggplot(aes(x=positions, y=H2overH1)) + geom_point() + theme_calc(base_size = 15) + ylab('H2/H1')

ggarrange(P_pi, P_pi_std, P_thetaW, P_Taj, P_Y, P_nHaplo, P_D, P_r2, P_pearson_r, P_pearson_pval, P_H1, P_H2, P_H12, P_H2H1, align='hv', commond.legend=T, ncol=3, nrow=5)

