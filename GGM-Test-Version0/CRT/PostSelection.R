

MethodSelection = function(setting) {
  return(switch(
    setting,
    'l_l_w' = list(
      ind = c(1, 2, 3, 4),
      names = c('LM-SST', 'LM-SSR', 'F-test', 'dCRT'),
      theta_lim = c(0, 2.5)
    ),
    'h_l_w' =  list(
      ind = c(1, 2,3, 5),
      names = c('LM-L1-R-SST', 'LM-L1-R-SSR', 'De-Lasso',  'dCRT'),
      theta_lim = c(0, 1.2)
    ),  
    
    'l_gl_w' = list(
      ind = c(1, 2, 3, 4),
      names = c('GLM-Dev', 'RF', 'Chi-square', 'dCRT'),
      theta_lim = c(0, 8)
    ),
    'h_gl_w' =  list(
      ind = c(1,2,3,4,5),
      names = c('GLM-L1-D', 'GLM-L1-R-SST', 'RF', 'CGM', 'dCRT'),
      theta_lim = c(0, 3)
    ),
    
    
    
    'l_l_m' =  list(
      ind = c(1,2,3 ,4, 6)  ,
      names = c('RF', 'RF-D' ,'RF-RR',  'F-test', 'dCRT (RF)'),
      theta_lim = c(0, 7)
    ),
    'h_l_m' =  list(
      ind = c(1,2,3, 4,7  )    ,  
      names = c('RF', 'RF-D' ,'RF-RR', 'De-Lasso', 'dCRT (RF)'),
      theta_lim = c(0, 7) 
    ), 
    
    
    
    
    
    'l_gl_m' =  list(
      ind = c(1,2,3,4,6 ), 
      names = c('RF','RF-D' ,'RF-RR',  'Chi-square', 'dCRT (RF)'),
      theta_lim = c(0, 7)
    ),
    'h_gl_m' =  list(
      ind = c(1,2,3, 4,6),  
      names = c('RF', 'RF-D' ,'RF-RR', 'CGM', 'dCRT (RF)'),
      theta_lim = c(0, 7)
    ),
  ))
}
