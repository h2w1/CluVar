library(pheatmap)
library(gridExtra)
library(readr)
library(Matrix)
library(DENDRO)



CluVar_simulation <- function(kprob = NULL, lprob = NULL,
                                          epi = 0, sampcells,
                                          sampgenes, count, k = NULL, subtype = NULL, verbose = FALSE) {

  
  
  alt_select <- count$alt[sampgenes, sampcells] # sampling에서 선택된 alt count 
  ref_select <- count$ref[sampgenes, sampcells] # sampling에서 선택된 ref count
  
  

  
  # 클레이드 수(k) 설정
  if (!is.null(kprob)) {
    k <- length(kprob)  # kprob이 주어지면 그 길이가 k
  } else if (is.null(k)) {
    k <- 3  # 기본값으로 k = 3
  }
  if (verbose) {
    cat('There are total ', k, ' clades in our simulation \n')
  }
  
  # 세포 수(m)와 유전자 수(n) 설정
  m <- m  # 세포의 수
  n <- n  # 유전자의 수
  l <- 2 * k - 1
  
  # 세포를 클레이드에 할당
  if (is.null(kprob)) {
    # kprob이 없으면 균등한 확률로 세포를 클레이드에 할당
    clades <- sample(seq(1, l), length(sampcells), replace = TRUE)
  } else {
    if (length(kprob) < 2 * k - 1) {
      stop('Error with kprob, kprob should have length  2 * k - 1')
    } else {
      if (!isTRUE(all.equal(sum(lprob), 1))) {
        stop('Error with iprob, lprob sum should be 1')
      }
      # kprob이 있으면 해당 확률에 따라 세포를 클레이드에 할당
      clades <- sample(seq(1, l), length(sampcells), prob = kprob, replace = TRUE)
    }
  }
  
  # 유전자(돌연변이 위치)를 분기에 할당
  if (is.null(lprob)) {
    # 계통수에서 가능한 분기의 수 (루트 포함)
    cladegene <- sample(seq(1, l), length(sampgenes), replace = TRUE)  # 유전자를 분기에 랜덤하게 할당
  } else {
    if (length(lprob) < 2 * k - 1) {
      stop('Error with lprob, lprob should have length 2 * k - 1')
    } else {
      # lprob or kprob sum should be 1 
      if (!isTRUE(all.equal(sum(lprob), 1))) {
        stop('Error with iprob, lprob sum should be 1')
      }
      cladegene <- sample(seq(1, l), length(sampgenes), prob = lprob, replace = TRUE)
    }
  }
  
  ALT <- TOTAL <- FRAC <- BI <- BI.true <- matrix(NA, ncol = m, nrow = n)
  
  
  
   if (k == 2) {
    for (i in 1:n) {  # 유전자(돌연변이 위치)마다 반복
      alt <- alt_select[i, ]  # 돌연변이 대립유전자
      ref <- ref_select[i, ]
      TOTAL[i, ] <- alt + ref 
      if (cladegene[i] == 1) { # cladegne 1에 해당하는 유전자들은 
        for (j in 1:m) { # 세포를 돌면서 
          if (!is.na(alt[j])) { # alt도 0이 아니고, ref도 0이 아니면   
            if (clades[j] == 1) { # clade가 1인 세포에 대해서 
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              } else {
                ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                BI[i, j] <- 1
              }
              } else {
              ALT[i, j] <- rbinom(1, TOTAL[i, j], epi)
              BI[i, j] <- 0
            }
          } else {  # alt[j]가 NA인 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
          }
        
          if (clades[j] == 1) {
              BI.true[i, j] <- 1
          }else{
              BI.true[i, j] <- 0 
            }
          }
        }else if (cladegene[i] == 2) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]가 NA가 아닐 경우
              if (clades[j] == 2) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                    }else{
                      ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                } else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                }
                } else {
                ALT[i, j] <- rbinom(1, TOTAL[i, j], epi)
                BI[i, j] <- 0
                }
              }  
             else {  # alt[j]가 NA인 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            
          if (clades[j] == 2) {
                BI.true[i, j] <- 1
            }else{
                BI.true[i,j] <- 0 
             }
           }
          }else if (cladegene[i] == 3) {
        for (j in 1:m) {
          if (!is.na(alt[j])&!is.na(ref[j])) {  # alt[j]가 NA가 아닐 경우
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              } else {
                ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                BI[i, j] <- 1
              }

          } else {  # alt[j]가 NA인 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
          }
          
          BI.true[i, j] <- 1
          }
        }
      }
   }      
   
  if (k == 3) {
    for (i in 1:n) {
      alt <- alt_select[i, ]  # 돌연변이 대립유전자
      ref <- ref_select[i, ]
      TOTAL[i, ] <- alt + ref
      if (cladegene[i] == 1) {
        for (j in 1:m) {
          if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
            if (clades[j] == 1) {
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              }
              else {
              ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
              BI[i, j] <- 1
              } 
            }
              else {
              ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
              BI[i, j] <- 0
            }
        } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
        }
          if (clades[j] == 1) {
            BI.true[i, j] <- 1
          }else{
            BI.true[i,j] <-0
          }
        }
    } else if (cladegene[i] == 2) {
        for (j in 1:m) {
          if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
            if (clades[j] == 2) {
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              }
              else {
                ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                BI[i, j] <- 1
              } 
            } else {
            ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
            BI[i, j] <- 0
            }
        } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
        }
          if (clades[j] == 2) {
            BI.true[i, j] <- 1
          }else{
            BI.true[i,j] <-0
          }
        }  
    } else if (cladegene[i] == 3) {
        for (j in 1:m) {
          if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
            if (clades[j] == 3) {
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              }
              else {
                ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                BI[i, j] <- 1
              } 
            
          } else {
             ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
              BI[i, j] <- 0
            }
        } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
        }
          if (clades[j] == 3) {
            BI.true[i, j] <- 1
          }else{
            BI.true[i,j] <-0
          }
        }
    } else if (cladegene[i] == 4) {
        for (j in 1:m) {
          if (!is.na(alt[j]) & !is.na(ref[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
            if (clades[j] %in% c(1, 2, 4)) {
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              }
              else {
                ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                BI[i, j] <- 1
              } 
              
            } else {
              ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
              BI[i, j] <- 0
            }
        } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
        }
          if (clades[j] %in% c(1, 2, 4)) {
            BI.true[i, j] <- 1
          }else{
            BI.true[i,j] <-0
          }
        }
          
      } else if (cladegene[i] == 5) {
        for (j in 1:m) {
          if (!is.na(alt[j]) & !is.na(ref[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
            if (alt[j] == 0) {
              if (runif(1) <= 1 - epi) {
                ALT[i, j] <- ref[j]
              }else{
                ALT[i, j] <- 0
              }
              BI[i, j] <- 1
            }
            else {
              ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
              BI[i, j] <- 1
            } 
            
          } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
          }
          BI.true[i, j] <- 1
          }
        }
      }
    }  
    
  
  if (k == 4) {
    if (subtype == 1) {
      for (i in 1:n) {
        alt <- alt_select[i, ]  # 돌연변이 대립유전자
        ref <- ref_select[i, ]
        TOTAL[i, ] <- alt + ref
        if (cladegene[i] == 1) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] == 1) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] == 1) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <-0
            }
          }
      } else if (cladegene[i] == 2) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] == 2) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] == 2) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <-0
            }
          }
     } else if (cladegene[i] == 3) {
        for (j in 1:m) {
          if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
            if (clades[j] == 3) {
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              }
              else {
                ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                BI[i, j] <- 1
              } 
              
            } else {
              ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
              BI[i, j] <- 0
            }
          } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
          }
          if (clades[j] == 3) {
            BI.true[i, j] <- 1
          }else{
            BI.true[i,j] <- 0
          } 
        }
    } else if (cladegene[i] == 4) {
        for (j in 1:m) {
          if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
            if (clades[j] == 4) {
              if (alt[j] == 0) {
                if (runif(1) <= 1 - epi) {
                  ALT[i, j] <- ref[j]
                }else{
                  ALT[i, j] <- 0
                }
                BI[i, j] <- 1
              }
              else {
                ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                BI[i, j] <- 1
              } 
              
            } else {
              ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
              BI[i, j] <- 0
            }
          } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
            ALT[i, j] <- NA
            BI[i, j] <- NA
          }
          if (clades[j] == 4) {
            BI.true[i, j] <- 1
          }else{
            BI.true[i,j] <- 0
          }
        }
   } else if (cladegene[i] == 5) {
       for (j in 1:m) {
         if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
           if (clades[j] %in% c(1, 2, 5)) {
             if (alt[j] == 0) {
               if (runif(1) <= 1 - epi) {
                 ALT[i, j] <- ref[j]
               }else{
                 ALT[i, j] <- 0
               }
               BI[i, j] <- 1
             }
             else {
               ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
               BI[i, j] <- 1
             } 
             
           } else {
             ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
             BI[i, j] <- 0
           }
        } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
          ALT[i, j] <- NA
          BI[i, j] <- NA
        }
         if (clades[j] %in% c(1, 2, 5)) {
           BI.true[i, j] <- 1
         }else{
           BI.true[i,j] <- 0
         }
      }
    } else if (cladegene[i] == 6) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] %in% c(3, 4, 6)) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] %in% c(3, 4, 6)) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <- 0
            }
          }
        } else if (cladegene[i] == 7) {
            for (j in 1:m) {
              if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
                ALT[i, j] <- NA
                BI[i, j] <- NA
              }
              BI.true[i, j] <- 1
            }
          
          }  
        }
      }else {
      for (i in 1:n) {
        alt <- alt_select[i, ]  # 돌연변이 대립유전자
        ref <- ref_select[i, ]
        TOTAL[i, ] <- alt + ref
        if (cladegene[i] == 1) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] == 1) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] == 1) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <- 0
            }
          }
      } else if (cladegene[i] == 2) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] == 2) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] == 2) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <- 0
            }
          }
      } else if (cladegene[i] == 3) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] == 3) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] == 3) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <- 0
            }
          }
      } else if (cladegene[i] == 4) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] == 4) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                  ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)  
                  BI[i, j] <- 0
                }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
                ALT[i, j] <- NA
                BI[i, j] <- NA
            }
            if (clades[j] == 4) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <- 0
            }
          }
       } else if (cladegene[i] == 5) {
          for (j in 1:m) {
            if (!is.na(alt[j]) & !is.na(ref[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] %in% c(1, 2, 5)) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] %in% c(1, 2, 5)) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <- 0
            }
          }
       } else if (cladegene[i] == 6) {
          for (j in 1:m) {
            if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
              if (clades[j] %in% c(1, 2, 3, 5, 6)) {
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {
                ALT[i, j] <- rbinom(1, TOTAL[i,j], epi)
                BI[i, j] <- 0
              }
            } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
              ALT[i, j] <- NA
              BI[i, j] <- NA
            }
            if (clades[j] %in% c(1, 2, 3, 5, 6)) {
              BI.true[i, j] <- 1
            }else{
              BI.true[i,j] <- 0
            }
          }
        } else if (cladegene[i] == 7) {
            for (j in 1:m) {
              if (!is.na(alt[j])) {  # alt[j]와 ref[j] 모두 NA가 아닐 경우
                if (alt[j] == 0) {
                  if (runif(1) <= 1 - epi) {
                    ALT[i, j] <- ref[j]
                  }else{
                    ALT[i, j] <- 0
                  }
                  BI[i, j] <- 1
                }
                else {
                  ALT[i, j] <- rbinom(1, alt[j], 1 - epi) + rbinom(1, ref[j], epi)
                  BI[i, j] <- 1
                } 
                
              } else {  # alt[j] 또는 ref[j]에 NA가 있을 경우
                ALT[i, j] <- NA
                BI[i, j] <- NA
              }
              BI.true[i,j] <- 1
            }  
          }  
        }
      }
    }
  
  # ALT와 REF에서 같은 위치에 0이 있을 때 해당 값을 NA로 변경

  FRAC <- ifelse(is.na(ALT) & is.na(TOTAL), NA, ALT / TOTAL)
  

        
  
  return(list(ALT = ALT, TOTAL = TOTAL, FRAC = FRAC, BI = BI , BI.true = BI.true, clades = clades))
}  


calculate_ratios <- function(ALT) {
  # NA 값의 비율을 계산
  na_ratio <- sum(is.na(ALT)) / length(ALT)
  
  # 0 값의 비율을 계산
  zero_ratio <- sum(ALT == 0, na.rm = TRUE) / length(ALT)
  
  # 0이 아닌 값의 비율을 계산
  non_zero_ratio <- sum(ALT != 0, na.rm = TRUE) / length(ALT)
  
  # 결과 출력
  cat("0 비율:", round(zero_ratio, 3), "\n")
  cat("변이(1)의 비율:", round(non_zero_ratio, 3), "\n")
  cat("missing 비율:", round(na_ratio, 3), "\n")
}


# 함수 호출 예시
# calculate_ratios(ALT)

# 함수 정의
select_rows_by_na <- function(sim, num_rows) {
  # 각 행에서 NA 값의 개수를 계산
  na_counts <- apply(sim$ALT , 1, function(x) sum(is.na(x)))
  
  # NA 값이 적은 순서대로 행의 인덱스를 정렬
  sorted_indices <- order(na_counts)
  
  # NA 값이 적은 순으로 상위 num_rows개의 행 선택
  selected_rows <- sorted_indices[1:num_rows]
  # 선택된 행들로 새로운 행렬 만들기
  sim$ALT <- sim$ALT[selected_rows, ]
  sim$TOTAL <- sim$TOTAL[selected_rows, ]
  sim$BI <- sim$BI[selected_rows, ]
  sim$BI.true <- sim$BI.true[selected_rows, ]
  sim$FRAC <- sim$FRAC[selected_rows, ]
  sim$clades <-sim$clades[selected_rows]
  
  
  # 리스트로 반환 (여러 값 반환)
  return(sim)
}

structure_heatmap <- function(BI_matrix) {
  
  # 행에서 1의 비율 계산 후, 1의 비율이 높은 순으로 정렬
  row_ones_ratio <- rowSums(BI_matrix == 1) / ncol(BI_matrix)
  BI_sorted_rows <- BI_matrix[order(-row_ones_ratio), ]
  
  # 열에서 1의 비율 계산 후, 1의 비율이 높은 순으로 정렬
  col_ones_ratio <- colSums(BI_sorted_rows == 1) / nrow(BI_sorted_rows)
  BI_sorted_final <- BI_sorted_rows[, order(-col_ones_ratio)]
  
  # 히트맵 생성
  p <- pheatmap(BI_sorted_final, cluster_rows = FALSE, cluster_cols = FALSE, 
                color = colorRampPalette(c("navy", "yellow"))(50),
                border_color = NA)
  
  # 결과 출력
  grid::grid.draw(p$gtable) # 경계선을 제거하고 다시 그리기
  
  # 히트맵 객체 반환 (필요한 경우)
  return(p)
}


# 0 값의 비율이 높은 순으로 행과 열을 모두 정렬하는 히트맵 함수
count_heatmap <- function(matrix_data, zero_color = "white", non_zero_palette = c("yellow", "red"), na_color = "gray") {
  
  # 각 행에서 0의 비율 계산 (NA 값 제외)
  zero_ratio_row <- rowSums(matrix_data == 0 & !is.na(matrix_data)) / rowSums(!is.na(matrix_data))
  
  # 0의 비율이 높은 순서대로 행 정렬 (낮은 비율이 먼저 오도록 오름차순 정렬)
  sorted_matrix <- matrix_data[order(zero_ratio_row), ]
  
  # 각 열에서 0의 비율 계산 (NA 값 제외)
  zero_ratio_col <- colSums(sorted_matrix == 0 & !is.na(sorted_matrix)) / colSums(!is.na(sorted_matrix))
  
  # 0의 비율이 높은 순서대로 열 정렬 (낮은 비율이 먼저 오도록 오름차순 정렬)
  sorted_matrix <- sorted_matrix[, order(zero_ratio_col)]
  
  # 데이터 값의 최소값과 최대값을 구함
  min_value <- min(sorted_matrix, na.rm = TRUE)
  max_value <- max(sorted_matrix, na.rm = TRUE)
  
  # breaks 설정: 0 근처 구간을 더 세밀하게 나누고, min_value가 0보다 큰지 확인
  if (min_value < 0) {
    breaks <- c(seq(min_value, 0, length.out = 25), seq(0, max_value, length.out = 26))  # 0을 기준으로 세밀하게 나눔
  } else {
    breaks <- seq(min_value, max_value, length.out = 51)  # min_value가 0 이상이면 일반적인 구간 설정
  }
  
  # 사용자 정의 색상 팔레트 (0은 zero_color, 그 외는 그라데이션)
  custom_colors <- c(zero_color, colorRampPalette(non_zero_palette)(50))
  
  # 히트맵 생성
  p <- pheatmap(sorted_matrix, cluster_rows = FALSE, cluster_cols = FALSE, 
                color = custom_colors,
                breaks = breaks,
                na_col = na_color,  # NA 값을 na_color로 표시
                border_color = NA)
  
  # 히트맵을 그리기
  grid::grid.draw(p$gtable)
  
  # 히트맵 객체 반환 (필요한 경우)
  return(p)
}

# 0 값의 비율이 높은 순으로 행과 열을 모두 정렬하는 히트맵 함수
count_heatmap <- function(matrix_data, zero_color = "white", non_zero_palette = c("yellow", "red"), na_color = "gray") {
  
  # 각 행에서 0의 비율 계산 (NA 값 제외)
  zero_ratio_row <- rowSums(matrix_data == 0 & !is.na(matrix_data)) / rowSums(!is.na(matrix_data))
  
  # 0의 비율이 높은 순서대로 행 정렬 (낮은 비율이 먼저 오도록 오름차순 정렬)
  sorted_matrix <- matrix_data[order(zero_ratio_row), ]
  
  # 각 열에서 0의 비율 계산 (NA 값 제외)
  zero_ratio_col <- colSums(sorted_matrix == 0 & !is.na(sorted_matrix)) / colSums(!is.na(sorted_matrix))
  
  # 0의 비율이 높은 순서대로 열 정렬 (낮은 비율이 먼저 오도록 오름차순 정렬)
  sorted_matrix <- sorted_matrix[, order(zero_ratio_col)]
  
  # 데이터 값의 최소값과 최대값을 구함
  min_value <- min(sorted_matrix, na.rm = TRUE)
  max_value <- max(sorted_matrix, na.rm = TRUE)
  
  # breaks 설정: 0 근처 구간을 더 세밀하게 나누고, min_value가 0보다 큰지 확인
  if (min_value < 0) {
    breaks <- c(seq(min_value, 0, length.out = 25), seq(0, max_value, length.out = 26))  # 0을 기준으로 세밀하게 나눔
  } else {
    breaks <- seq(min_value, max_value, length.out = 51)  # min_value가 0 이상이면 일반적인 구간 설정
  }
  
  # 사용자 정의 색상 팔레트 (0은 zero_color, 그 외는 그라데이션)
  custom_colors <- c(zero_color, colorRampPalette(non_zero_palette)(50))
  
  # 히트맵 생성
  p <- pheatmap(sorted_matrix, cluster_rows = FALSE, cluster_cols = FALSE, 
                color = custom_colors,
                breaks = breaks,
                na_col = na_color,  # NA 값을 na_color로 표시
                border_color = NA)
  
  # 히트맵을 그리기
  grid::grid.draw(p$gtable)
  
  # 히트맵 객체 반환 (필요한 경우)
  return(p)
}

select_rows_by_na_real <- function(count, num_rows) {
  # 각 행에서 NA 값의 개수를 계산
  na_counts <- apply(count$alt , 1, function(x) sum(is.na(x)))
  
  # NA 값이 적은 순서대로 행의 인덱스를 정렬
  sorted_indices <- order(na_counts)
  
  # NA 값이 적은 순으로 상위 num_rows개의 행 선택
  selected_rows <- sorted_indices[1:num_rows]
  
  # 선택된 행들로 새로운 행렬 만들기
  count$alt <- count$alt[selected_rows, ]
  count$ref <- count$ref[selected_rows, ]
  
  
  
  # 리스트로 반환 (여러 값 반환)
  return(count)
}






######################## 1. 데이터 불러오기 ############################



# 두 개의 TSV 파일 경로
file1 <- "./merged_firtst_filter_alt.tsv"
file2 <- "./merged_firtst_filter_ref.tsv"

# TSV 파일을 헤더 없이 불러오고, integer matrix로 변환
data1 <- as.matrix(read.table(file1, sep = "\t", header = FALSE, row.names = NULL))
data2 <- as.matrix(read.table(file2, sep = "\t", header = FALSE, row.names = NULL))

# 데이터가 numeric일 경우, 정수로 변환
# NA가 있어도 경고 없이 진행되도록 처리
alt <- apply(data1, 2, function(x) as.integer(ifelse(is.na(x), NA, as.numeric(x))))
ref <- apply(data2, 2, function(x) as.integer(ifelse(is.na(x), NA, as.numeric(x))))

# 두 개의 integer matrix를 count라는 리스트에 저장
count <- list(alt = alt,ref = ref)


######################## 2. 파라미터 설정, simulation 실행 ############################


kprob <- NULL # - kprob: 각 클레이드에 세포를 할당할 확률 벡터 (클레이드 수에 따라 길이 결정)
#lprob=NULL # - lprob: 각 분기(edge)에 돌연변이를 할당할 확률 벡터
#lprob <-c(0.11, 0.11, 1 - epi8) # k=2
#lprob <-c(0.06, 0.06, 0.06,0.08,1 - epi4) # k=3
lprob <-c(0.03,0.03,0.03,0.03, 0.05,0.09,1 - epi4) #k=4, subtype=1
#lprob <-c(0.03,0.03,0.03,0.03, 0.06,0.08,1 - epi4) #k=4, subtype=2
m=10000 #세포
n=1335 # 우선 다 사용하고 후에 필터링 하는 방식 
epi=0 # - epi: 시퀀싱 오류율,noise (기본값 0.001) # FPR, FNR 
k=4 # - k:= 2,3,4 클레이드의 수 (계통수에서 분기된 그룹의 수)
subtype=1 # - subtype: k가 4일 때 계통수의 하위 유형을 지정하는 파라미터
sampgenes=sample(seq(1,nrow(count$alt)),n,replace = FALSE) # - sampcells: 샘플링된 세포의 인덱스
sampcells=sample(seq(1,ncol(count$alt)),m,replace = FALSE) # - sampgenes: 샘플링된 유전자의 인덱스 (돌연변이 위치)   



sim=  CluVar_simulation(kprob,lprob,epi,
                                  sampcells,sampgenes,
                                  count,k,subtype,verbose=TRUE)


# 시뮬레이션 결과 
# sim$TOTAL : simulation 결과 변이 + 참조 count
# sim$ALT : simulation 결과 변이 count
# sim$FRAC : simulation 결과  변이 count /총(변이 + 참조) count
# sim$BI : simulation 결과 변이 1, 참조 0
# sim#BI.true: simulation 결과 변이 1,참조 0, NA 값 없는 버전 
# sim$clades : 세포들이 어디 그룹에 속해 있는가, 정답 

####################### simulation 결과 확인 ############################

#####fitering 전
# missing, 변이, 0비율 파악하기
#calculate_ratios(sim$ALT)
#calculate_ratios(sim$TOTAL)
#calculate_ratios(sim$FRAC)
#calculate_ratios(sim$BI)
#calculate_ratios(sim$BI.true)


#####fitering 후 
#원하는 기준 300개 뽑기, missing 낮은순   

filtered_sim = select_rows_by_na(sim, 300)


calculate_ratios(filtered_sim$ALT)
calculate_ratios(filtered_sim$TOTAL)
calculate_ratios(filtered_sim$FRAC)
calculate_ratios(filtered_sim$BI)
calculate_ratios(filtered_sim$BI.true)
#FRAC 에서 0,1이 아닌 값의 비율 확인 
non_zero_one_ratio <- sum(filtered_sim$FRAC != 0 & filtered_sim$FRAC != 1, na.rm = TRUE) / length(filtered_sim$FRAC)
print(non_zero_one_ratio)
#ALT 에서 0,1이 아닌 값의 비율 확인 
non_zero_one_ratio <- sum(filtered_sim$ALT != 0 & filtered_sim$ALT != 1, na.rm = TRUE) / length(filtered_sim$ALT)
print(non_zero_one_ratio)

#히스토그램으로 확인 

# 매트릭스 값을 벡터로 변환 후 히스토그램 생성
hist(as.vector(filtered_sim$ALT), main="ALT Distribution", xlab="Values", ylab="Frequency", col="lightblue")
hist(as.vector(filtered_sim$FRAC), main="FRAC Distribution", xlab="Values", ylab="Frequency", col="lightblue")


########################결과 heatmap############################


# 구조 히트맵 그리기 
structure_heatmap(filtered_sim$BI.true)

BI_ALT <- filtered_sim$ALT
BI_ALT[BI_ALT >= 1 ] <- 1 

# 카운트 구조 히트맵 그리기 
count_heatmap(filtered_sim$BI)
count_heatmap(BI_ALT) 


#######################파일 다운로드 ################################


# 파일을 TSV 형식으로 저장
write.table(filtered_sim$ALT, file = "sim_ALT.tsv", sep = "\t", row.names = FALSE, col.names =FALSE , quote = FALSE)
write.table(filtered_sim$FRAC, file = "sim_FRAC.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(filtered_sim$clades, file = "sim_clades.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




# 1 값의 비율 계산 (NA 제외)
one_ratio <- sum(filtered_sim$ALT == 1, na.rm = TRUE) / sum(!is.na(filtered_sim$ALT))

# 결과 확인
one_ratio




