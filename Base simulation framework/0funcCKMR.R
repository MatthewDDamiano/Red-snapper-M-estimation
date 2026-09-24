#####################################################################################
# Create dummy data for CKMR sub-model
# dimensions
ckmr_fsamp <- 1
ckmr_lsamp <- 50 # number of samples for each cohort

# vector of cohorts (years of birth)
c1 <- c(rep(60,5), rep(62,5), rep(63,5), rep(61,2), rep(64,2), rep(59,1), rep(65,5), rep(66,2), rep(55,1), rep(67,2),
        rep(68, 10), rep(63,2), rep(61,2), rep(58,1), rep(59,5))
c2 <- c(rep(63,5), rep(65,5), rep(65,5), rep(66,2), rep(67,2), rep(61,1), rep(67,5), rep(68,2), rep(65,1), rep(69,2),
        rep(70, 10), rep(70,2), rep(70,2), rep(70,1), rep(65,5))
# vector of sampling years
s <- c(rep(68,15), rep(69,15), rep(70,20))
# vectors of sibships
# Half-sibling pairs
hsp <- c(rep(0,10),rep(4,1),rep(0,9),rep(12,1),rep(6,1),rep(2,1),rep(0,7),rep(15,1),rep(3,1),rep(0,8),rep(2,1),rep(3,1),
         rep(16,1),rep(0,7))
# Parent-offspring pairs
pop <- c(rep(0,15),rep(1,1),rep(0,4),rep(2,1),rep(0,4),rep(1,1),rep(0,8),rep(2,1),rep(0,10),rep(2,1),rep(0,3),rep(1,1))
# neither (not) hsps or pops
# should be on order of 100s to 1000s of fish
# not <- round(runif(n = 50, min = 10, max = 300))
not<- c(256, 154, 122, 81,  42, 123, 176, 73, 139,  73, 156, 113, 198, 119, 113, 165, 225,  74, 130,  87, 193,
        63, 260, 227, 204, 189, 118, 164, 264, 179, 254, 101, 215,  87, 182, 150, 87, 174, 275, 272,  90, 103,
        296, 190, 282, 145, 128, 201, 54, 176)

ckmr.dat <- list(ckmr_fsamp = ckmr_fsamp, ckmr_lsamp = ckmr_lsamp, c1 = c1, c2 = c2, s = s, hsp = hsp,
                 pop = pop, not = not)
