#
# The test case from Lange, chapter 5
#
id <- 1:6
momid <- c(0,0,2,2,4,4)
dadid <- c(0,0,1,1,3,3)

xx <- kinship(id, momid, dadid)
aeq(xx, c(4,0,2,2,2,2, 0,4,2,2,2,2, 2,2,4,2,3,3, 
          2,2,2,4,3,3, 2,2,3,3,5,3, 2,2,3,3,3,5) /8)


# And here is an an odd one with cross marriages, but no inbreeding
#
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
round(8*kinship(test1$id, test1$mom, test1$dad))
