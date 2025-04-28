
effect <- function(thetas) {
    thetas[["IDE(,gmjo(0))"]] <- thetas[["Y(1,gmjo(0))"]] - thetas[["Y(0,gmjo(0))"]]
    thetas[["IIE_Mjo(1,)"]] <- thetas[["Y(1,gmjo(1))"]] - thetas[["Y(1,gmjo(0))"]]
    thetas[["IIE_M1(1,,gm2(0))"]] <- thetas[["Y(1,gm1(1),gm2(0))"]] - thetas[["Y(1,gm1(0),gm2(0))"]]
    thetas[["IIE_M2(1,gm1(1),)"]] <- thetas[["Y(1,gm1(1),gm2(1))"]] - thetas[["Y(1,gm1(1),gm2(0))"]]
    thetas[["IIE_Mdep"]] <- thetas[["IIE_Mjo(1,)"]] - thetas[["IIE_M1(1,,gm2(0))"]] - thetas[["IIE_M2(1,gm1(1),)"]]
    
    thetas
}
