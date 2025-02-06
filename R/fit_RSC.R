#' Condition an operating model for Giant Red Sea Cucumber
#'
#' Uses the RCM model of OpenMSE to fit an operating model to data
#'
#' @param In A list object of class ('In') that includes the operating model parameters (slot 2, class OM) and data (slot 3, class RCMinput) for conditioning
#' @param sims Integer or vector of integers - the number of simulations for the operating model (e.g. 96) or the specific vector of simulations for the operating model (e.g. 13:48)
#' @param max_F Positive real number - the maximum instantaneous mortality rate in any historical time step
#' @param resample Logical - should parameters be drawn from the variance covariance matrix of the MLE fit (T) or just the MLE parameter values (F)
#' @param silent Logical - should all RCM messages be repressed (ie not read out onto the console)?
#' @examples
#' cond.RSC("In.RSC.7F")
#' @author T. Carruthers
#' @seealso \link{RCM}
#' @export
cond.RSC = function(RCMinput, sims = 12, max_F=0.5, resample = T, silent=T){
  setup()
  if(length(sims) == 1) simy = 1:sims
  if(length(sims) > 1) simy = sims
  OM = SubCpars(RCMinput$OM, simy)


  RCMfit = RCM(OM, RCMinput[[3]], s_selectivity = c("B","B"),
             max_F = max_F, mean_fit = T, condition = "catch", cores = 8,
             drop_nonconv=T, drop_highF=T,resample = T, silent=silent)

  RCMfit@OM@Name = paste0("Giant Red Sea Cucumber QMA ",RCMinput[[1]])
  RCMfit

}


