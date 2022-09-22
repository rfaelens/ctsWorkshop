library(mrgsolve)

simulationModelBase <- "
$PLUGIN autodec nm-vars
$CMT  A0 A1 A2 CAUC
$PARAM MOD=1

$OMEGA @annotated
EKA: 0 : ETA on absorption rate
EVc: 0 : ETA on volume
ECL: 0 : ETA on clearance
$CAPTURE CONC
"

tvSimulationModel <- "
$PK
KA=0.126;
Vc=2700;
CL=34.1;
Vp=774;
Q=0.688;
ALAG_A0 = 0.527;

$DES
dxdt_A0 = -KA*A0;
dxdt_A1 = KA*A0 - CL/Vc * A1 - Q/Vc * A1 + Q/Vp*A2;
dxdt_A2 = Q/Vc*A1 - Q/Vp*A2;
CONC = A1 / Vc * 1000;
dxdt_CAUC = CONC;
"
iivModel <- "
$PK
KA=0.126*exp(EKA);
Vc=2700*exp(EVc);
CL=34.1*exp(ECL);
Vp=774;
Q=0.688;
TLag = 0.527;

$DES
dxdt_A0 = -KA*A0;
dxdt_A1 = KA*A0 - CL/Vc * A1 - Q/Vc * A1 + Q/Vp*A2;
dxdt_A2 = Q/Vc*A1 - Q/Vp*A2;
CONC = A1 / Vc * 1000;
dxdt_CAUC = CONC;
"
iivModelMod <- "
$PK
KA=0.126*exp(EKA);
Vc=2700*exp(EVc);
CL=34.1*exp(ECL)*MOD;
Vp=774;
Q=0.688;
TLag = 0.527;

$DES
dxdt_A0 = -KA*A0;
dxdt_A1 = KA*A0 - CL/Vc * A1 - Q/Vc * A1 + Q/Vp*A2;
dxdt_A2 = Q/Vc*A1 - Q/Vp*A2;
CONC = A1 / Vc * 1000;
dxdt_CAUC = CONC;
"

buildModel <- function(text) {
  mrgsolve::mread_cache(model=digest::sha1(text), code=paste0(simulationModelBase, text))
}

simulate <- function(model, dose, nbr.doses, dosing.interval, observationTimes, omega=NULL, nSub=1L, params=NULL) {
  cat("Simulating...")
  cat("(model=", model$model, ", dose=", dose, ", nbr.doses=", nbr.doses, ", interval=", dosing.interval, ", nsub=", nSub, ")")
  events <- mrgsolve::ev(time=0, amt=dose, addl=nbr.doses-1, ii=dosing.interval)
  if(!is.null(omega)) {
    oNames <- as.list(omat(model))$... %>% rownames()
    omega2 <- omega[oNames, ][, oNames]
    model <- omat(model, omega2)
  }
  if(!is.null(params)) {
    model <- param(model, params)
  }
  result <- mrgsim_df(model, 
                      events=events, 
                      nid=nSub, 
                      tgrid=tgrid(start=observationTimes[1], end=observationTimes[1], add=observationTimes[-1]),
                      obsonly=T) %>%
    mutate(sim.id = ID)
  #result <- rxSolve(model, events=events, omega = omega, nSub=nSub, params=params)
  cat("DONE\n")
  result
}

